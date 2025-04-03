import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.mstats import winsorize
from sklearn.metrics import roc_auc_score
import seaborn as sns
import svgutils.transform as sg
import sys
import shutil
import configparser
import argparse
import subprocess
import csv
import time
import psutil
import json

# Local libraries
import append_svgs
import trim_svgs

# Initialize variables to store the start time and previous call time
start_time = time.time()
previous_call_time = start_time

LEGEND_SVG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "data", "colorbar.svg"))

def log_timing_and_memory(msg):
    global previous_call_time
    current_time = time.time()
    elapsed_time = current_time - previous_call_time
    process_id = os.getpid()
    process = psutil.Process(process_id)
    memory_info = process.memory_info()
    memory_usage_mb = memory_info.rss / (1024 * 1024)
    print(f"{msg} - Elapsed time since the previous call: {elapsed_time:.2f} seconds")
    print(f"{msg} - Current memory usage: {memory_usage_mb:.2f} MB")
    previous_call_time = current_time


def main():
    log_timing_and_memory("first")
    # Commandline argument parsing
    parser = argparse.ArgumentParser(description="Specify the configuration file")
    parser.add_argument('--config', '-c', type=str, required=True, help='The configuration file to use')
    parser.add_argument('--outdir', '-o', type=str, required=False, help='The output directory to use')
    args = parser.parse_args()

    # Read the TOML configuration file
    config = configparser.ConfigParser()
    config.read(args.config)

    # Ensure the binaries we will need are installed and in the path.
    fold_bin_path = abs_binary_path("Fold")
    efn2_bin_path = abs_binary_path("efn2")
    ct2dot_bin_path = abs_binary_path("ct2dot")

    varna_path = os.path.abspath('./bin/VARNAv3-93.jar')
    if not os.path.exists(varna_path):
        print(
            "ERR: VARNA jar file not found. Please download it from http://varna.lri.fr/ and place it in the root directory of this project.")
        sys.exit(1)

    # create variables for the paths to different files

    # arg outdir or
    if args.outdir:
        wkdir = os.path.abspath(args.outdir)
    else:
        wkdir = os.path.abspath('./out/')
    fastadir = wkdir + '/fasta'
    datdir = wkdir + '/fold_constraints'
    ctdir = wkdir + '/ct'
    csvdir = wkdir + '/csv'
    fold_dbndir = wkdir + '/fold_dbn'
    final_image = wkdir + '/final_image'
    fold_FEdir = wkdir + '/fold_FE'
    barplotdir = wkdir + '/barplot'
    jsondir = wkdir + '/json'
    aucdir = wkdir + '/auc'

    directories = [wkdir, fastadir, datdir, ctdir,
                   csvdir, fold_dbndir, final_image, fold_FEdir, barplotdir, jsondir, aucdir]

    ensure_directories_exist(directories)

    data_dir = os.path.abspath('./data')
    bgdir = data_dir + '/chrom'
    gene_names_bed_path = data_dir + '/genename.bed'
    seq_fasta_dir = data_dir + '/fasta/'

    log_timing_and_memory("setup complete")


    # Get the input variables from the configuration file if they are set, otherwise use user input
    if config.has_option('general', 'AOI'):
        AOI = config.get('general', 'AOI')
        if AOI == '':
            print('ERR: Please enter an AOI in the config file. Exiting.')
            sys.exit()
    else:
        print('ERR: Please enter an AOI in the config file. Exiting.')
        sys.exit()

    scrubbed_aoi = scrub_filename(AOI)

    # obtain list of possible chromosome inputs from user
    chr_lst = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
           '20', '21', '22', 'X', 'Y', 'M']

    # take in chromosome number or letter from user
    if config.has_option('chromosome', 'user_chr'):
        user_chr = config.get('chromosome', 'user_chr')
    else:
        print('ERR: Please enter a chromosome in the config file. Exiting.')
        sys.exit()
    # set an if statement if the user inputs a lowercase letter for x, y, m. create list of letters
    if user_chr in ['x', 'y', 'm']:
        # turn lowercase letter into uppercase letter
        user_chr = user_chr.upper()
    # create an elif statement for an invalid input 
    elif user_chr not in chr_lst:
        print(f"ERR: Invalid chromosome input: {user_chr}. Please enter a valid chromosome number or letter. Exiting.")
        sys.exit()

    # concantenate the user input with chr to resemble the chrx format in the file with x being the selected chromosme
    chrom = 'chr' + user_chr
    if user_chr == '1':
        cmax = 248956422
    elif user_chr == '2':
        cmax = 242193529
    elif user_chr == '3':
        cmax = 198295559
    elif user_chr == '4':
        cmax = 190214555
    elif user_chr == '5':
        cmax = 181538259
    elif user_chr == '6':
        cmax = 170805979
    elif user_chr == '7':
        cmax = 159345973
    elif user_chr == '8':
        cmax = 145138636
    elif user_chr == '9':
        cmax = 138394717
    elif user_chr == '10':
        cmax = 133797422
    elif user_chr == '11':
        cmax = 135086622
    elif user_chr == '12':
        cmax = 133275309
    elif user_chr == '13':
        cmax = 114364328
    elif user_chr == '14':
        cmax = 107043718
    elif user_chr == '15':
        cmax = 101991189
    elif user_chr == '16':
        cmax = 90338345
    elif user_chr == '17':
        cmax = 83257441
    elif user_chr == '18':
        cmax = 80373285
    elif user_chr == '19':
        cmax = 58617616
    elif user_chr == '20':
        cmax = 64444167
    elif user_chr == '21':
        cmax = 46709983
    elif user_chr == '22':
        cmax = 50818468
    elif user_chr == 'X':
        cmax = 156040895
    elif user_chr == 'Y':
        cmax = 57227415
    elif user_chr == 'M':
        cmax = 16569

    # get argument for number of start/end coords ("coord_opt") from config file
    if config.has_option('coordinates', 'coord_opt'):
        coord_opt = config.get('coordinates', 'coord_opt')
    else:
        print('ERR: Please enter the number of coordinate sets to enter in the coord_opt line in the config file. Exiting.')
        sys.exit()
    if not coord_opt.isnumeric():
        print('ERR: Invalid input for coord_opt. Please enter an integer. Exiting.')
        sys.exit()
    coord_opt = int(coord_opt)
    if coord_opt == 0:
        print('ERR: Invalid input for coord_opt. Please enter an integer > 0. Exiting.')
        sys.exit()
    
    # now get the actual coordinate arguments, appending values to these lists for each coordinate set
    start_coords = []
    end_coords = []
    
    for i in range(1, coord_opt + 1): # loop through each set of start/end coords
        if config.has_option('coordinates', 'start_coord' + str(i)):
            curr_start_coord = config.get('coordinates', 'start_coord' + str(i))
        else:
            print('ERR: Please enter the same number of start/end coordinates as expected by the selected number of coordinate sets. Exiting.')
            sys.exit()
        # create a list of individual strings for each character
        s_nums = [x for x in curr_start_coord]
        # remove commas from user input
        while ',' in s_nums:
            s_nums.remove(',')
        # join user input string without commas
        s_str = ''.join(s_nums)
        # error out when input does not contain only numbers (after comma removal). Also excludes decimals
        if (s_str.isnumeric() == False):
            print('ERR: Coordinate inputs must consist of only numbers. Exiting.')
            sys.exit()
        elif (s_str.isnumeric() == True):
            int_curr_start = int(s_str)
            if int_curr_start == 0:
                print('ERR: Coordinate inputs must be integers greater than 0. Exiting.')
                sys.exit()
            if cmax <= int_curr_start:
                print('ERR: Coordinate input is too large for the selected chromosome. Exiting.')
                sys.exit()
            start_coords.append(int_curr_start)
            if max(start_coords) != int_curr_start:
                print('ERR: Coordinate sets must be entered in ascending order (even for negatively-stranded genes). Please re-enter in ascending order. Exiting.')
                sys.exit()

        if config.has_option('coordinates', 'end_coord' + str(i)):
            curr_end_coord = config.get('coordinates', 'end_coord' + str(i))
        else:
            print('ERR: Please enter the same number of start/end coordinates as expected by the selected number of coordinate sets. Exiting.')
            sys.exit()
        # create a list of individual strings for each character
        e_nums = [x for x in curr_end_coord]
        # remove commas from user input
        while ',' in e_nums:
            e_nums.remove(',')
        e_str = ''.join(e_nums)
        # error out when input does not contain only numbers (after comma removal). Also excludes decimals
        if (e_str.isnumeric() == False):
            print('ERR: Coordinate inputs must consist of only numbers. Exiting.')
            sys.exit()
        elif (e_str.isnumeric() == True):
            int_curr_end = int(e_str)
            # create an if statement if the user input is larger than the maximum chromosome number
            if int_curr_end == 0:
                print('ERR: Coordinate inputs must be integers greater than 0. Exiting.')
                sys.exit()
            if cmax < int_curr_end:
                print('ERR: Coordinate input is too large for the selected chromosome. Exiting.')
                sys.exit()
            elif int_curr_end <= int_curr_start:
                print('ERR: End coordinate is less than or equal to the start coordinate. Exiting.')
                sys.exit()
            end_coords.append(int_curr_end)
            if max(end_coords) != int_curr_end:
                print('ERR: Coordinate sets must be entered in ascending order. Exiting.')
                sys.exit()
    
    # now have lists of start and end coords
    # check for overlapping regions
    all_coord_list = []
    for i in range(len(start_coords)):
        all_coord_list.extend(list(range(start_coords[i], end_coords[i] + 1)))
    if len(all_coord_list) != len(set(all_coord_list)): # checks for duplicate coordinates from overlapping coord sets
        print('ERR: Coordinate sets are overlapping. Please re-enter the coordinate sets with non-overlapping coordinates. Exiting.')
        sys.exit()    

    # take in user input for strandedness
    if config.has_option('strand', 'choice'):
        choice = config.get('strand', 'choice')
    else:
        print('ERR: Please enter a strandedness option (+ or -) in the config file. Exiting.')
        sys.exit()
    
    # define strand variables for appropriate data file lookup later
    if (choice == "+"):
        strand = 'pos'
    elif (choice == "-"):
        strand = 'neg'

    # exit if incorrect strand input
    else:
        print(f"ERR: Invalid input {choice}. Please enter a valid strand (either + or -). Exiting.")
        sys.exit()
    
    # DONE WITH INPUTS!

    log_timing_and_memory("config complete")

    # read in canonical annotation bed file 
    anno = pd.read_csv(gene_names_bed_path, sep='\t')
    log_timing_and_memory("read annotations")
    gene_list = []
    
    for i in range(len(start_coords)):
        # for each set of start/end coords, get all intersecting regions from the annotation and get a list of gene names from those regions  
        anno_current = anno[((((anno['start'] < start_coords[i]) & (anno['end'] >= start_coords[i])) | # where the region start coord falls within an annot region (exon) OR
                    ((anno['start'] < end_coords[i]) & (anno['end'] >= end_coords[i]))) & # where the region end coord falls within an exon
                    (anno['chrom'] == chrom))] # AND  where the chromosome matches. Note that this is not a strand-specific selection.
        gene_strand_vals = anno_current['gene'] + '(' + anno_current['strand'] + ')'
        gene_list.extend(list(gene_strand_vals.unique()))
    gene_list = [x for x in gene_list if str(x) != 'nan'] # remove nan values; uncommon      
    gene_list = list(dict.fromkeys(gene_list)) # removes any duplicates
    if len(gene_list) == 1:
        add_message('Annotation', 'info', 'Based on canonical annotations, the following gene is in your area of interest: ' + gene_list[0])
        print('Based on canonical annotations, the following gene is in your area of interest: ' + gene_list[0])
    elif len(gene_list) > 1:
        total_genes = ', '.join(gene_list)
        add_message('Annotation', 'info', 'Based on canonical annotations, the following genes are in your area of interest: ' + total_genes)
        print('Based on canonical annotations, the following genes are in your area of interest: ' + total_genes)
    elif len(gene_list) == 0:
        add_message('Annotation', 'warn', 'Based on canonical annotations, there are not any genes in your area of interest.')

    # change directory to folder with bedGraph files containing DMS mismatch rates
    os.chdir(bgdir)

    gene_strand = pd.DataFrame(columns=['pos', 'D'])
    gene_strand = gene_strand.astype({'pos' : int, 'D' : float})
    for i in range(len(start_coords)):
        gene_strand_current = stream_filter_bedgraph_file(f'Davg_{chrom}_mmRate_{strand}.bedGraph', start_coords[i], end_coords[i])
        if len(gene_strand_current) == 0:
            continue # prevents error in next line changing datatype in no-coverage regions
        gene_strand_current = gene_strand_current.astype({'pos' : int, 'D' : float})
        gene_strand = pd.concat([gene_strand, gene_strand_current], ignore_index = True)
    

    # read in fasta file for sequence and add ref column to df. NOTE: this script requires the fasta sequence to be written without line breaks.
    abs_seq_file_path = os.path.abspath(f"{seq_fasta_dir}{chrom}.fasta")

    seq_file = open(abs_seq_file_path)

    def read_specific_position(file_path, start_pos, end_pos):
        with open(file_path, 'rb') as file:
            file.readline()  # Skip the first line
            file.seek(start_pos, 1)  # Skip to the start position from the current position (beginning of second line)
            data = file.read(end_pos - start_pos)  # Read up to the end position
            return [char for char in data.decode()]

    log_timing_and_memory("write fasta")
    
    read_specific_position_seq = lambda start_pos, end_pos: read_specific_position(abs_seq_file_path, start_pos, end_pos)
    seq_list = []
    for i in range(len(start_coords)):
        seq_list.extend(read_specific_position_seq(start_coords[i] - 1, end_coords[i]))

    # this function removes the top 5% of outliers and scales that output to a max of 1, as typical for DMS data processing
    def winsor_scale(col):
        arr = np.array(col)
        win_arr = winsorize(arr, limits=[0.00, 0.05]).data  # winsorize with limits of 0%, 95%.
        norm_arr = (win_arr - win_arr.min()) / (win_arr.max() - win_arr.min())  # scale to 1
        return norm_arr

    gene_strand.rename(columns={'D' : 'mmRate'}, inplace=True)
    # winsorize & scale raw mismatch rates
    if len(gene_strand) > 0: #skip winsorization for empty dfs (no-cov regions), would throw an error
        gene_strand['D'] = winsor_scale(gene_strand['mmRate'])
        
    elif len(gene_strand) == 0:
        gene_strand['D'] = np.nan # fill with NaN, will be replaced with -999 later (no-constraint value in Fold command)
        add_message('Structures', 'warn', 'NO SIGNAL COVERAGE IN THIS AREA. STRUCTURE PREDICTION WILL BE BASED ON SEQUENCE ALONE.')
        print('NO SIGNAL COVERAGE IN THIS AREA. STRUCTURE PREDICTION WILL BE BASED ON SEQUENCE ALONE.')
    
    if len(gene_strand) < 10:
        add_message('Structures', 'warn', 'Too few A/C datapoints for AUC calculation. Recommended >= 10. Region contains only ' + str(len(gene_strand)) + '.')
        print('Too few A/C datapoints for AUC calculation. Recommended >= 10. Region contains only ' + str(len(gene_strand)) + '.')
    
    pos_list = []
    for i in range(len(start_coords)):
        pos_list.extend(list(range(start_coords[i], end_coords[i] + 1)))
        
    RNAstructure = pd.DataFrame()
    RNAstructure['pos'] = pos_list #list of coordinates based on user input
    RNAstructure['coord'] = ''
    RNAstructure['ref'] = seq_list
    # merge RNAstructure with gene_strand to add DMS values to RNAstructure
    RNAstructure = RNAstructure.merge(gene_strand, on=['pos'], how='left')
    
    # for reversely stranded genes, reverse the order of the bases and replace them with their complements
    if strand == 'neg':
        RNAstructure.sort_values('pos', ascending=False, inplace=True)
        RNAstructure.reset_index(drop=True, inplace=True)  # reverses the order and resets index in this order
        for i in range(0, len(RNAstructure)):  # replace with their complements
            if RNAstructure.loc[i, 'ref'] == "A":
                RNAstructure.loc[i, 'ref'] = 'U'
            elif RNAstructure.loc[i, 'ref'] == 'T':
                RNAstructure.loc[i, 'ref'] = 'A'
            elif RNAstructure.loc[i, 'ref'] == 'G':
                RNAstructure.loc[i, 'ref'] = 'C'
            elif RNAstructure.loc[i, 'ref'] == 'C':
                RNAstructure.loc[i, 'ref'] = 'G'
    elif strand == 'pos':  # for the pos strand, the order and complement of the sequence is fine, but T must be replaced with U
        RNAstructure.reset_index(inplace=True, drop=True)
        for i in range(0, len(RNAstructure)):
            if RNAstructure.loc[i, 'ref'] == 'T':
                RNAstructure.loc[i, 'ref'] = 'U'

    # fill coord column starting from 1, now that strandedness has been corrected where necessary
    RNAstructure['coord'] = np.arange(1, len(RNAstructure) + 1) 

    # create output messages for region length
    region_len = len(RNAstructure)
    if region_len <= 500:
        add_message('Length', 'info', 'Length of region: ' + str(region_len) + ' nt.')
        print('Length of region: ' + str(region_len) + ' nt.')
    else:
        add_message('Length', 'error', 'Input region length: ' + str(region_len) + ' nt. Maximum length: 500 nt. For longer regions, download and execute the code locally.')
        print('Input region length: ' + str(region_len) + ' nt. Maximum length: 500 nt when running in online database. For longer regions, download and execute the code locally.')
    
    if len(gene_strand) > 0:
    # create data frame for plot only including A's and C's, only when there is any valid signal in the defined region
        plot = pd.DataFrame()
        plot = RNAstructure[(RNAstructure['ref'] == 'A') | (RNAstructure['ref'] == 'C')]
        plot2 = plot.rename(columns={'D' : 'DMS_signal', 'pos' : 'chrom_coord', 'coord' : 'local_coord'})
        plot2.to_csv(csvdir + '/' + scrubbed_aoi + '.csv', index=False)


        # plot using RNAstructure frame with dms signal as y and local_coord as x
        coordinates = plot2['local_coord']
        mismatch = plot2['DMS_signal']
        plt.figure(figsize=(10, 7))
        plt.bar(coordinates, mismatch)
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylim(0, 1)
        plt.xlabel("coordinates", fontsize=14)
        plt.ylabel("DMS signal", fontsize=14)
        plt.title(AOI, fontsize=18)
        plt.savefig(barplotdir + '/' + scrubbed_aoi + '_plot.svg')


    RNAstructure.reset_index(drop=True, inplace=True)
    RNAstructure['D'] = RNAstructure['D'].fillna(-999) # bases without sufficient coverage and G/Ts will be ignored in Fold constraints using -999 designation

    seq = ''.join(RNAstructure.ref.tolist())  # creates a string with all bases included in the region

    AC_count = seq.count('A') + seq.count('C')
    AC_cov = round(len(gene_strand)/AC_count, 3) * 100

    # set output messages for coverage 
    if AC_cov >= 70:
        add_message('Coverage', 'info', str(AC_cov) + "% of the region's A/C bases included in the filtered DMS dataset. Ideally, this number should be as high as possible for an experimenally-accurate prediction, and over 70% is recommended.")
        print(str(AC_cov) + "% of the region's A/C bases included in the filtered DMS dataset. Ideally, this number should be as high as possible for an experimenally-accurate prediction, and over 70% is recommended.")

    else:
        add_message('Coverage', 'warn', str(AC_cov) + "% of the region's A/C bases included in the filtered DMS dataset. This is lower than the recommended 70%. While structure prediction is still possible even without any coverage, the structures will be based on the sequence alone, rather than based on experimental constraints.")
        print(str(AC_cov) + "% of the region's A/C bases included in the filtered DMS dataset. This is lower than the recommended 70%. While structure prediction is still possible even without any coverage, the structures will be based on the sequence alone, rather than based on experimental constraints.")

    os.chdir(fastadir)

    # writes the temp fasta file containing the header and sequence
    text_file = open(scrubbed_aoi + '.fasta', 'w')
    text_file.write('>' + scrubbed_aoi + '\n')
    text_file.write(seq)
    text_file.close()

    # change to the .dat directory
    os.chdir(datdir)
    RNAstructure.to_csv(scrubbed_aoi + '.dat', sep='\t', header=False, index=False, columns=['coord', 'D'])
    log_timing_and_memory("write dat")

    # attach sample name to all files
    # fold command line
    run(f"{fold_bin_path} \"{os.path.join(fastadir, f'{scrubbed_aoi}.fasta')}\" \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" --SHAPE \"{os.path.join(datdir, f'{scrubbed_aoi}.dat')}\"")
    log_timing_and_memory("fold")

    # run command to obtain free folding energy
    run(f"{efn2_bin_path} \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" \"{os.path.join(fold_FEdir, f'{scrubbed_aoi}.txt')}\" -sh \"{os.path.join(datdir, f'{scrubbed_aoi}.dat')}\"")
    log_timing_and_memory("efn2")

    # read the number of lines in the energy file to display the amount of structures to the user
    FE_lines = open(fold_FEdir + '/' + scrubbed_aoi + '.txt')
    struct_num = len(FE_lines.readlines())

    add_message('Structures', 'info', 'This region has ' + str(struct_num) + ' maximum predicted structures. 5 is the maximum structures visible on this site. For up to 20 predictions per region, download and run the code locally.')
    print('This region has ' + str(struct_num) + ' predicted structures. 5 is the maximum structures visible on the web server. For up to 20 predictions per region, download and run the code locally.')

    write_messages_to_json(filename = jsondir + '/' + scrubbed_aoi + '.json')
    
    AUC_list = []

    for x in range(1, struct_num + 1): # loops through to convert ct to dbn for each structure, and then calculate the respective AUC value. this structure order is temporary and will be resorted by AUC value in the next loop
        struct_ct_num = str(x)
        temp_struct_num = 'temp' + str(x)
        # obtain dbn file for specified structure
        # ct2dot converts individual structure in ct to dbn 
        run(f"{ct2dot_bin_path} \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" {struct_ct_num} \"{os.path.join(fold_dbndir, f'{scrubbed_aoi}_{temp_struct_num}.dbn')}\"")
        log_timing_and_memory("ct2dot")

        # open the new dbn file
        fold_dbn_file = open(fold_dbndir + '/' + scrubbed_aoi + '_' + temp_struct_num + '.dbn')
        # read third line for dbn notation
        # -1 for extra space at the end of line
        fold_dbn = fold_dbn_file.readlines()[2][:-1]
        
        # get dbn binary for AUC calc
        dbn_binary_string = fold_dbn.replace('(', '0').replace(')', '0').replace('.', '1') 
        dbn_binary_list = [int(x) for x in list(dbn_binary_string)] 
        
        RNAstructure['dbn_binary'] = dbn_binary_list
        
        RNAstructure_dataonly = RNAstructure.drop(index=RNAstructure[RNAstructure['D'] == -999].index) # only valid datapoints, with their corresponding dbn prediction, used in AUC calc
        if len(RNAstructure_dataonly) >= 10:
            AUC = round(roc_auc_score(RNAstructure_dataonly['dbn_binary'].astype('bool').to_list(),
                                RNAstructure_dataonly['D'].to_list()), 3)
            AUC_list.append(AUC)
        
    for x in range(1, struct_num + 1): #loops through to visualize all structures
        
        AUC = AUC_list[x-1]
        
        if len(AUC_list) > 0:
            AUC_struct_num = str(sorted(AUC_list, reverse=True).index(AUC) + 1)
        
        else: 
            AUC_struct_num = str(x)
        
        old_struct_num = 'temp' + str(x)
        
        os.system('mv ' + fold_dbndir + '/' + scrubbed_aoi + '_' + old_struct_num + '.dbn ' + fold_dbndir + '/' + scrubbed_aoi + '_' + AUC_struct_num + '.dbn')
        
        FE_file = open(fold_FEdir + '/' + scrubbed_aoi + '.txt')
        FE_readlines = FE_file.readlines()
        current_FE_line = FE_readlines[x - 1]
        current_FE_line_list = current_FE_line.split(' ')
        energy = ' '.join(current_FE_line_list[6:])

        # open the dbn file from previous loop
        fold_dbn_file = open(fold_dbndir + '/' + scrubbed_aoi + '_' + AUC_struct_num + '.dbn')
        # read third line for dbn notation
        # -1 for extra space at the end of line
        fold_dbn = fold_dbn_file.readlines()[2][:-1]


        # create a lst1 for bases without DMS data (Gs/Ts, low-coverage bases) and lst 2 for bases with valid data (high-cov As/Cs), for colormap purposes
        lst1 = []
        lst2 = []

        # iterate though list of individual letters
        for val in range(len(RNAstructure)):
            if RNAstructure['D'][val] == -999: #from fillna earlier
                lst1.append(val + 1) # read by VARNA starting at pos 1
            else:
                lst2.append(val) # used to index df in a moment, so no adding 1 yet 
        # create list low, mid and high for low, mid, and high signal values, for colormap purposes
        low = []
        mid = []
        high = []
        # append each value to their respective lists based on coverage
        for val in lst2:
            if (0.00 <= RNAstructure['D'][val] < 0.25): 
                low.append(val + 1) # read by VARNA starting at pos 1
            elif (0.25 <= RNAstructure['D'][val] < 0.5):
                mid.append(val + 1)
            elif (0.5 <= RNAstructure['D'][val] <= 1.0):
                high.append(val + 1)

        # change directory for varna imaging
        os.chdir(wkdir)

        # get image for fold dbn
        from varnaapi import VARNA
        import varnaapi
        varnaapi.set_VARNA(varna_path)
        v = VARNA(structure=fold_dbn, sequence=seq)
        style1 = varnaapi.param.BasesStyle(fill="#FFFFFF") # white; bases without data due to being G/U or due to low coverage
        style2 = varnaapi.param.BasesStyle(fill="#00FFFF") # blue; low signal
        style3 = varnaapi.param.BasesStyle(fill="#FFFF00") # yellow; moderate signal
        style4 = varnaapi.param.BasesStyle(fill="#FF0000") # red; high signal
        v.add_bases_style(style1, lst1)
        v.add_bases_style(style2, low)
        v.add_bases_style(style3, mid)
        v.add_bases_style(style4, high)
        
        v.set_title(AOI + ', ' + energy + ' kcal/mol, AUC: ' + str(AUC), color='#000000', size=12) 

        preliminary_figure_path = final_image + '/' + scrubbed_aoi + '_fold_' + AUC_struct_num + '.svg'
        v.savefig(preliminary_figure_path)

        fig = sg.fromfile(preliminary_figure_path)
        fig.set_size(('2000', '2000'))

        final_svg_path = final_image + '/' + scrubbed_aoi + '_fold_final_' + AUC_struct_num + '.svg'
        fig.save(final_svg_path)

        # delete the preliminary figure
        os.remove(preliminary_figure_path)

        # Trim svg (inplace)
        trim_svgs.trim_svg_file(final_svg_path)

        # Append legend (inplace)
        append_svgs.append_svg_legend(final_svg_path, LEGEND_SVG_PATH)
    
    os.chdir(aucdir)   
    auc_file = open(scrubbed_aoi + '_AUCs.txt', 'w')
    if len(AUC_list) > 0:
        for i in range(len(AUC_list)):
            auc_file.write('Structure ' + str(i + 1) +': ' + str(sorted(AUC_list, reverse=True)[i]) + '\n')
    else:
        auc_file.write('AUC ROCs could not be calculated, likely due to lack of coverage in defined region.')
    auc_file.close()

def stream_filter_bedgraph_file(filename, start, end):
    filtered_data = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            pos = int(row[2])
            if start <= pos <= end:
                filtered_data.append({'pos': pos, 'D': float(row[3])})
    return pd.DataFrame(filtered_data)

def abs_binary_path(name):
    which_path = shutil.which(name)
    if which_path is None:
        print(
            f"Error: The {name} binary was not found. Please ensure it is installed and available in the system's PATH.")
        sys.exit(1)

    return os.path.abspath(which_path)


def ensure_directories_exist(directories):
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)


# Run method (system run plus exit error catching)
def run(cmd):
    try:
        # replace the cwd in this command with blank.
        short_cmd = cmd.replace(os.getcwd() + "/", "")
        print(f"Running command: {short_cmd}\n")
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        sys.exit(1)

def scrub_filename(filename):
    filename = re.sub(r'\W+', '_', filename)  # Replace one or more non-alnum characters with a single underscore
    return filename

# Global output message array
messages = []

def add_message(category, level, message):
    global messages
    messages.append({
        "category": category,
        "level": level,
        "message": message
    })

def write_messages_to_json(filename="messages.json"):
    global messages
    with open(filename, 'w') as f:
        json.dump(messages, f, indent=2)

# Run if called directly.
if __name__ == "__main__":
    main()
