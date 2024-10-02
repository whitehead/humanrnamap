import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.mstats import winsorize
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
            "VARNA jar file not found. Please download it from http://varna.lri.fr/ and place it in the root directory of this project.")
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

    directories = [wkdir, fastadir, datdir, ctdir,
                   csvdir, fold_dbndir, final_image, fold_FEdir, barplotdir, jsondir]

    ensure_directories_exist(directories)

    data_dir = os.path.abspath('./data')
    bgdir = data_dir + '/chrom'
    gene_names_bed_path = data_dir + '/genename.bed'
    seq_fasta_dir = data_dir + '/fasta/'

    log_timing_and_memory("setup complete")


    # Get the input variables from the configuration file if they are set, otherwise use user input
    if config.has_option('general', 'AOI'):
        AOI = config.get('general', 'AOI')
    else:
        AOI = input('input a title for your area of interest: ')

    scrubbed_aoi = scrub_filename(AOI)

    # obtain list of possible chromosome inputs from user
    lst = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
           '20', '21', '22', 'X', 'Y', 'M']
    # set active condition as true for the while loop
    active = True
    while active:
        # take in chromosome number or letter from user
        if config.has_option('chromosome', 'user_chr'):
            user_chr = config.get('chromosome', 'user_chr')
        else:
            user_chr = input('input the chromosome number, X, Y, or M (hg38): ')
        # set an if statement if the user inputs a lowercase letter for x, y, m. create list of letters
        if user_chr in ['x', 'y', 'm']:
            # turn lowercase letter into uppercase letter
            user_chr = user_chr.upper()
            # set active to false to break loop
            active = False
        # create an elif statement for an invalid input to reloop
        elif user_chr not in lst:
            print(f"Invalid input {user_chr}. Please enter a valid chromosome number or letter.")
        # create an statement to break loop
        else:
            active = False
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

    # loop for coordinate entries
    working = True
    while working:  # for coordinate set options
        if config.has_option('coordinates', 'coord_opt'):
            coord_opt = config.get('coordinates', 'coord_opt')
        else:
            coord_opt = input('would you like to enter 1 or 2 sets of coordinates: ')
        if coord_opt == '1':
            running = True
            while running:  # option 1 starting coordinates
                if config.has_option('coordinates', 'start_coord1'):
                    start_coord = config.get('coordinates', 'start_coord1')
                else:
                    start_coord = input('input the starting coordinate (hg38): ')
                # create a list of individual strings for each character
                s_nums = [x for x in start_coord]
                # create a while loop to remove commas from user input
                while ',' in s_nums:
                    s_nums.remove(',')
                # join user input string without commas
                s_str = ''.join(s_nums)
                # create an if statement for the user inputting a string that doesnt consist of only numbers
                if (s_str.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    running
                # create elif statement where input is only numbers
                elif (s_str.isnumeric() == True):
                    f_s = int(s_str)
                    if cmax <= f_s:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        running
                    # break loop if input is adequate
                    else:
                        running = False
            executing = True
            while executing:  # option 1 ending coordinates
                # take in coordinate as a string from user
                if config.has_option('coordinates', 'end_coord1'):
                    end_coord = config.get('coordinates', 'end_coord1')
                else:
                    end_coord = input('input the ending coordinate (hg38): ')
                # create a list of individual strings for each character
                e_nums = [x for x in end_coord]
                # remove commas from user input
                while ',' in e_nums:
                    e_nums.remove(',')
                e_str = ''.join(e_nums)
                # create an if statement for the user inputting a string that doesnt consist of only numbers
                if (e_str.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    executing
                # create elif statement where input is only numbers
                elif (e_str.isnumeric() == True):
                    # convert string to integer
                    f_e = int(e_str)
                    # create an if statement if the user input is larger than the maximum chromosome number
                    if cmax < f_e:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        executing
                    elif f_e == f_s:
                        print('this number is equal to the starting coordinate. try again.')
                        # reloop through while to promt user to input another value
                        executing
                    # create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e < f_s:
                        print('the end coordinate must be larger than the start coordiate. try again.')
                        executing
                    # create an else statemnt to break loop when input is adequate
                    else:
                        executing = False
                        working = False

        elif coord_opt == '2':
            playing = True
            while playing: 
                if config.has_option('coordinates', 'start_coord1'):
                    start_coord1 = config.get('coordinates', 'start_coord1')
                else:
                    start_coord1 = input('input the starting coordinate (hg38): ')
                s_nums1 = [x for x in start_coord1]
                # create a while loop to remove commas from user input
                while ',' in s_nums1:
                    s_nums1.remove(',')
                s_str1 = ''.join(s_nums1)
                # create an if statement for the user inputting a string that doesnt consist of only numbers
                if (s_str1.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    playing
                elif (s_str1.isnumeric() == True):
                    # convert string to integer
                    f_s1 = int(s_str1)
                    if cmax <= f_s1:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        playing
                    # break loop if input is adequate
                    else:
                        playing = False
            active = True
            while active:  # option 2 first set of ending coordinates
                if config.has_option('coordinates', 'end_coord1'):
                    end_coord1 = config.get('coordinates', 'end_coord1')
                else:
                    end_coord1 = input('input the ending coordinate (hg38): ')
                e_nums1 = [x for x in end_coord1]
                while ',' in e_nums1:
                    e_nums1.remove(',')
                e_str1 = ''.join(e_nums1)
                if (e_str1.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    active
                elif (e_str1.isnumeric() == True):
                    f_e1 = int(e_str1)
                    if cmax < f_e1:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        active
                    # create an elif statement if user input is equal to staarting coordinate
                    elif f_e1 == f_s1:
                        print('this number is equal to the starting coordinate. try again.')
                        # reloop through while to promt user to input another value
                        active
                    # create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e1 < f_s1:
                        print('the end coordinate must be larger than the start coordinate. try again.')
                        active
                    # create an else statemnt to break loop when input is adequate
                    else:
                        active = False
                        working = False
            alert = True
            while alert:  # option 2 second starting coordinate
                if config.has_option('coordinates', 'start_coord1'):
                    start_coord2 = config.get('coordinates', 'start_coord2')
                else:
                    start_coord2 = input('input the second starting coordinate: ')
                s_nums2 = [x for x in start_coord2]
                while ',' in s_nums2:
                    s_nums2.remove(',')
                s_str2 = ''.join(s_nums2)
                if (s_str2.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    alert
                elif (s_str2.isnumeric() == True):
                    f_s2 = int(s_str2)
                    if cmax <= f_s2:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        alert
                    # break loop if input is adequate
                    else:
                        alert = False
            conducting = True
            while conducting:  # option 2 second set ending coordinate
                if config.has_option('coordinates', 'end_coord2'):
                    end_coord2 = config.get('coordinates', 'end_coord2')
                else:
                    end_coord2 = input('input the second ending coordinate: ')
                e_nums2 = [x for x in end_coord2]
                while ',' in e_nums2:
                    e_nums2.remove(',')
                e_str2 = ''.join(e_nums2)
                if (e_str2.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    # reloop though while to prompt user to enter another value
                    conducting
                elif (e_str2.isnumeric() == True):
                    f_e2 = int(e_str2)
                    if cmax < f_e2:
                        print('this number is too large for the selected chromosome. try again.')
                        # reloop through while to prompt user to input another value
                        conducting
                    # create an elif statement if user input is equal to staarting coordinate
                    elif f_e2 == f_s2:
                        print('this number is equal to the starting coordinate. try again.')
                        # reloop through while to promt user to input another value
                        conducting
                    # create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e2 < f_s2:
                        print('the end coordinate must be larger than the start coordinate. try again.')
                        conducting
                    # create an else statemnt to break loop when input is adequate
                    else:
                        conducting = False
                        working = False
        # create else for invalid inputs
        else:
            print("invalid input. try again.")
            working

    managing = True
    while managing:  # option two stranded options
        # take in user input
        if config.has_option('strand', 'choice'):
            choice = config.get('strand', 'choice')
        else:
            choice = input('enter a "+" for a positive strand and a "-" for a negative strand: ')
        # create an if statement for positive strand
        if (choice == "+"):
            strand = 'pos'
            # break loop for correct input
            managing = False
        # create an elif statement for negative strand
        elif (choice == "-"):
            strand = 'neg'
            # break loop for correct input
            managing = False
        # create an else statement for incorrect input
        else:
            print(f"Invalid input {choice}. Please enter a valid strand. (Either + or -)")
            # reloop through while to prompt user to input another value
            managing
    
    # DONE WITH INPUTS!

    log_timing_and_memory("config complete")

    # read in canonical annotation bed file 
    anno = pd.read_csv(gene_names_bed_path, sep='\t')
    log_timing_and_memory("read annotations")
    gene_list = []
    # create an if condition for gene name if there is only 1 set of coordinates
    if coord_opt == '1':
        # iterate through the data frame to find gene that satisfies all user inputs.
        for i in range(0, len(anno)):
            if (((((anno['start'][i] < f_s) & (anno['end'][i] >= f_s)) | (
                    (anno['start'][i] < f_e) & (anno['end'][i] >= f_e)))) & (
                    (anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                # append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
        # set conditions to display gene(s) to the user
        gene_list = [x for x in gene_list if str(x) != 'nan']      

    # create an if condition for gene name if there are 2 sets of coordinates
    if coord_opt == '2':
        # iterate through the data frame to find gene that satisfies all user inputs.
        for i in range(0, len(anno)):
            if (((((anno['start'][i] < f_s1) & (anno['end'][i] >= f_s1)) | (
                    (anno['start'][i] < f_e1) & (anno['end'][i] >= f_e1)))) & (
                    (anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                # append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
            if (((((anno['start'][i] < f_s2) & (anno['end'][i] >= f_s2)) | (
                    (anno['start'][i] < f_e2) & (anno['end'][i] >= f_e2)))) & (
                    (anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                # append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
        gene_list = [x for x in gene_list if str(x) != 'nan']
        
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

    # create an if when user selects only one set of coordinates
    if coord_opt == '1':

        gene_strand = stream_filter_bedgraph_file(f'Davg_{chrom}_mmRate_{strand}.bedGraph', f_s, f_e)

    # create an elif for when user slects to have two sets of coordinates
    elif coord_opt == '2':
        gene_strand = pd.concat([
            stream_filter_bedgraph_file(f'Davg_{chrom}_mmRate_{strand}.bedGraph', f_s1, f_e1),
            stream_filter_bedgraph_file(f'Davg_{chrom}_mmRate_{strand}.bedGraph', f_s2, f_e2)
        ])

    # read in fasta file for sequence and add ref column to df. NOTE: this script requires the fasta sequence to be written without line breaks.
    abs_seq_file_path = os.path.abspath(f"{seq_fasta_dir}{chrom}.fasta")

    seq_file = open(abs_seq_file_path)

    def read_specific_position(file_path, start_pos, end_pos):
        with open(file_path, 'rb') as file:
            file.readline()  # Skip the first line
            file.seek(start_pos, 1)  # Skip to the start position from the current position (beginning of second line)
            data = file.read(end_pos - start_pos)  # Read up to the end position
            return [char for char in data.decode()]

    read_specific_position_seq = lambda start_pos, end_pos: read_specific_position(abs_seq_file_path, start_pos, end_pos)

    log_timing_and_memory("write fasta")

    # this function removes the top 5% of outliers and scales that output to a max of 1, as typical for DMS data processing
    def winsor_scale(col):
        arr = np.array(col)
        win_arr = winsorize(arr, limits=[0.00, 0.05]).data  # winsorize with limits of 0%, 95%.
        norm_arr = (win_arr - win_arr.min()) / (win_arr.max() - win_arr.min())  # scale to 1
        return norm_arr

    gene_strand.rename(columns={'D' : 'mmRate'}, inplace=True)
    # winsorize & scale raw mismatch rates
    gene_strand['D'] = winsor_scale(gene_strand['mmRate'])

    if coord_opt == '1':
        RNAstructure = pd.DataFrame()
        RNAstructure['pos'] = range(f_s, f_e + 1) #list of coordinates based on user input
        RNAstructure['coord'] = ''
        RNAstructure['ref'] = read_specific_position_seq(f_s - 1, f_e)
        # merge RNAstructure with gene_strand to add DMS values to RNAstructure
        RNAstructure = RNAstructure.merge(gene_strand, on=['pos'], how='left')
    # create if statement when the user input is two coordinates
    if coord_opt == '2':
        pos_list = list(range(f_s1, f_e1 + 1)) + list(range(f_s2, f_e2 + 1))
        RNAstructure = pd.DataFrame()
        RNAstructure['pos'] = pos_list
        RNAstructure['coord'] = ''
        two_seq = read_specific_position_seq(f_s1 - 1,f_e1) + read_specific_position_seq(f_s2 - 1,f_e2)
        RNAstructure['ref'] = two_seq
        # merge RNAstructure with gene_strand to append coverage values to RNAstructure
        RNAstructure = RNAstructure.merge(gene_strand, on=['pos'], how='left')

    # for reversely stranded genes, reverse the order of the bases and replace them with their complements
    if strand == 'neg':
        RNAstructure.sort_values('pos', ascending=False, inplace=True)
        RNAstructure.reset_index(drop=True, inplace=True)  # reverses the order
        for i in range(0, len(RNAstructure)):  # replace with thier complements
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
            if RNAstructure['ref'][i] == 'T':
                RNAstructure['ref'][i] = 'U'

    # fill coord column starting from 1, now that strandedness has been corrected where necessary
    RNAstructure['coord'] = np.arange(1, len(RNAstructure) + 1) 
    RNAstructure.head()

    # create output messages for region length
    region_len = len(RNAstructure)
    if region_len <= 500:
        add_message('Length', 'info', 'Length of region: ' + str(region_len) + ' nt.')
        print('Length of region: ' + str(region_len) + ' nt.')
    else:
        add_message('Length', 'error', 'Input region length: ' + str(region_len) + ' nt. Maximum length: 500 nt. For longer regions, download and execute the code locally.')
        print('Input region length: ' + str(region_len) + ' nt. Maximum length: 500 nt when running in online database. For longer regions, download and execute the code locally.')
    

    # create data frame for plot only including A's and C's
    plot = pd.DataFrame()
    plot = RNAstructure[(RNAstructure['ref'] == 'A') | (RNAstructure['ref'] == 'C')]
    plot.rename(columns={'D' : 'DMS_signal', 'pos' : 'chrom_coord', 'coord' : 'local_coord'}, inplace = True)
    plot.to_csv(csvdir + '/' + scrubbed_aoi + '.csv', index=False)


    # plot using RNAstructure frame with dms signal as y and local_coord as x
    coordinates = plot['local_coord']
    mismatch = plot['DMS_signal']
    plt.figure(figsize=(10, 7))
    plt.bar(coordinates, mismatch)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylim(0, 1)
    plt.xlabel("coordinates", fontsize=14)
    plt.ylabel("mismatch rates", fontsize=14)
    plt.title(AOI, fontsize=18)
    plt.savefig(barplotdir + '/' + scrubbed_aoi + '_plot.svg')


    RNAstructure.reset_index(drop=True, inplace=True)
    RNAstructure['D'] = RNAstructure['D'].fillna(-999) # bases without sufficient coverage and G/Ts should not be considered in calculation

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
    run(f"{fold_bin_path} \"{os.path.join(fastadir, f'{scrubbed_aoi}.fasta')}\" \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" -dms \"{os.path.join(datdir, f'{scrubbed_aoi}.dat')}\"")
    log_timing_and_memory("fold")

    # run command to obtain free folding energy
    run(f"{efn2_bin_path} \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" \"{os.path.join(fold_FEdir, f'{scrubbed_aoi}.txt')}\" -dms \"{os.path.join(datdir, f'{scrubbed_aoi}.dat')}\"")
    log_timing_and_memory("efn2")

    # read the number of lines in the energy file to display the amount of structures to the user
    FE_lines = open(fold_FEdir + '/' + scrubbed_aoi + '.txt')
    struct_num = len(FE_lines.readlines())

    add_message('Structures', 'info', 'This region has ' + str(struct_num) + ' maximum predicted structures. 5 is the maximum structures visible on this site. For up to 20 predictions per region, download and run the code locally.')
    print('Structures', 'info', 'This region has ' + str(struct_num) + ' maximum predicted structures. 5 is the maximum structures visible on the web server. For up to 20 predictions per region, download and run the code locally.')

    write_messages_to_json(filename = jsondir + '/' + scrubbed_aoi + '.json')

    for x in range(1, struct_num + 1): #loops through to visualize all structures
        user_struct = str(x)
        FE_lines2 = open(fold_FEdir + '/' + scrubbed_aoi + '.txt')
        ttl = FE_lines2.readlines()
        y = ttl[x - 1]
        ttl2 = y.split(' ')
        s_num = ttl2[1]
        energy = ' '.join(ttl2[4:])
        # obtain dbn file for specified structure
        # ct2dot converts individual structure in ct to dbn 
        run(f"{ct2dot_bin_path} \"{os.path.join(ctdir, f'{scrubbed_aoi}.ct')}\" {user_struct} \"{os.path.join(fold_dbndir, f'{scrubbed_aoi}_{user_struct}.dbn')}\"")
        log_timing_and_memory("ct2dot")

        # open the new dbn file
        fold_dbn = open(fold_dbndir + '/' + scrubbed_aoi + '_' + user_struct + '.dbn')
        # read third line for dbn notation as list
        # -1 for extra space at the end of line
        struct1 = fold_dbn.readlines()[2][:-1]

        # create a lst1 for bases without DMS data (Gs/Ts, low-coverage bases) and lst 2 for bases with valid data (high-cov As/Cs)
        lst1 = []
        lst2 = []

        # iterate though list of individual letters
        for val in range(len(RNAstructure)):
            if RNAstructure['D'][val] == -999: #from fillna earlier
                lst1.append(val)
            else:
                lst2.append(val)
        # create list low, mid and high for low, mid, and high signal values, for colormap purposes
        low = []
        mid = []
        high = []
        # append each value to their respective lists based on coverage
        for val in lst2:
            if (0.00 <= RNAstructure['D'][val] < 0.25): 
                low.append(val)
            elif (0.25 <= RNAstructure['D'][val] < 0.5):
                mid.append(val)
            elif (0.5 <= RNAstructure['D'][val] <= 1.0):
                high.append(val)

        # change directory for varna imaging
        os.chdir(wkdir)

        # get image for fold dbn
        from varnaapi import VARNA
        import varnaapi
        varnaapi.set_VARNA(varna_path)
        v = VARNA(structure=struct1, sequence=seq)
        style1 = varnaapi.param.BasesStyle(fill="#FFFFFF") # white; bases without data
        style2 = varnaapi.param.BasesStyle(fill="#00FFFF") # blue; low signal
        style3 = varnaapi.param.BasesStyle(fill="#FFFF00") # yellow; moderate signal
        style4 = varnaapi.param.BasesStyle(fill="#FF0000") # red; high signal
        v.add_bases_style(style1, lst1)
        v.add_bases_style(style2, low)
        v.add_bases_style(style3, mid)
        v.add_bases_style(style4, high)

        v.set_title(AOI + ', structure ' + s_num + ', ' + energy + ' kcal/mol', color='#000000', size=12)

        preliminary_figure_path = final_image + '/' + scrubbed_aoi + '_fold_' + user_struct + '.svg'
        v.savefig(preliminary_figure_path)

        fig = sg.fromfile(preliminary_figure_path)
        fig.set_size(('2000', '2000'))

        final_svg_path = final_image + '/' + scrubbed_aoi + '_fold_final_' + user_struct + '.svg'
        fig.save(final_svg_path)

        # delete the preliminary figure
        os.remove(preliminary_figure_path)

        # Trim svg (inplace)
        trim_svgs.trim_svg_file(final_svg_path)

        # Append legend (inplace)
        append_svgs.append_svg_legend(final_svg_path, LEGEND_SVG_PATH)

        # break to avoid infinite loop
        finishing = False

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
