import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.mstats import winsorize
import seaborn as sns
import svgutils.transform as sg
import sys
import shutil

def main():

    # Ensure the binaries we will need are installed and in the path.
    fold_bin_path = abs_binary_path("Fold")
    efn2_bin_path = abs_binary_path("efn2")
    ct2dot_bin_path = abs_binary_path("ct2dot")

    #create variables for the paths to different files 
    wkdir = os.path.abspath('./out/')
    fastadir = wkdir + '/fasta'
    datdir = wkdir + '/dat'
    pfsdir = wkdir + '/pfs'
    ctdir = wkdir + '/ct'
    cmapdir = wkdir + '/cmap'
    fold_dbndir = wkdir + '/fold_dbn'
    final_image = wkdir + '/final_image'
    fold_FEdir = wkdir + '/fold_FE'
    barplotdir = wkdir + '/barplot'


    directories = [ wkdir, fastadir, datdir, pfsdir, ctdir,
        cmapdir, fold_dbndir, final_image, fold_FEdir, barplotdir ]

    ensure_directories_exist(directories)


    data_dir = os.path.abspath('./data')
    bgdir = data_dir + '/chrom'
    #gene_names_bed_path = wrkdir + '/bed/genename.bed'
    gene_names_bed_path = data_dir + '/genename.bed'
    #seq_fasta_dir = '/lab/jain_imaging/genomes/hg40.primary_STAR/fasta_stripped/'
    seq_fasta_dir = data_dir + '/fasta/'

    #Input variables and coordinates. 
    #Ideally we'd want the option for at least 2 start/end coords so that they can look at structures over intron/exon junctions.
    #obtain string of title from user
    AOI = input('input a title for your area of interest: ')

    #obtain list of possible chromosome inputs from user
    lst = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
    #set active condition as true for the while loop
    active = True
    while active:
        #take in chromosome number or letter from user
        user_chr = input('input the chromosome number, X, Y, or M (hg38): ')
        #set an if statement if the user inputs a lowercase letter for x, y, m. create list of letters
        if user_chr in ['x', 'y', 'm']:
            #turn lowercase letter into uppercase letter
            user_chr = user_chr.upper()
            #set active to false to break loop
            active = False
        #create an elif statement for an invalid input to reloop
        elif user_chr not in lst:
            print('invalid input.')
            active
    #create an statement to break loop
        else:
            active = False
    #concantenate the user input with chr to resemble the chrx format in the file with x being the selected chromosme
    chrom = 'chr' + user_chr
    #create an if statement with several elif statements to assign the length of the cromosome as the max value
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
        cmax =  80373285
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

    #set working equal to true for while loop
    working = True
    while working: #for coordinate set options
        #take in user input for coordinate set option
        coord_opt = input('would you like to enter 1 or 2 sets of coordinates: ')
        #create an if statement for 1 set of coordinate 
        if coord_opt == '1':
            #set running equal to true for while loop
            running = True
            while running: #option 1 starting coordinates
                #take in coordinate as a string from user
                start_coord = input('input the starting coordinate (hg38): ')
                #create a list of individual strings for each character
                s_nums = [x for x in start_coord]
                #create a while loop to remove commas from user input
                while ',' in s_nums:
                    s_nums.remove(',')
                #join user input string without commas
                s_str = ''.join(s_nums)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (s_str.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    running
                #create elif statement where input is only numbers 
                elif (s_str.isnumeric() == True):
                    #convert string to integer
                    f_s = int(s_str)
                    #create an if statement if the user input is larger or equal to the maximum chromosome number
                    if cmax <= f_s:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        running
                    #break loop if input is adequate
                    else:
                        running = False
            #set executing equal to true for while loop
            executing = True
            while executing: #option 1 ending coordinates
                #take in coordinate as a string from user
                end_coord = input('input the ending coordinate (hg38): ')
                #create a list of individual strings for each character
                e_nums = [x for x in end_coord]
                #create a while loop to remove commas from user input
                while ',' in e_nums:
                    e_nums.remove(',')
                #join user input string without commas
                e_str = ''.join(e_nums)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (e_str.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    executing
                #create elif statement where input is only numbers
                elif (e_str.isnumeric() == True):
                    #convert string to integer
                    f_e = int(e_str)
                    #create an if statement if the user input is larger than the maximum chromosome number
                    if cmax < f_e:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        executing
                    #create an elif statement if user input is equal to staarting coordinate
                    elif f_e == f_s:
                        print('this number is equal to the starting coordinate. try again.')
                        #reloop through while to promt user to input another value
                        executing
                    #create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e < f_s:
                        print('this number is smaller than the starting coordinate. try again.')
                        executing
                    #create an else statemnt to break loop when input is adequate
                    else:
                        executing = False
                        working = False

        #create an elif statement if the user wants to input two sets of coordinates
        elif coord_opt == '2':
            playing = True
            while playing: #option 2 first set of starting coordinates
                #take in coordinate as a string from user
                start_coord1 = input('input the first starting coordinate (hg38): ')
                #create a list of individual strings for each character
                s_nums1 = [x for x in start_coord1]
                #create a while loop to remove commas from user input
                while ',' in s_nums1:
                    s_nums1.remove(',')
                #join user input string without commas
                s_str1 = ''.join(s_nums1)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (s_str1.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    playing
                #create elif statement where input is only numbers 
                elif (s_str1.isnumeric() == True):
                    #convert string to integer
                    f_s1 = int(s_str1)
                    #create an if statement if the user input is larger or equal to the maximum chromosome number
                    if cmax <= f_s1:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        playing
                    #break loop if input is adequate
                    else:
                        playing = False
            #set executing equal to true for while loop
            active = True
            while active: #option 2 first set of ending coordinates
                #take in coordinate as a string from user
                end_coord1 = input('input the first ending coordinate (hg38): ')
                #create a list of individual strings for each character
                e_nums1 = [x for x in end_coord1]
                #create a while loop to remove commas from user input
                while ',' in e_nums1:
                    e_nums1.remove(',')
                #join user input string without commas
                e_str1 = ''.join(e_nums1)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (e_str1.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    active
                #create elif statement where input is only numbers
                elif (e_str1.isnumeric() == True):
                    #convert string to integer
                    f_e1 = int(e_str1)
                    #create an if statement if the user input is larger than the maximum chromosome number
                    if cmax < f_e1:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        active
                    #create an elif statement if user input is equal to staarting coordinate
                    elif f_e1 == f_s1:
                        print('this number is equal to the starting coordinate. try again.')
                        #reloop through while to promt user to input another value
                        active
                    #create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e1 < f_s1:
                        print('this number is smaller than the starting coordinate. try again.')
                        active
                    #create an else statemnt to break loop when input is adequate
                    else:
                        active = False
                        working = False
            alert = True
            while alert: #option 2 second starting coordinate
                #take in coordinate as a string from user
                start_coord2 = input('input the second starting coordinate: ')
                #create a list of individual strings for each character
                s_nums2 = [x for x in start_coord2]
                #create a while loop to remove commas from user input
                while ',' in s_nums2:
                    s_nums2.remove(',')
                #join user input string without commas
                s_str2 = ''.join(s_nums2)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (s_str2.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    alert
                #create elif statement where input is only numbers 
                elif (s_str2.isnumeric() == True):
                    #convert string to integer
                    f_s2 = int(s_str2)
                    #create an if statement if the user input is larger or equal to the maximum chromosome number
                    if cmax <= f_s2:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        alert
                    #break loop if input is adequate
                    else:
                        alert = False
            conducting = True
            while conducting: #option 2 second set ending coordinate
                #take in coordinate as a string from user
                end_coord2 = input('input the second ending coordinate: ')
                #create a list of individual strings for each character
                e_nums2 = [x for x in end_coord2]
                #create a while loop to remove commas from user input
                while ',' in e_nums2:
                    e_nums2.remove(',')
                #join user input string without commas
                e_str2 = ''.join(e_nums2)
                #create an if statement for the user inputting a string that doesnt consist of only numbers
                if (e_str2.isnumeric() == False):
                    print('input must consist of only numbers. try again.')
                    #reloop though while to prompt user to enter another value
                    conducting
                #create elif statement where input is only numbers
                elif (e_str2.isnumeric() == True):
                    #convert string to integer
                    f_e2 = int(e_str2)
                    #create an if statement if the user input is larger than the maximum chromosome number
                    if cmax < f_e2:
                        print('this number is too large. try again.')
                        #reloop through while to prompt user to input another value
                        conducting
                    #create an elif statement if user input is equal to staarting coordinate
                    elif f_e2 == f_s2:
                        print('this number is equal to the starting coordinate. try again.')
                        #reloop through while to promt user to input another value
                        conducting
                    #create an elif statement if the ending coordinate is smaller than the starting coordinate
                    elif f_e2 < f_s2:
                        print('this number is smaller than the starting coordinate. try again.')
                        conducting
                    #create an else statemnt to break loop when input is adequate
                    else:
                        conducting = False
                        working = False
        #create else for invalid inputs
        else:
            print("invalid input. try again.")
            working

    managing = True
    while managing: #option two stranded options
        #take in user input 
        choice = input('enter a "+" for a positive strand and a "-" for a negative strand: ') #pos or neg strand
        #create an if statement for positive strand 
        if (choice == "+"):
            strand = 'pos'
            #break loop for correct input
            managing = False
        #create an elif statement for negative strand
        elif (choice == "-"):
            strand = 'neg'
            #break loop for correct input
            managing = False
        #create an else statement for incorrect input
        else:
            print('invalid input. try again.')
            #reloop through while to prompt user to input another value
            managing

    #read canonical bed file under the variable anno
    anno = pd.read_csv(gene_names_bed_path, sep = '\t')
    #create an empty gene list
    gene_list = []
    #create an if condition for gene name if there is only 1 set of coordinates
    if coord_opt == '1':
        #iterate through the data frame to find gene that satisfies all user inputs.
        for i in range(0, len(anno)):
            if (((((anno['start'][i] < f_s) & (anno['end'][i] >= f_s)) | ((anno['start'][i] < f_e) & (anno['end'][i] >= f_e)))) & ((anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                #append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
        #set conditions to display gene(s) to the user
        gene_list = [x for x in gene_list if str(x) != 'nan']
        if len(gene_list) == 1:
            for i in range(len(gene_list)):
                print('based on canonical annotations, the following gene is in your area of interest: ', gene_list[i])
        elif len(gene_list) > 1:
            total = ', '.join(gene_list)
            print('based on canonical annotations, the following genes are in your area of interest: ', total)
        #set condition when there aren't any genes in the area of interest
        elif len(gene_list) == 0:
            print('based on canonical annotations, there are not any genes in your area of interest.')

    #create an if condition for gene name if there are 2 sets of coordinates
    if coord_opt == '2':
        #iterate through the data frame to find gene that satisfies all user inputs.
        for i in range(0, len(anno)):
            if (((((anno['start'][i] < f_s1) & (anno['end'][i] >= f_s1)) | ((anno['start'][i] < f_e1) & (anno['end'][i] >= f_e1)))) & ((anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                #append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
            if (((((anno['start'][i] < f_s2) & (anno['end'][i] >= f_s2)) | ((anno['start'][i] < f_e2) & (anno['end'][i] >= f_e2)))) & ((anno['chrom'][i] == chrom) & (anno['strand'][i] == choice))):
                #append all possible genes to the gene list
                gene_list.append(anno['gene'][i])
        gene_list = [x for x in gene_list if str(x) != 'nan']
        #set conditions to display gene(s) to the user
        if len(gene_list) == 1:
            for i in range(len(gene_list)):
                print('based on canonical annotations, the following gene is in your area of interest: ', gene_list[i])
        elif len(gene_list) > 1:
            total = ','.join(gene_list)
            print('based on canonical annotations, the following genes are in your area of interest: ', total)
        #set condition when there aren't any genes in the area of interest
        elif len(gene_list) == 0:
            print('based on canonical annotations, there are not any genes in your area of interest.')

    #change directory to folder with pileups
    os.chdir(bgdir)

    #create an if when user selects only one set of coordinates
    if coord_opt == '1':
        #set if statement for when the strand is positive 
        if strand == 'pos':
            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it 
            #drop the start column and only include the pos for the coordinates and D1 column for coverage 
            D1 = (pd.read_csv('D1_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D1'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D1 = D1.loc[D1['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D1

            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D2 column for coverage 
            D2 = (pd.read_csv('D2_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D2'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D2 = D2.loc[D2['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D2

            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D3 column for coverage 
            D3 = (pd.read_csv('D3_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D3'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D3 = D3.loc[D3['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D3

        #create an if statement for when the stand is negative
        if strand == 'neg':
            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D1 column for coverage
            D1 = (pd.read_csv('D1_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D1'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D1 = D1.loc[D1['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D1

            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D2 column for coverage
            D2 = (pd.read_csv('D2_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D2'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D2 = D2.loc[D2['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D2

            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D3 column for coverage
            D3 = (pd.read_csv('D3_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D3'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            gene_D3 = D3.loc[D3['pos'].between(f_s, f_e)].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D3

    #create an elif for when user slects to have two sets of coordinates
    elif coord_opt == '2':
        #set if statement for when the strand is positive
        if strand == 'pos':
            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D1 column for coverage
            D1 = (pd.read_csv('D1_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D1'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D1 = D1.loc[(D1['pos'].between(f_s1, f_e1)) | (D1['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            del D1

            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D2 column for coverage
            D2 = (pd.read_csv('D2_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D2'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D2 = D2.loc[(D2['pos'].between(f_s1, f_e1)) | (D2['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D2

            #make the positive bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D3 column for coverage
            D3 = (pd.read_csv('D3_' + chrom + '_mmRate_pos.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D3'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D3 = D3.loc[(D3['pos'].between(f_s1, f_e1)) | (D3['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D3

        #create an if statement for when the strand is negative
        if strand == 'neg':
            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D1 column for coverage
            D1 = (pd.read_csv('D1_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D1'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D1 = D1.loc[(D1['pos'].between(f_s1, f_e1)) | (D1['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D1

            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D2 column for coverage
            D2 = (pd.read_csv('D1_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D2'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D2 = D2.loc[(D2['pos'].between(f_s1, f_e1)) | (D2['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D2

            #make the negative bedgraph file for the respective chromosome a data frame and iterate through it
            #drop the start column and only include the pos for the coordinates and D3 column for coverage
            D3 = (pd.read_csv('D3_' + chrom + '_mmRate_neg.bedGraph', sep = '\t', header = None, names = ['start', 'pos','D3'])).drop(columns=['start'])
            #splice the data frame to only include the region the user entered
            #for two corrdinates, ensure and operator splices for each coordinate set
            gene_D3 = D3.loc[(D3['pos'].between(f_s1, f_e1)) | (D3['pos'].between(f_s2, f_e2))].reset_index(drop=True)
            #delete the unspliced sample from the data drame
            del D3
    #put them all in one data frame gene_D1, gene_D2, gene_D3
    gene_strand = pd.DataFrame()
    #assign pos based on the pos of the first spliced sample
    gene_strand['pos'] = gene_D1['pos']
    #obtaiin the average of the 3 samples to represent D
    gene_strand['D'] = (gene_D1['D1'] + gene_D2['D2'] + gene_D3['D3']) / 3

    # read in fasta file for sequence and add ref column to df
    abs_seq_file_path = os.path.abspath(f"{seq_fasta_dir}{chrom}.fasta")
    seq_file = open(abs_seq_file_path)
    seq_lines = seq_file.readlines()
    seq = seq_lines[1]
    seq = seq.strip()
    seq_list = [letter for letter in seq]
    seq_file.close()

    #this function removes the top 5% of outliers and scales that output to a max of 1
    def winsor_scale(col):
        arr = np.array(col)
        win_arr = winsorize(arr, limits = [0.00, 0.05]).data #winsorize with limits of 0%, 95%.
        #FOR VERY SMALL WINDOWS, WE WILL NEED TO SHOW A WARNING IF THE HIGHEST VALLUE IS AN OUTLIER.
        norm_arr = (win_arr - win_arr.min())/(win_arr.max() - win_arr.min()) #scale to 1
        return norm_arr

    #use the windor function on column D
    gene_strand['D'] = winsor_scale(gene_strand['D'])

    #create if statement when user only inputs 1 coordinate
    if coord_opt == '1':
        #create the data frame RNAstructure
        RNAstructure = pd.DataFrame()
        #iterate through the user input for start and end coordinates and append them as column pos in the df
        RNAstructure['pos'] = range(f_s, f_e + 1)
        #empty coord column
        RNAstructure['coord'] = ''
        #append ref column to RNAstructure to assign base letter to each coordinate
        RNAstructure['ref'] = seq_list[f_s-1:f_e]
        #merge RNAstructure with gene_strand to append coverage values to RNAstructure
        RNAstructure = RNAstructure.merge(gene_strand, on=['pos'], how='left')
    #create if statement when the user input is two coordinates
    if coord_opt == '2':
        #create a list that includes coordinates from both sets of user inputs
        pos_list = list(range(f_s1, f_e1 + 1)) + list(range(f_s2, f_e2 + 1))
        #create the data frame RNAstructure
        RNAstructure = pd.DataFrame()
        #create column pos from pos_list values to append all coordinates to the df
        RNAstructure['pos'] = pos_list
        #empty coord column
        RNAstructure['coord'] = ''
        #create two_seq to append the sequence from both sets of coordinates
        two_seq = seq_list[f_s1-1:f_e1] + seq_list[f_s2-1:f_e2]
        #append two_seq to the df to assign a base letter to each coordinate
        RNAstructure['ref'] = two_seq
        #merge RNAstructure with gene_strand to append coverage values to RNAstructure
        RNAstructure = RNAstructure.merge(gene_strand, on=['pos'], how='left')

    #for reversely stranded genes, reverse the order of the bases and replace them with their complements
    def rev_strand(df):
        df.sort_values('pos', ascending = False, inplace=True).reset_index(drop=True, inplace=True) #reverses the order
        for i in range(0, len(df)): #replace with their complements
            if df['ref'][i] == "A":
                df['ref'][i] = 'U'
            elif df['ref'][i] == 'T':
                df['ref'][i] = 'A'
            elif df['ref'][i] == 'G':
                df['ref'][i] = 'C'
            elif df['ref'][i] == 'C':
                df['ref'][i] = 'G'
        return df

    #create an if statement to reverse stand if strand is negative
    if strand == 'neg':
        RNAstructure.sort_values('pos', ascending = False, inplace=True)
        RNAstructure.reset_index(drop=True, inplace=True) #reverses the order
        for i in range(0, len(RNAstructure)): # replace with thier complements
            if RNAstructure.loc[i, 'ref'] == "A":
                RNAstructure.loc[i, 'ref'] = 'U'
            elif RNAstructure.loc[i, 'ref'] == 'T':
                RNAstructure.loc[i, 'ref'] = 'A'
            elif RNAstructure.loc[i, 'ref'] == 'G':
                RNAstructure.loc[i, 'ref'] = 'C'
            elif RNAstructure.loc[i, 'ref'] == 'C':
                RNAstructure.loc[i, 'ref'] = 'G'
    elif strand == 'pos': #for the pos strand, the order and complement of the sequence is fine, but T must be replaced with U
        RNAstructure.reset_index(inplace=True, drop=True)
        for i in range (0, len(RNAstructure)):
            if RNAstructure['ref'][i] == 'T':
                RNAstructure['ref'][i] = 'U'

    #fill coord column starting from 1
    RNAstructure['coord'] = np.arange(1, len(RNAstructure) + 1)
    RNAstructure.head()

    #create data frame for plot only including A's and C's
    plot = pd.DataFrame()
    plot = RNAstructure[(RNAstructure['ref'] == 'A') | (RNAstructure['ref'] == 'C')]

    #replace underscore used in files for a space
    if '_' in AOI:
        AOI = original.replace('_', ' ')

    #plot using RNAstructure frame with d as y and coord as x
    coordinates = plot['coord']
    mismatch = plot['D']
    plt.figure(figsize = (10,7))
    plt.bar(coordinates, mismatch)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylim(0,1)
    plt.xlabel("coordinates", fontsize=14)
    plt.ylabel("mismatch rates", fontsize=14)
    plt.title(AOI, fontsize=18)
    #save plot in files
    plt.savefig(barplotdir + '/' + 'plot.svg')

    #anything above 25%, warning that their might be an outlier
    #at least 40 valid data points
    if max(RNAstructure['D']) >= 0.2 and len(gene_strand) < 40:
        print('warning: there might be an outlier in this data. this may affect normalization for small areas of interest.')

    RNAstructure.reset_index(drop=True, inplace=True)

    #use fill na to change values to -999
    RNAstructure['D'] = RNAstructure['D'].fillna(-999)

    seq = ''.join(RNAstructure.ref.tolist()) #creates a string with all bases included in the region

    #change to the fasta directory
    os.chdir(fastadir)

    #place the underscore back in the title
    if ' ' in AOI:
        AOI = AOI.replace(' ', '_')

    #writes the fasta file containing the header and sequence
    text_file = open(AOI + '.fasta', 'w')
    text_file.write('>' + AOI + '\n')
    text_file.write(seq)
    text_file.close()

    #change to the .dat directory
    os.chdir(datdir)
    #save the .dat file(s)
    #only uses the coord and condition-average reactivity columns
    #use rnastrusture to create a .dat file using RNAstructure df
    RNAstructure.to_csv(AOI + '_D.dat', sep='\t', header=False, index=False, columns = ['coord', 'D'])

    #attach sample name to all files
    #fold command line
    os.system(fold_bin_path + fastadir + '/' + AOI + '.fasta ' + ctdir + '/' + AOI + '_fold.ct ' + '--SHAPE ' + datdir + '/' + AOI + '_D.dat')

    #run command to obtain free folding energy

    os.system(efn2_bin_path + ' ' + ctdir + '/' + AOI + '_fold.ct ' + fold_FEdir + '/' + AOI + '.txt --SHAPE ' + datdir + '/' + AOI + '_D.dat')

    #read the number of lines in the energy file to display the amount of structures to the user
    FE_lines = open(fold_FEdir + '/' + AOI + '.txt')
    struct_num = len(FE_lines.readlines())
    print('this area of interest has', struct_num, 'structures. how many would you like to visualize? : ')
    user_int = int(input())
    for x in range(1, user_int + 1):
        user_struct = str(x)
        FE_lines2 = open(fold_FEdir + '/' + AOI + '.txt')
        ttl = FE_lines2.readlines()
        y = ttl[x-1]
        ttl2 = y.split(' ')
        s_num = ttl2[1]
        energy = ' '.join(ttl2[4:])
        #obtain dbn file for specified structure
        #ct2dot ct to djbn (fold)
        os.system(ct2dot_bin_path + ' ' + ctdir + '/' + AOI + '_fold.ct ' + user_struct + ' ' + fold_dbndir + '/' + AOI + '_' + user_struct + '.dbn')

        #obtain coverage percentage 
        d = len(gene_strand)/len(RNAstructure)
        percentage = round(d,2)
        #multiply by 100 to obtain percentage
        cov = 100 * percentage

        #begin while for coverage
        finishing = True
        while finishing:
            #create if statement when coverage is less than 30%
            if cov < 30:
                #display coverage to the user 
                print('your DMS coverage is', cov, '%')
                #take in user input yes/no to detrmine whether or not to proceed with imaging
                choice2 = input('since your DMS coverage is less than 30 %, do you wish to proceed?: ')
                #make user input all uppercase to eliminate case sensitivity
                choice2 = choice2.upper()
                if choice2 == 'YES':
                    finishing = False
                elif choice2 == 'NO':
                    print('ending program...')
                    sys.exit(0)
                else:
                    print('invalid input. try again.')
                    finishing = True
            #create an elif statemnt when the coverage is greater than 30% or the user inputs yes
            elif cov >= 30:
                print('your DMS coverage is', cov, '%')
                finishing = False

        #open the new dbn file
        fold_dbn = open(fold_dbndir + '/' + AOI + '_' + user_struct + '.dbn')
        #read third line for dbn notation as list
        #-1 for extra space at the end of line
        struct1 = fold_dbn.readlines() [2][:-1]

        #create a lst1 for g's and u's and lst 2 for a's and c's
        lst1 = []
        lst2 = []

        #iterate though list of individual letters
        for val in range(len(RNAstructure)):
        #append the base number of g's and u's for'
        #final version shouldn't have to have the T line
            if RNAstructure['D'][val] == -999:
                lst1.append(val)
            else:
                lst2.append(val)
        #create list low, mid and high for low, mid, and high coverage values
        low = []
        mid = []
        high = []
        #append each value to their respective lists based on coverage
        for val in lst2:
            if (0.00 <= RNAstructure['D'][val] < 0.25):
                low.append(val)
            elif (0.25 <= RNAstructure['D'][val] < 0.5):
                mid.append(val)
            elif (0.5 <= RNAstructure['D'][val] <= 1.0):
                high.append(val)

        #change directory for varna imaging
        os.chdir(wkdir)

        #get image for fold dbn
        from varnaapi import VARNA
        import varnaapi
        v = VARNA(structure = struct1, sequence = seq)
        style1 = varnaapi.param.BasesStyle(fill="#FFFFFF")
        style2 = varnaapi.param.BasesStyle(fill="#00FFFF")
        style3 = varnaapi.param.BasesStyle(fill="#FFFF00")
        style4 = varnaapi.param.BasesStyle(fill="#FF0000")
        v.add_bases_style(style1, lst1)
        v.add_bases_style(style2, low)
        v.add_bases_style(style3, mid)
        v.add_bases_style(style4, high)
        if '_' in AOI:
            AOI = AOI.replace('_', ' ')
        v.set_title(AOI + ', structure ' + s_num + ', ' + energy + ' kcal/mol', color='#000000', size=12)
        if ' ' in AOI:
            AOI = AOI.replace(' ', '_')
        v.savefig(final_image + '/' + AOI + '_fold_' + user_struct + '.svg')
        fig = sg.fromfile(final_image + '/' + AOI + '_fold_' + user_struct + '.svg')
        fig.set_size(('2000', '2000'))
        fig.save(final_image + '/' + AOI + '_fold_final_' + user_struct + '.svg')
        #break to avoid infinite loop
        finishing = False

def abs_binary_path(name):
    which_path = shutil.which(name)
    if which_path is None:
        print(f"Error: The {name} binary was not found. Please ensure it is installed and available in the system's PATH.")
        sys.exit(1)

    return os.path.abspath(which_path)

def ensure_directories_exist(directories):
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

# Run if called directly.
if __name__ == "__main__":
    main()
