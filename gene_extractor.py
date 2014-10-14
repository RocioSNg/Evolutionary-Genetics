#-------------------------------------------------------------------------------
# Name:         Gene Sequence Extractor
# Purpose:      Extracts gene sequences from Genome FASTQ files using chromosome positions
#               Places all sequences into a single file that can be used for further analysis
#               *Keeps quality scores intact*
#               *Use with Python Version 2.7*
# Author:       Rocio Ng
# Created:      15/09/2014
#-------------------------------------------------------------------------------


##-----IMPORT NEEDED MODULES------##

import os
import glob

##---DEFINE NEEDED VARIABLES------##


input_folder ="ChrX_Genome"
gene_folder = "3_Extracted_Sequences"
query = open("2_Query_Files\\tan_CDS_query.txt", "r")
#print "What is name of the folder containing the input files?"
#input_folder = raw_input()
#print "Where would you like the output file to go?"
#gene_folder = raw_input()
#print "What is the name of the Query File (include extension)?"
#query = open(raw_input(), 'r')
#query = open("sample_query.txt")

#---READ the parameters in the QUERY FILE----#
#   TITLE
#   GENE NAME
#   DATA TYPE - .txt or .fq
#   OUTPUT TYPE - gene or CDS
#   START POSITION
#   END POSTION
#--------------------------------------------#
query_lines = query.readlines()
gene = query_lines[1]
gene = gene.strip()
file_type = query_lines[2]      #.txt or .fq
file_type = file_type.strip()   #strip gets rid of /n
output_type = query_lines[3]    # gene or CDS
output_type= output_type.strip()

#----------------WHOLE GENE SEQUENCES-------------------------#
if output_type == 'gene':
    s = query_lines[4]
    s = eval(s)
    e = query_lines[5]
    e = eval(e)

#-----------------CDS SEQUENCE ONLY--------------------------#

if output_type == 'CDS':

    positions = list()
    for i in range(4, len(query_lines)):  #put all the numbers in a neat list
        line = query_lines[i]
        line= eval(line)
        positions.append(line)

#print "What is the name of the gene?"
#gene = raw_input()
#print "Where does the gene start?"
#s = int(raw_input())
#print "Where does the gene end?"
#e = int(raw_input())
#print ".txt or .fq?"
#file_type = raw_input()
#print "Gene or CDS?"
#output_type = raw_input()


#-------IDENTIFYING WHERE THE INPUT AND OUTPUT FOLDERS ARE---------------

data_dir = os.getcwd()      #string that defines working directory
input_dir = data_dir + '\\' + input_folder + '\\'  #path to where input files are held
gene_dir = data_dir + '\\' + gene_folder + '\\' + gene + '\\'  #path to where the files will be imported into.  Includes gene folder name

##--check to make sure paths are correctly defined
print input_dir
print gene_dir


file_list = os.listdir(input_dir)

"""
        this will create a list of all file names
        Note: Indexing starts at 0
        txt_list[0] reutrns the first file in the folder
"""

##---FUNCTION that takes gene coordinates and indexes the FAST Q files to only include it

def gene_extract(file_name,gene, start, end):
    "gene is the name of the gene"
    "start and end of the gene location as found in flybase"

    name = gene + '_' + file_name       #--determines how to name the final file
    name = name.replace(file_type, '') #--gets rid of the extension

    #---convert positions to python indexing methods----#
    t_start = start -1
    t_end = end     # b/c of range rules

    data = open(input_dir + file_name, 'r') #opens the data file
    txt_lines = data.readlines()    #puts all the lines in a variable

    nucleotide_seq = txt_lines[1]    #extracts the lines with the nucleotide seqeunces
    gene_seq = nucleotide_seq[t_start: t_end]  #indexes based on gene location

    quality_scores = txt_lines[3]       #extracts the line with the Quality scores
    gene_qscore = quality_scores[t_start: t_end]    #indexes based on gene location

    #Create new files with extracted genes and place into the proper folder
    new_file = open(gene_dir + name + '.txt','w')

    new_file.write(txt_lines[0])  #name of the line and chrom
    new_file.write(gene_seq)    #subsetted sequence
    new_file.write('\n')
    new_file.write(txt_lines[2])  #name of line
    new_file.write(gene_qscore)  #substted quality scores
    new_file.write('\n')
    new_file.close()

    print name + " is created!"


def CDS_extract(file_name):
    "gene is the name of the gene"
    "start and end of the gene location as found in flybase"

    name = gene + '_CDS_' + file_name       #--determines how to name the final file
    name = name.replace(file_type, '') #--gets rid of the extension


    data = open(input_dir + file_name, 'r') #opens the GENOME data file
    txt_lines = data.readlines()    #puts all the lines in a variable- FOUR lines in total

    nucleotide_seq = txt_lines[1]    #extracts the lines with the nucleotide seqeunces
    quality_scores = txt_lines[3]       #extracts the line with the Quality scores

    CDS_seq = ""       # will collect nucleotides
    CDS_qscore = ""    #will collect Qscores



    k = 0   #start with the first fragment

    for j in range(1, len(positions)/2 + 1):    #range is not all inclusive
        s = positions[k] -1        # start of fragment --subtract one of of python's illogical indexing rules
        e = positions[k+1]        # end of fragment
        CDS_seq += nucleotide_seq[s:e]
        CDS_qscore += quality_scores[s:e]
        k += 2      #skips to the next pair of position coordinates
        #--***Positions defined above (list of all numbers)**

    print name + " is created!"

    #Create new files with extracted genes and place into the proper folder
    new_file = open(gene_dir + name + '.txt','w')

    new_file.write(txt_lines[0])  #name of the line and chrom
    new_file.write(CDS_seq)    #collected CDS fragments
    new_file.write('\n')
    new_file.write(txt_lines[2])  #name of line again
    new_file.write(CDS_qscore)  #collected Qscore fragments
    new_file.write('\n')
    new_file.close()

    print name + " is created!"




##----LOOPS that carries out extract functions for every file in a folder


if output_type == "gene":
    print "Whole gene sequences will be extracted!"
    for file in file_list:
        gene_extract(file, gene, s, e)

if output_type == "CDS":
    print "CDS gene sequences will be extracted!"
    for file in file_list:
        CDS_extract(file)

#---------COMBINE ALL OUTPUT FILES INTO ONE *.TXT file---------#

files = glob.glob(gene_dir + '*.txt' )

with open(gene + '_' + output_type + '.txt', 'w' ) as result:
    for file_ in files:
        for line in open( file_, 'r' ):
            result.write( line )