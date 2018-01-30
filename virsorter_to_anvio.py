#!/usr/bin/env python

import collections
import argparse
import sys

"""
Parses viral contig prediction results from VirSorter and generates two files usable by Anvi'o:
1. "virsorter_additional_info.txt" which you can import as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.
2. "virsorter_collection.txt" which is a collection file that will automatically group splits belonging to each VirSorter phage prediction into its own bin.
"""
parser = argparse.ArgumentParser(description="Parses VirSorter predictions for Anvi\'o.")
parser.add_argument("-a","--affi-file", help = "VIRSorter_affi-contigs.tab file. REQUIRED.")
parser.add_argument("-g","--global-file", help="VIRSorter_global_signal.csv file. REQUIRED.")
parser.add_argument("-s","--splits-info", help = "splits_basic_info.txt file. REQUIRED.")
parser.add_argument("-l","--min-phage-length", type = int, default = 1000, help = "Specify the minimum phage length to report. Default is 1000 bp.")
parser.add_argument("--exclude-cat3", action = "store_true", help = "Excludes all category 3 predictions from both output files.")
parser.add_argument("--exclude-prophages", action = "store_true", help = "Exclude all prophage predictions.")
parser.add_argument("-A","--addl-info", default = "virsorter_additional_info.txt", help = "Additional info output file. The default file name is \"virsorter_additional_info.txt\". You can import this as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.")
parser.add_argument("-C","--phage-collection", default = "virsorter_collection.txt", help = "Outputs an Anvi\'o collections file with splits for each phage gathered into a separate bin. The default name is \"virsorter_collection.txt\".")
args = parser.parse_args()
arg_dict = vars(args)
#print(arg_dict)

if arg_dict['affi_file'] == None or arg_dict['global_file'] == None or arg_dict['splits_info'] == None:
	print("\n***A required option is missing. Try again.***\n\n")
	parser.print_help()
	sys.exit()

#PART ONE
#First, I need to read in virsorter_affi, the input to Step3 in Virsorter.
#I need this because I need to translate the gene boundaries of each predicted phage
#to start and stop positions (bp) along the contig.

#columns of VIRSorter_affi-contigs.tab:
#     0  | 1   | 2  |  3   |  4   |    5     |  6  |   7  |   8     |   9    | 10   | 11
# gene_id|start|stop|length|strand|affi_phage|score|evalue|category|affi_pfam|score|evalue|

virsorter_affi = open(arg_dict['affi_file'], "r")
affi_dict = {}
affi_dict = collections.defaultdict(dict)
line = virsorter_affi.readline()
while line != "" and line[0] == ">":
    header = line[1:].strip().split('|')
    line = virsorter_affi.readline()
    while line != "" and line[0] != ">":
        var = line.strip().split('|')
        gene = var[0].split('-')
        gene = gene[-1]
        gene_start = var[1]
        gene_stop = var[2]
        header2 = header[0].replace("VIRSorter_","").replace("-circular","")
        #do something about if the gene wraps around... circular... gah...
        affi_dict[header2][gene] = {}
        affi_dict[header2][gene]['start'] = gene_start
        affi_dict[header2][gene]['stop'] = gene_stop
        line = virsorter_affi.readline()
virsorter_affi.close()

#PART TWO
#Now we read in the final output of VirSorter, which is VIRSorter_global_signal.csv
#and save it as a dictionary.

#From this, we need: (0) contig_id, (1) genes in contig, (2) Fragment, (3) Nb genes, (4) Category, (5) nb_hallmark
#Column types for virsorter_global:
#    0   /       1        /    2       /  3  /     4          /    5      /      6        /     7        /     8       /       9        /     10   /      11     /    
# Contig / Total Nb Genes /  Fragment / Size / Type detection / Category /  Enrich Phage / Enrich Pfam / Enrich Unch / Enrich Switch / Avg_g_size / Nb Hallmark

virsorter_global = open(arg_dict['global_file'], "r")

global_dict = {}
global_dict = collections.defaultdict(dict)
m = 1
n = 1
line = virsorter_global.readline()
while line != "":
    if line[0] == "#":
        line = virsorter_global.readline()
    else:
        line = line.strip().split(',')
        #print line
        contig = line[0].replace("VIRSorter_","").replace("-circular","")
        genes_in_contig = line[1]
        fragment = line[2].replace("VIRSorter_","").replace("-circular","")
        num_fragment_genes = int(line[3])
        
        nb_hallmark = line[5]
        if nb_hallmark == "":
            nb_hallmark = 0
        nb_hallmark = int(nb_hallmark)

        global_dict[contig]['num_genes_in_contig'] = int(genes_in_contig)
        global_dict[contig]['fragment'] = fragment
        global_dict[contig]['num_fragment_genes'] = int(num_fragment_genes)
        global_dict[contig]['nb_hallmark'] = nb_hallmark
        global_dict[contig]['length'] = 0
        
        if "-circular" in line[0]:
            global_dict[contig]['circular'] = True
        else:
            global_dict[contig]['circular'] = False
        
        if contig == fragment:
            #print category
            category = int(line[4])
            global_dict[contig]['category'] = str(category)  #This is messed up somehow... refreshing helped things? Dunno.
            global_dict[contig]['phage_num'] = "phage_"+str(m)
            m += 1
        if contig != fragment:
            #print category
            category = int(line[4])
            global_dict[contig]['category'] = str(category)
            genes_in_fragment = fragment.split('-')
            start_gene = genes_in_fragment[-2]
            stop_gene = genes_in_fragment[-1]
            global_dict[contig]['start_gene'] = start_gene
            global_dict[contig]['stop_gene'] = stop_gene
            global_dict[contig]['start_gene_pos'] = int(affi_dict[contig][start_gene]['start'])
            global_dict[contig]['stop_gene_pos'] = int(affi_dict[contig][stop_gene]['stop'])
            global_dict[contig]['phage_num'] = "prophage_"+str(n)
            n += 1

        line = virsorter_global.readline()
virsorter_global.close()

#PART THREE
#This writes all phages and prophages, cat1-3, to additional-info and collections files.
splits_input = open(arg_dict['splits_info'], "r")
splits_output = open(arg_dict['addl_info'],'w')
collection_output = open(arg_dict['phage_collection'],'w')

#splits_basic_info column format:
#  0            1             2     3      4          5               6              7
#split   order_in_parent   start   end   length   gc_content   gc_content_parent   parent

#is_prev_phage = False
#n = 1

splits_output.write("split\tphage_name\tphage_category\tphage_length\tnum_genes_in_phage\tnum_phage_hallmark_genes_in_phage\n")
line = splits_input.readline()
line = splits_input.readline()

split_length_dict = collections.defaultdict(dict)
while line != "":
    line = line.strip().split('\t')
    if 'length' not in split_length_dict[line[7]]:
        split_length_dict[line[7]]['length'] = 0
    split_length_dict[line[7]]['length'] += int(line[4])
    line = splits_input.readline()
splits_input.close()

splits_input = open(arg_dict['splits_info'], "r")
line = splits_input.readline()
line = splits_input.readline()
while line != "":
    line = line.strip().split('\t')
#    print line[7]
    split_name = line[0]
    split_start = int(line[2])
    split_stop = int(line[3])
    split_parent = line[7]
    parent_length = split_length_dict[split_parent]['length']
    if split_parent in global_dict:
        nb_hallmark = global_dict[split_parent]['nb_hallmark']
        nb_genes = global_dict[split_parent]['num_fragment_genes']
        phage_name = global_dict[split_parent]['phage_num']
        if split_parent == global_dict[split_parent]['fragment']:
            category = "cat"+str(global_dict[split_parent]['category'])+"_phage"
            global_dict[split_parent]['length'] = int(parent_length)
            phage_length = global_dict[split_parent]['length']
        if split_parent != global_dict[split_parent]['fragment']:
            #This is where it gets hairy....
            start_gene_pos = global_dict[split_parent]['start_gene_pos']
            stop_gene_pos = global_dict[split_parent]['stop_gene_pos']
            if start_gene_pos <= split_start <= stop_gene_pos or start_gene_pos <= split_stop <= stop_gene_pos:
                category = "cat"+str(global_dict[split_parent]['category'])+"_prophage"
                global_dict[split_parent]['length'] = stop_gene_pos - start_gene_pos
                phage_length = global_dict[split_parent]['length']
            else:
                category = ""
                nb_genes = 0
                nb_hallmark = 0
                phage_name = ""
                phage_length = 0
    else:
        nb_genes = 0
        nb_hallmark = 0
        category = ""
        phage_name = ""
        phage_length = 0
    
    #Output format: split_name  phage_name  category  phage_length  nb_genes   nb_hallmark
    if phage_length >= int(arg_dict['min_phage_length']):
        if arg_dict['exclude_cat3'] == True:
            if category == "cat1_phage" or category == "cat2_phage" or category == "cat1_prophage" or category == "cat2_prophage":
                if arg_dict['exclude_prophages'] == True:
                    if category == "cat1_phage" or category == "cat2_phage":
                        splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            split_name, phage_name, category, str(phage_length), str(nb_genes), str(nb_hallmark)))
                        collection_output.write("%s\t%s\n" %(split_name, phage_name))
                else:
                    splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        split_name, phage_name, category, str(phage_length), str(nb_genes), str(nb_hallmark)))
                    collection_output.write("%s\t%s\n" %(split_name, phage_name))
            ### THIS PART NEEDS LOTS OF WORK TO GET ALL ARGS TO WORK CORRECTLY. ###
        
        else:
            if arg_dict['exclude_prophages'] == True:
                if category == "cat1_phage" or category == "cat2_phage" or category == "cat3_phage":
                        splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            split_name, phage_name, category, str(phage_length), str(nb_genes), str(nb_hallmark)))
                        collection_output.write("%s\t%s\n" %(split_name, phage_name))
            else:
                splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    split_name, phage_name, category, str(phage_length), str(nb_genes), str(nb_hallmark)))
                collection_output.write("%s\t%s\n" %(split_name, phage_name))

    line = splits_input.readline()

splits_input.close()
splits_output.close()
collection_output.close()

