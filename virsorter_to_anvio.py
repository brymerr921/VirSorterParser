#!/usr/bin/env python

import collections
import argparse
import sys

"""
Parses viral contig prediction results from VirSorter and generates three files usable by Anvi'o:
1. "virsorter_additional_info.txt" which you can import as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.
2. "virsorter_collection.txt" which is a collection file that will automatically group splits belonging to each VirSorter phage prediction into its own bin. This can be imported using anvi-import-collection.
3. "virsorter_annotations.txt" which is an annotations file containing predicted functions for hallmark genes. This can be imported using anvi-import-functions.
"""
"""
Parses viral contig prediction results from VirSorter and generates three files usable by Anvi'o:
1. "virsorter_additional_info.txt" which you can import as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.
2. "virsorter_collection.txt" which is a collection file that will automatically group splits belonging to each VirSorter phage prediction into its own bin. This can be imported using anvi-import-collection.
3. "virsorter_annotations.txt" which is an annotations file containing predicted functions for hallmark genes. This can be imported using anvi-import-functions.
"""
parser = argparse.ArgumentParser(description="Parses VirSorter predictions for Anvi\'o.")
parser.add_argument("-a","--affi-file", help = "REQUIRED INPUT. VIRSorter_affi-contigs.tab file.")
parser.add_argument("-g","--global-file", help="REQUIRED INPUT. VIRSorter_global_signal.csv file.")
parser.add_argument("--db", help="REQUIRED INPUT. Which VirSorter database you used. Enter 1 for RefSeq, 2 for RefSeq+viromes.")
parser.add_argument("-s","--splits-info", help = "REQUIRED INPUT. splits_basic_info.txt file.\n\n")
parser.add_argument("-n","--anvio-gene-calls", help = "OPTIONAL INPUT, but REQUIRED if you want an output functions file. Anvi'o gene calls exported from the contigs database")
parser.add_argument("-f","--hallmark-functions", help = "OPTIONAL INPUT, but REQUIRED if you want an output functions file. Functional annotation for protein clusters of hallmark genes. You'll want to specify the db1 file if you ran VirSorter against db1, and the db2 file if you ran VirSorter against db2.")
parser.add_argument("-A","--addl-info", default = "virsorter_additional_info.txt", help = "OUTPUT. Additional info file. Default = \"virsorter_additional_info.txt\". You can import this as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.")
parser.add_argument("-C","--phage-collection", default = "virsorter_collection.txt", help = "OUTPUT. An Anvi\'o collections file with splits for each phage gathered into a separate bin. The default name is \"virsorter_collection.txt\".")
parser.add_argument("-F","--output-functions", default = "virsorter_annotations.txt", help = "OPTIONAL OUTPUT. Annotations for hallmark genes in VirSorter-predicted phages / prophages. This file is importable using anvi-import-functions.\n\n")
parser.add_argument("-L","--min-phage-length", type = int, default = 1000, help = "PARAMETER, AFFECTS ALL OUTPUTS. Specify the minimum phage length to report. Default is 1000 bp.")
parser.add_argument("--exclude-cat3", action = "store_true", help = "PARAMETER, AFFECTS ALL OUTPUTS. Excludes all category 3 predictions.")
parser.add_argument("--exclude-prophages", action = "store_true", help = "PARAMETER, AFFECTS ALL OUTPUTS. Exclude all prophage predictions.")
args = parser.parse_args()
arg_dict = vars(args)
#print(arg_dict)

if arg_dict['affi_file'] == None or arg_dict['global_file'] == None or arg_dict['splits_info'] == None or arg_dict['db'] == None:
	print("\n***A required option is missing. Try again.***\n\n")
	parser.print_help()
	sys.exit()
if (arg_dict['anvio_gene_calls'] == None and arg_dict['hallmark_functions'] != None) or (arg_dict['anvio_gene_calls'] != None and arg_dict['hallmark_functions'] == None):
    print("\n***You only specified one of the files required to generate an annotations file for predicted phages.\n")
    parser.print_help()
    sys.exit()


#PART ZERO
#Here I need to read in functional annotations for each hallmark phage cluster
f_dict = collections.defaultdict(dict)
if 'hallmark_functions' in arg_dict:
    f_in = open(arg_dict['hallmark_functions'],"r")
    line = f_in.readline()
    while line != "":
        line = line.strip().split('\t')
        f_dict[line[0]] = line[1]
        line = f_in.readline()
f_in.close()


#PART ONE

contig_set = set()

#Populate contig_set with all of the contigs we're interested in:
virsorter_global = open(arg_dict['global_file'], "r")
line = virsorter_global.readline()
while line != "":
    if line[0] == "#":
        line = virsorter_global.readline()
    else:
        line = line.strip().split(',')
        #print line
        contig = line[0].replace("VIRSorter_","").replace("-circular","")
        contig_set.add(contig)
        line = virsorter_global.readline()
virsorter_global.close()

#First, I need to read in virsorter_affi, the input to Step3 in Virsorter.
#I need this because I need to translate the gene boundaries of each predicted phage
#to start and stop positions (bp) along the contig.

#columns of VIRSorter_affi-contigs.tab:
#     0  | 1   | 2  |  3   |  4   |    5     |  6  |   7  |   8     |   9    | 10   | 11
# gene_id|start|stop|length|strand|affi_phage|score|evalue|category|affi_pfam|score|evalue|

virsorter_affi = open(arg_dict['affi_file'], "r")
affi_dict = {}
affi_dict = collections.defaultdict(dict)
hallmark_dict = collections.defaultdict(dict)

line = virsorter_affi.readline()
while line != "" and line[0] == ">":
    header = line[1:].strip().split('|')
    line = virsorter_affi.readline()
    while line != "" and line[0] != ">":
        var = line.strip().split('|')
        header2 = header[0].replace("VIRSorter_","").replace("-circular","")
        if header2 in contig_set:
            gene = var[0].split('-')
            gene = gene[-1]
            gene_start = var[1]
            gene_stop = var[2]
            hit = var[5]
            cat = var[8]
            evalue = var[7]
            pfam_eval = var[11]
            affi_dict[header2][gene] = {}
            affi_dict[header2][gene]['start'] = gene_start
            affi_dict[header2][gene]['stop'] = gene_stop
            affi_dict[header2][gene]['cat'] = cat
            if cat == str(0) or cat == str(3):
                #This next if statement deals with the fact that VirSorter doesn't count hallmark genes
                #if a hallmark gene has a PFAM hit with a better e-value than the phage cluster hit e-value.
                if pfam_eval == "-" or float(pfam_eval) > float(evalue):
                    affi_dict[header2][gene]['hallmark_function'] = f_dict[hit]
                    hallmark_dict[header2][gene] = {}
                    hallmark_dict[header2][gene]['hallmark_function'] = f_dict[hit]
                    hallmark_dict[header2][gene]['phage_cluster'] = hit
                    hallmark_dict[header2][gene]['start'] = gene_start
                    hallmark_dict[header2][gene]['stop'] = gene_stop
                    hallmark_dict[header2][gene]['cat'] = cat
                    hallmark_dict[header2][gene]['evalue'] = evalue
                    hallmark_dict[header2][gene]['pfam_eval'] = pfam_eval
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
            global_dict[contig]['category'] = str(category)
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
            #print(split_name+" prophage start/stop "+str(start_gene_pos)+" "+str(stop_gene_pos)+" split start/stop "+str(split_start)+" "+str(split_stop))
            if (start_gene_pos <= split_start <= stop_gene_pos or start_gene_pos <= split_stop <= stop_gene_pos) or \
            (split_start <= start_gene_pos <= stop_gene_pos <= split_stop):
                category = "cat"+str(global_dict[split_parent]['category'])+"_prophage"
                global_dict[split_parent]['length'] = stop_gene_pos - start_gene_pos
                phage_length = global_dict[split_parent]['length']
                #print(split_parent)
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


#Part 4a:
#Here I need to read in gene calls from Anvi'o
#File Format:
#gene_callers_id contig  start   stop    direction       partial source  version
#We just read in contigs to gene_dict that were annotated as being all or partially viral.
gene_dict = collections.defaultdict(dict)
if 'anvio_gene_calls' in arg_dict:
    g_in = open(arg_dict['anvio_gene_calls'],'r')
    line = g_in.readline()
    line = g_in.readline()
    while line != "":
        line = line.strip().split('\t')
        if line[1] in contig_set:
            gcid = str(line[0])
            contig = line[1]
            start = line[2]
            stop = line[3]
            direction = line[4]
            gene_dict[contig][gcid] = {}
            gene_dict[contig][gcid]['start'] = start
            gene_dict[contig][gcid]['stop'] = stop
            gene_dict[contig][gcid]['direction'] = direction          
        line = g_in.readline()
    g_in.close
    
    
#Part 4b
#I need to read in the splits file I just generated to get a set of just the contigs we're reporting on here
#rather than integrating the below code into the many cases in part 3 above.
#I'll only export functions for contigs/regions declared interesting in the command line.

#open additional data file, ead in first column splitting on 'split', add to set
if 'anvio_gene_calls' in arg_dict:
    contigs_set = set()
    vir_contigs = open(arg_dict['phage_collection'],"r")
    line = vir_contigs.readline()
    while line != "":
        line = line.strip().split('\t')
        contig = line[0].split('_split_')[0]
        contigs_set.add(contig)
        line = vir_contigs.readline()
    vir_contigs.close()
    
    #write output functions file
    func_out = open(arg_dict['output_functions'],"w")
    func_out.write("gene_callers_id\tsource\taccession\tfunction\te_value\n")

    count = 0
    matched = 0
    for contig in hallmark_dict:
        #if contig in additional data file:... tab the stuff below after turhing this on
        if 'start_gene_pos' in global_dict[contig]:
            prophage_start = global_dict[contig]['start_gene_pos']
            prophage_end = global_dict[contig]['stop_gene_pos']
            for gene in hallmark_dict[contig]:
                h_gene_start = int(hallmark_dict[contig][gene]['start'])
                h_gene_stop = int(hallmark_dict[contig][gene]['stop'])
                if prophage_start <= h_gene_start <= prophage_end or \
                prophage_start <= h_gene_stop <= prophage_end:
                    matched = 0
                    count += 1
                    #print(contig+"\t"+gene)
                    for gcid in gene_dict[contig]:
                        #if h_gene_start == int(gene_dict[contig][gcid]['start'])+1 or h_gene_stop == int(gene_dict[contig][gcid]['stop']):
                        if abs(h_gene_start-int(gene_dict[contig][gcid]['start']))<2 or abs(h_gene_stop-int(gene_dict[contig][gcid]['stop']))<2:
                            #print(contig+",  anvi-gene = "+gcid+", start = "+str(h_gene_start)+", gcid start = "+str(gene_dict[contig][gcid]['start'])+", phage function = "+hallmark_dict[contig][gene]['hallmark_function'])
                            matched = 1
                            pc = hallmark_dict[contig][gene]['phage_cluster']
                            hf = hallmark_dict[contig][gene]['hallmark_function']
                            ev = hallmark_dict[contig][gene]['evalue']
                            #print(hallmark_dict[contig][gene]['cat'])
                            if hallmark_dict[contig][gene]['cat'] == str(3):
                                hf = hf+(" (non-Caudovirales)")
                            func_out.write("%s\t%s\t%s\t%s\t%s\n" % (
                                          str(gcid), "VirSorter (db"+str(arg_dict['db'])+")", pc, hf, str(ev)))
                            break
                    if matched == 0:
                        print("No match for "+contig+"\t"+gene+"\t"+str(h_gene_start)+"\t"+str(h_gene_stop))
        else:
            for gene in hallmark_dict[contig]:
                h_gene_start = int(hallmark_dict[contig][gene]['start'])
                h_gene_stop = int(hallmark_dict[contig][gene]['stop'])
                matched = 0
                count += 1
                #print("    hstart = "+str(h_gene_start)+", hstop = "+str(h_gene_stop))
                for gcid in gene_dict[contig]:
                #    print(gene_dict[contig][gcid]['start']+", "+gene_dict[contig][gcid]['stop'])
                    if abs(h_gene_start-int(gene_dict[contig][gcid]['start']))<2 or abs(h_gene_stop-int(gene_dict[contig][gcid]['stop']))<2:
                    #if h_gene_start == int(gene_dict[contig][gcid]['start'])+1 or h_gene_stop == int(gene_dict[contig][gcid]['stop']):
                        #print("anvi-gene = "+gcid+", "+direction+" , start = "+str(h_gene_start)+", gcid start = "+str(gene_dict[contig][gcid]['start'])+", stop = "+str(h_gene_stop)+", gcid stop = "+str(gene_dict[contig][gcid]['stop'])+", phage function = "+hallmark_dict[contig][gene]['hallmark_function'])
                        matched += 1
                        pc = hallmark_dict[contig][gene]['phage_cluster']
                        #print(hallmark_dict[contig][gene]['cat'])
                        hf = hallmark_dict[contig][gene]['hallmark_function']
                        ev = hallmark_dict[contig][gene]['evalue'] 
                        if hallmark_dict[contig][gene]['cat'] == str(3):
                            hf = hf+(" (non-Caudovirales)")                       
                        func_out.write("%s\t%s\t%s\t%s\t%s\n" % (
                                          str(gcid), "VirSorter (db"+str(arg_dict['db'])+")", pc, hf, str(ev)))
                        break
                if matched == 0:
                    #print("No match for "+contig+"\t"+gene+"\t"+str(h_gene_start)+"\t"+str(h_gene_stop))
                    pass
    func_out.close()