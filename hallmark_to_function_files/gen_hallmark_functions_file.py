#!/usr/bin/env python

import collections
import re
import argparse

"""
This script parses the file \"Phage_Clusters_current.tab\". For each hallmark gene in this file, 
(category 0 or 3), it exports a two-column table where each line represents a phage cluster
and the cluster\'s predicted function:

Phage cluster      Function
Phage_cluster_1    large terminase

Writes output file \"hallmark_functions.txt\".
"""

parser = argparse.ArgumentParser(description="Extracts functions from Phage_Clusters_current.tab into a mapping file used by virsorter_to_anvio.py parser.")
parser.add_argument('clusters', help="Specify the file \"Phage_Clusters_current.tab\" from one of the database directories inside your virsorter-data pack.")
args = parser.parse_args()
arg_dict = vars(args)

clusts = open(arg_dict["clusters"],"r")
out = open("hallmark_functions.txt","w")

words_to_include = "spike|terminase|terminase large|large terminase|portal|coat|major capsid|tail tube|tail sheath|tail protein|spike|virion structural|tail tube|minor capsid|encapsidation|major tail|capsid protein|tubular tail|tail tubular|baseplate hub|tail component|phage capsid"
words_to_exclude = "unnamed|endopeptidase|stabilizer|assembly|tail fiber|minor tail protein|tailspike|tail spike|chaperonin|endolysin|tape measure|tape-measure|transglycosylase|protease|stabilization|transcriptional|receptor|completion|chaperone|lysozyme|hypothetical"

clust_functions = collections.defaultdict(dict)

line = clusts.readline()
while line != "":
    line = line.strip().split('|')
    if int(line[1]) == 0 or int(line[1]) == 3:
        functs = line[2].split(';')
        keep = ""
        for item in functs:
            item = item.split(':')
            if bool(re.search(words_to_include, item[0])) == True and len(item[0]) >= 5: #and bool(re.search(words_to_exclude, item[0])) == False:
                keep = item[0]
                break
        if keep == "":
            for item in functs:
                item = item.split(':')
                if bool(re.search(words_to_exclude, item[0])) == False and len(item[0]) >= 5:
                    keep = item[0]
                    break
        if keep[0] == " ":
            keep = keep[1:]
        out.write((line[0]+"\t"+keep+"\n"))
    line = clusts.readline()

clusts.close()
out.close()