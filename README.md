# VirSorterParser
Display VirSorter annotations in Anvi'o

Metagenome studies are a great way to explore complex communities. Many algorithsms and tools have made it possible to reconstruct bacterial, archaeal, and even eukaryotic genome bins for diverse organisms from assembled metagenome sequencing data. However, these algorithms are not well equipped to deal with the most abundant biological entity on earth -- bacteriophages. There are several tools designed to predict viral contigs from metagenome assemblies, and prophages from bacterial genomes.  This tutorial will walk you through the steps needed to (1) use VirSorter to predict which contigs in your assembly are phages, and (2) visualize these results in Anvi'o using anvi-interactive and anvi-refine. The only thing better than binning with Anvi'o is _phage-aware_ binning with Anvi'o. Here we go!

## How VirSorter works
Virsorter was published in Roux et al (2015) - https://peerj.com/articles/985/. Figure 1 explains the pipeline well. Briefly:
- The VirSorter input is a single FASTA file. In our case this will be the same metagenome assembly you used to make your Anvi'o contigs database.
- VirSorter annotates the FASTA file using MetaGeneAnnotator, and then uses hmmsearch to predict (1) PFAM domains and (2) viral domains using pre-computed, downloadable HMM databases on the annotated genes.
- For each contig, VirSorter uses a sliding window analysis to identify regions of several genes that: (1) contain one or more viral "hallmark" genes (capsid, large terminase, etc.), (2) are enriched for viral domains, (3) have few PFAM domains, (4) have many uncharacterized genes, (5) have many short genes, and (6) have many genes encoded sequentially on the same strand. Overlapping gene regions predicted by each criterion are combined. Depending on which criteria are met, the contig is assigned a category number, where 1 is high confidence that it the contig is a phage, and 3 is possible, but low confidence.
- If a category 1, 2, or 3 prediction on a conting encompasses > 80% of the contig, the whole contig is annotated as being a "phage". If a category 1, 2, or 3 prediction on a contig encompasses <= 80% of the contig, then a subset of the contig is annotated as being a "prophage". For metagenome assemblies, note that a contig annotated as being a "phage" could actually be part of a prophage but isn't called as such simply becasue more than 80% of the genes are "phage-like".

## Using VirSorter
You can run VirSorter on CyVerse after signing up for an account here - http://user.cyverse.org/. You'll need to make an account anyways in order to download the 3.7 GB "virsorter-data" archive.
The virsorter-data archive contains:
1. PFAM-A and PFAM-B HMM models
2. Two HMM databases computed from (1) all phages in RefSeq prior to 2014, and (2) those same phages plus curated phages from several viromes.
3. Other files needed to run VirSorter

To install VirSorter locally, do the following. I have only tested it on Linux (Ubuntu and CentOS).
  - I prefer to install the VirSorter dependencies using a conda environment. If you have Anaconda or Miniconda installed, do the following steps from the directory where you wish to install VirSorter:
  - perl-parallel-forkmanager and diamond=0.9.14 are also required for a version of VirSorter I've forked and edited but haven't made available yet
  
```
conda create --name virsorter -c bioconda mcl=14.137 muscle blast hmmer perl-bioperl perl-file-which
git clone https://github.com/simroux/VirSorter.git
cd VirSorter
git checkout bd7a4b7d7d28691bed44caa6d8f07093c882293f
```
  - To run VirSorter from any directory, you can make symbolic links to VirSorter/wrapper_phage_contigs_sorter_iPlant.pl and VirSorter/Scripts and place them in the "bin" folder for your "virsorter" conda environment.
```
ln -s ~/Applications/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl ~/miniconda/envs/virsorter/bin
ln -s ~/Applications/VirSorter/Scripts ~/miniconda/envs/virsorter/bin
```
You'll also need to download MetaGeneAnnotator. I like to just put this in the virsorter environment's bin folder alongside the VirSorter symbolic links.
```
cd ~/miniconda/envs/virsorter/bin
wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf metagene/mga_x86_64.tar.gz
```

After you download the "virsorter-data" archive from CyVerse, you'll need to edit the wrapper script to point to the virsorter-data directory on your computer. On line 46 of `wrapper_phage_contigs_sorter_iPlant.pl`, you'll need to change `'/data';` to the path to the virsorter-data archive you've downloaded and extracted, e.g. `'/path/to/virsorter-data'`.

## Running VirSorter
Once you have the conda environment created, VirSorter downloaded, and the symbolic links created, running VirSorter is as easy as:
```
source activate virsorter
wrapper_phage_contigs_sorter_iPlant.pl -f assembly.fasta --db 1 --wdir output_directory --ncpu 4
```
It is _critical_ that VirSorter is run on the _exact same assembly FASTA file_ that was used to generate the Anvi'o contigs database, map and profile reads, etc.

## VirSorter outputs
VirSorter writes several outputs to a working directory that you specify when run VirSorter. The files we need for importing VirSorter annotations into Anvi'o include:
- VIRSorter_global-phage-signal.csv
  - This file contains one line for each phage prediction. Often, this results in one line per contig, though very large contigs might have two or more prophage predictions.
- Metric_files/VIRSorter_affi-contigs.tab
  - This file contains one line per gene and includes any PFAM or phage domain annotations. The lines for all genes for a given contig are preceded by a line containing ">Contig_name".

## Anvi'o files needed
Once you have your Anvi'o contigs database generated, run the following command to generate the final file we need to import VirSorter annotations into Anvi'o.

`anvi-export-table --table splits_basic_info CONTIGS.db`

This will generate a file called "splits_basic_info.txt" file.

## Rationale of the parser
Because VirSorter uses MetaGeneAnnotator and doesn't accept custom gene calls, the predicted genes will likely not line up with what Anvi'o has predicted on your contigs database. At this time, the parser does not convert VirSorter gene annotations into functions to import into an Anvi'o contigs database. This isn't as important during the binning process though, as binning focuses on contig splits.

The parser generates an additional data file that can be visualized when running anvi-interactive or anvi-refine. The column headers of the additional data file are as follows:
`split	phage_name	category	num_genes_in_phage	num_phage_hallmark_genes_in_phage`
For each VirSorter phage or prophage prediction that spans several splits, all other columns besides "split" are identical across splits. These four metrics are reported on the predicted phage or prophage, not on a given split. For example, if there are 86 genes and 3 splits in the phage contig, `num_genes_in_phage` will report `86` for `split_00001`, `split_00002`, and `split_00003`. The same is true for `num_phage_hallmark_genes_in_phage`.

Reported categories in the `category` column include:
cat1_phage
cat2_phage
cat3_phage
cat1_prophage
cat2_prophage
cat3_prophage

For the `phage_name` column, the first phage predicted by VirSorter is named "phage_1" and increments by 1 up to "phage_n". Similarly, the first prophage predicted by VirSorter is named "prophage_1" and incrememts by 1 up to "prophage_n".

## Running the parser
