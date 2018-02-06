# VirSorterParser
Display VirSorter annotations in Anvi'o

Metagenome studies are a great way to explore complex communities. Many algorithms and tools have made it possible to reconstruct bacterial, archaeal, and even eukaryotic genome bins for diverse organisms from assembled metagenome sequencing data. However, these algorithms are not well equipped to deal with the most abundant biological entity on earth -- bacteriophages. There are several tools designed to predict viral contigs from metagenome assemblies, and prophages from bacterial genomes.  

This tutorial will walk you through the steps needed to (1) use VirSorter to predict which contigs in your assembly are phages, and (2) visualize these results in [Anvi'o](https://doi.org/10.7717/peerj.1319) using anvi-interactive and anvi-refine. The only thing better than binning with Anvi'o is _phage-aware_ binning with Anvi'o. Here we go!

## How VirSorter works
Virsorter was published in [Roux et al (2015)](https://doi.org/10.7717/peerj.985 "VirSorter: mining viral signal from microbial genomic data"). The source code is housed at https://github.com/simroux/VirSorter. [Figure 1A](https://doi.org/10.7717/peerj.985/fig-1) explains the VirSorter pipeline. Briefly:
- The VirSorter input is a single FASTA file. In our case this will be the same metagenome assembly you used to make your Anvi'o contigs database.
- VirSorter annotates the FASTA file using MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)), and then uses hmmsearch ([Eddy et al, 2011](10.1371/journal.pcbi.1002195))to predict (1) PFAM domains (Version 27, [Finn et al, 2014](https://doi.org/10.1093/nar/gkt1223)) and (2) viral domains using pre-computed, downloadable HMM databases on the annotated genes.
- For each contig, VirSorter uses a sliding window analysis to identify regions of several genes that: (1) contain one or more viral "hallmark" genes (capsid, large terminase, etc.), (2) are enriched for viral domains, (3) have few PFAM domains, (4) have many uncharacterized genes, (5) have many short genes, and (6) have many genes encoded sequentially on the same strand. Overlapping gene regions predicted by each criterion are combined. Depending on which criteria are met, the contig is assigned a category number, where 1 is high confidence that it the contig is a phage, and 3 is possible, but low confidence.
- If a category 1, 2, or 3 prediction on a conting encompasses > 80% of the contig, the whole contig is annotated as being a "phage". If a category 1, 2, or 3 prediction on a contig encompasses <= 80% of the contig, then a subset of the contig is annotated as being a "prophage". For metagenome assemblies, note that a contig annotated as being a "phage" could actually be part of a prophage but isn't called as such simply becasue more than 80% of the genes are "phage-like".

## Getting VirSorter
If you don't want to install anything, you can run VirSorter on CyVerse after signing up for an account here - http://user.cyverse.org/.  

You can also run VirSorter locally [using Docker](https://github.com/simroux/VirSorter), or by manually installing the codebase. A local installation requires you to download additional data files (3.7 GB downloaded, 12 GB when uncompressed) that VirSorter needs to run. This data pack was originally made available on CyVerse, but it has also been uploaded [here](https://doi.org/10.5281/zenodo.1165210 "Data Pack from VirSorter: mining viral signal from microbial genomic data") which will be downloaded later using a simple wget command. 

### Local installation instructions
To manually install the VirSorter codebase locally, do the following. I have only tested it on Linux (Ubuntu and CentOS).
- I prefer to install the VirSorter dependencies using a conda environment. If you have [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html) installed, do the following steps from the directory where you wish to install VirSorter:

```
conda create --name virsorter -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer
git clone https://github.com/simroux/VirSorter.git
cd VirSorter
git checkout bd7a4b7d7d28691bed44caa6d8f07093c882293f
```
  - `perl-parallel-forkmanager` and `diamond=0.9.14` are also required at the end of the `conda create` command for a version of VirSorter I've forked and edited but haven't made available yet.
  - To run VirSorter from any directory, you can make symbolic links to `VirSorter/wrapper_phage_contigs_sorter_iPlant.pl` and `VirSorter/Scripts` and place them in the `bin` folder for your "virsorter" conda environment. An example location of this `bin` folder is `~/miniconda/envs/virsorter/bin`. Substitute this path with the path to the `bin` folder for your newly created "virsorter" environment.
```
ln -s ~/Applications/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl ~/miniconda/envs/virsorter/bin
ln -s ~/Applications/VirSorter/Scripts ~/miniconda/envs/virsorter/bin
```
You'll also need to download MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)). I like to just put this in the virsorter environment's `bin` folder alongside the VirSorter symbolic links.
```
cd ~/miniconda/envs/virsorter/bin
wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf metagene/mga_x86_64.tar.gz
```

### Installing the VirSorter data pack

The virsorter-data archive contains:  
1. PFAM-A and PFAM-B HMM models (Version 27, [Finn et al, 2014](https://doi.org/10.1093/nar/gkt1223))
2. Two HMM databases computed from (1) all phages in RefSeq prior to 2014, and (2) those same phages plus curated phages from several viromes
3. Other files VirSorter needs in order to run

Navigate to a directory where you want the data pack to live, and run the following. It doesn't have to be the same location where you downloaded VirSorter. Run the following commands.  

```
wget https://zenodo.org/record/1165210/files/virsorter-data.tar.gz
md5sum virsorter-data.tar.gz
#m5sum should return d063a8f91038181f0267258249f43345
tar -xvzf virsorter-data.tar.gz
```

The HMM models in `virsorter-data.tar.gz` are built for HMMER version 3.0. At this point, you have two options:
Option A (recommended): Install hmmer as part of the `conda create` command (see below). Once you download and extract `virsorter-data.tar.gz`, you'll need to change to the `virsorter-data` directory (containing the folders PFAM, Phage_gene_catalog, and Phage_gene_catalog_plus_viromes). Next, activate your virsorter virtual environment (or make sure the version of hmmer that VirSorter will be using is in your `PATH`) and run the following command:  

```
for i in */*.hmm; do echo "Converting ${i}..."; hmmconvert ${i} > ${i}.new; mv ${i}.new ${i}; hmmpress -f ${i}; done
```

Option B: Manually Install HMMER 3.0 instead of the latest version through conda. If you're using a virtual environment like conda (as described above), omit "hmmer" from the conda create command and manually install HMMER 3.0 instead(http://hmmer.org/download.html).  
  
You'll need to manually tell VirSorter where you installed the data pack. To do this, you'll need to edit the wrapper script to point to the virsorter-data directory on your computer. Navigate to your VirSorter installation directory, and open the file `wrapper_phage_contigs_sorter_iPlant.pl` in your favorite text editor. On line 46 of `wrapper_phage_contigs_sorter_iPlant.pl`, you'll need to change `'/data';` to the path to the virsorter-data archive you've downloaded and extracted, e.g. `'/path/to/virsorter-data'`. Save and exit!   

## Running VirSorter
Of course, you can always run VirSorter on CyVerse, or using Docker. If you've managed it install it yourself, you can run it that way too. If you installed it yourself, once you have the conda environment created, VirSorter downloaded, and the symbolic links made, running VirSorter is as easy as:  

```
source activate virsorter
wrapper_phage_contigs_sorter_iPlant.pl -f assembly.fasta --db 1 --wdir output_directory --ncpu 4
```

It is _critical_ that VirSorter is run on the _exact same FASTA file_ that was used to generate the Anvi'o contigs database, map and profile reads, etc. This point can't be overstated!  
  
The --db argument can either be 1 or 2.  
If set to 1, then VirSorter will use phage HMMs computed from RefSeq phages published before January 2014.  
If set to 2, then Virsorter will use all of the HMMs from option 1, plus additional HMMs from phage contigs identified in curated virome datasets.  

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
Because VirSorter uses MetaGeneAnnotator ([Noguchi et al, 2006](https://doi.org/10.1093/nar/gkl723)) and doesn't accept custom gene calls, the predicted genes will likely not line up with what Anvi'o has predicted on your contigs database. At this time, the parser does not convert VirSorter gene annotations into functions to import into an Anvi'o contigs database. This isn't as important during the binning process though, as binning focuses on contig splits.  
  
The parser generates an additional data file that can be visualized when running anvi-interactive or anvi-refine. The column headers of the additional data file are as follows:  
`split | phage_name | category | num_genes_in_phage | num_phage_hallmark_genes_in_phage`  
For each VirSorter phage or prophage prediction that spans several splits, all other columns besides "split" are identical across splits. These four metrics are reported on the predicted phage or prophage, not on a given split. For example, if there are 86 genes and 3 splits in the phage contig, `num_genes_in_phage` will report `86` for `split_00001`, `split_00002`, and `split_00003`. The same is true for `num_phage_hallmark_genes_in_phage`.  
  
Reported categories in the `category` column include the following (see [Figure 1B](https://doi.org/10.7717/peerj.985/fig-1):  
- cat1_phage
- cat2_phage
- cat3_phage
- cat1_prophage
- cat2_prophage
- cat3_prophage

For the `phage_name` column, the first phage predicted by VirSorter is named "phage_1" and increments by 1 up to "phage_n". Similarly, the first prophage predicted by VirSorter is named "prophage_1" and incrememts by 1 up to "prophage_n".  
  
## Running the parser
To run the parser you just need the python script `virsorter_to_anvio.py`. Arguments are described at the command-line when passing the `-h` flag. Briefly, the script takes as input the two output files from VirSorter (VIRSorter_affi-contigs.tab and VIRSorter_global_signal.csv) and the splits_basic_info.txt file from Anvi'o. These are required inputs.  If you want to test out the parser, samples of each of the required files are provided in the "sample_files" directory of this repository.  
  
You can control which VirSorter predictions are prepared for importing into Anvi'o. The `--exclude-cat3` flag will skip over all "category 3" predictions. The `--exclude-prophages` flag will skip over all prophages. the `-L` flag can be used to specify the minimum phage length to report in the output files. For example, `-L 5000` means that all phage predictions shorter than 5000 bp will be not reported in the output files. These flags can be used in any combination with each other.  
  
## Output files
The parser, by default, generates two output files. The first is "virsorter_additional_info.txt" which can be imported as additional data into Anvi'o. For Anvi'o version 3, you use this file when running `anvi-interactive` or `anvi-refine` by using the flag `-A virsorter_additional_info.txt`. For higher versions of Anvi'o, you [run the following command](http://merenlab.org/2017/12/11/additional-data-tables/).  
`anvi-import-misc-data virsorter_additional_info.txt -p profile.db --target-data-table items`  
  
The second output file is "virsorter_collection.txt". This file can be imported into Anvi'o using `anvi-import-collection` and specifying "virsorter_collection.txt" as the tab-delimited file to import and will generate a collection where each phage prediction will become a bin containing all of the splits for that phage.  
  
## Happy binning!
Congratulations, you can now enjoy an even better "phage-aware" binning experience in Anvi'o!  
  
If you have any questions or problems, please don't hesitate to contact Bryan or submit an issue on GitHub.  
