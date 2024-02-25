# ACAD aeDNA Workshop Commands Eukaryotes

Login into the cloud server as the other days using ssh

```
ssh user@acadworkshop.uoa.cloud
```

Then lets activate your node
```
bash /apps/scripts/cluster_allocation.sh 1
```

As we will be going into detail on each command/programme lets open an interactive bash shell on the slurm
```
srun --export=ALL --ntasks-per-node 15  --nodes 1 --mem 120G  -t 02:00:00 --pty bash
```

For checking your ressources and jobs running use `htop` / `top`
```
htop
```
press q for quit to exit the screen


**OBS remember to deactivate your node when done!!!!!**
```
bash /apps/scripts/cluster_allocation.sh 0
```

## Fork github repo and clone it to your home directory

To fork a Git repository means to create a copy of it in your own GitHub account (or any other Git hosting service). Here are the general steps to fork a Git repository:

Log in to GitHub: Go to https://github.com and log in to your GitHub account. If you don't have one, you'll need to sign up first.

Find the Repository to Fork : Navigate to the repository you want to fork. You can do this by searching for the repository in the GitHub search bar or by accessing it through a direct link.

Fork the Repository: Once you're on the repository's page, you'll see a button labeled "Fork" at the top right corner of the page. Click on this button. GitHub will then create a copy of the repository under your GitHub account.

Clone Your Forked Repository: After forking, you'll have your own copy of the repository on your GitHub account. To work with the repository locally on your computer, you need to clone it. To do this, click on the green "Code" button on your forked repository's page, copy the HTTPS or SSH URL, then use Git to clone the repository to your local machine.

Example:

```
git clone https://github.com/your-username/forked-repository.git
```
Add a Remote (Optional): By default, Git will add a remote named "origin" that points to your forked repository on GitHub. If you want to keep track of the original repository that you forked from, you can add a remote with a different name.

csharp
Copy code
git remote add upstream https://github.com/original-owner/original-repository.git
Syncing with the Original Repository (Optional): If you want to keep your forked repository up-to-date with changes made to the original repository, you can fetch the changes from the original repository and merge them into your local repository.

sql
Copy code
git fetch upstream
git merge upstream/main
Make Changes and Push: Now that you have your own forked repository, you can make changes to the code, commit them, and push them back to your forked repository on GitHub.


## Setting up your conda environments and other dependencies
Create two (acad-euks_1 and acad-euks_2) environment files and here after the corresponding conda environment.
First use a text editor in the terminal, it could be vi/vim/nano, copy and paste the content in the box below and save the file as acad-euks_1.yaml.

**acad-euks_1**
```
name: acad-euks_1
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.8,<=3.9
  - Cython>=0.29.24
  - pip
  - pip:
      - pyrle>=0.0.31
      - pyranges>=0.0.112
      - kneed>=0.8.1
      - matplotlib>=3.6.0
  - samtools
  - bowtie2
  - r-base
  - r-ggplot2
  - fastp
  - vsearch
  - sga
  - seqtk
  - fasta-splitter
  - datamash
  - seqkit
```

Now repeat the same here but calling the files acad-euks_2.yaml.

**acad-euks_2**
```
name: acad-euks_2
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
- beast2=2.6.3=hf1b8bbb_0
- mafft=7.475=h516909a_0
- raxml-ng=1.0.1=h7447c1b_0
- tabview
- angsd
- samtools
- epa-ng
- gappa
- seqkit
```

Lets create both environments (NOTE! that each will take a bit of time ~5 min)!
```
conda env create -f acad-euks_1.yaml
conda env create -f acad-euks_2.yaml
```
and then activate the first environment
```
conda activate acad-euks_1.yaml
```

Lastly, we will install bam-filter in the acad-euks_1 using pip
```
pip install bam-filter
```

For your information we have also installed the following programmes that are not part of the `conda` packages:
```
pathPhynder https://github.com/ruidlpm/pathPhynder
Bamcov https://github.com/fbreitwieser/bamcov
metaDMG-cpp https://github.com/metaDMG-dev/metaDMG-cpp
ngsLCA https://github.com/miwipe/ngsLCA
Euka https://github.com/grenaud/vgan
```

The dataset we are going to play with is downscaled from raw data, and 10 mill. reads were extracted from a QCed, dereplicated, low complexity trimmed fastq file. It has further been spiced up with a few surprises. Raw data can be downloaded here. But please don
```
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/077/ERR10493277/ERR10493277_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/077/ERR10493277/ERR10493277_2.fastq.gz
```

This was just for you information, and if you wanted to play around with the full size dataset later. So do not download. And as James have been going through adaptor trimming and QC I will not go into much detail about this here. How I did this can be found here https://github.com/miwipe/KapCopenhagen


# Lets get some data analysis done

If not already on the node log in as described above, then activate the acad-euks_1 conda environment
```
conda activate acad-euks_1
```

**Quality control of trimmed and merged sequences**
When handling large data and mapping against large reference genome collections, it can be important to remove duplicates, to save cpu and run time. For this, I use `vsearch` https://github.com/torognes/vsearch which is a fast tool that screens for 100% identical sequences (most likely caused by PCR duplication). You can use `vsearch --help` to familiarize yourself with its options.
```
time vsearch --fastx_uniques ERR10493277_small-FINAL.fq.gz --fastqout ERR10493277_small-FINAL.vs.fq --minseqlength 30 --strand both
```

Another important aspect is to clean out low complex sequences again several tools can do this, I use sga / bbduk in which dust ranges from 1-4 where 1 is the most stringent.
```
time sga preprocess -m 30 --dust-threshold=1 ERR10493277_small-FINAL.vs.fq  -o ERR10493277_small-FINAL.vs.d1.fq
```

It is good practice to double-check that your output is good, and what consequences the filters have had for the data. Lets check what sequences was filtered out.
```
# prints the readIDs of each sequence that parsed filters
grep 'M_A00706' ERR10493277_small-FINAL.vs.d1.fq

# we can make this into a text file with all readIDs that parsed
for file in ERR10493277_small-FINAL.vs.d1.fq; do grep 'M_A00706' $file | cut -f2 -d@ > $file.readID.tmp ; done

# we can then inverse grep these readIDs in the original fastq file, to get the readIDs of the filtered sequences
grep 'M_A00706' ERR10493277_small-FINAL.vs.fq | grep -f readID.tmp -v | cut -f2 -d@ > $file.readID_lowcom.tmp

# how many reads did we filter out, and does it correspond to the stdout sga printed?
wc -l ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp

# we can use seqtk to grep out all sequences that where categorized as low complex, what do you see?
seqtk subseq ERR10493277_small-FINAL.vs.fq ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp
```

First validation step extract and plot the read lengths of your QCed fasta file.


# extract readlength distribution from fastq files
```
for file in *.vs.fq
do
cat $file | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.read_length.txt &
done
```


Rename your readlength distribution files, by copying, pasting and saving this bash script (like we did with the environment files) as rename_readlength_files.sh

```
#!/bin/bash

for file in *.read_length.txt
do
    if [[ -f $file ]]; then
        new_file="$(echo $file | cut -d'.' -f1).txt"
        mv "$file" "$new_file"
        echo "Renamed $file to $new_file"
    fi
done
```

Now change the permissions of the file for it to be executable using `chmod +x`

```
chmod +x rename_readlength_files.sh
```

Now lets run a script that takes one or more read_length.txt files and plots them using ggplot. Make sure that the R conda environment is activated or there is a system wide R installation

```
/shared/mikkelpedersen/acad_test/tutorials/acad_workshop2024/scripts/readlengthPLOT_fastq.sh
```

Lets download the pdf and have a look at it. Is there differences between the fqs read length distributions?

Next, what is it `bowtie2` can do and what options have we set? Look at the command below and check the options and settings.

```
bowtie2 --help
```

Lets map the reads to our database (database was generated by indexing the multi reference fasta file refseq211_small_dedup.fa with `bowtie2-build -@ 12 refseq211_small_dedup.fa`). Depending on the size of the file it will run for short and longer time. We skip this step, due to limited time. A bowtie2 database consists of 6 subfiles ending .bt2l. Lets first check they are there.

```
ls -lh  /shared/mikkelpedersen/acad_test/euks_database/refseq211_small_dedup.fa*
```

Map reads against database
```
DB=/shared/mikkelpedersen/acad_test/euks_database/refseq211_small_dedup.fa
time bowtie2 --threads 15 -k 1000 -x $DB -U ERR10493277_small-FINAL.vs.d1.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.bam

```

The stdout should look like this, time might vary a bit
```
Mapping ERR10493277_small-FINAL.vs.d1.fq against /shared/mikkelpedersen/acad_test/euks_database/refseq211_small_dedup.fa
10045794 reads; of these:
  10045794 (100.00%) were unpaired; of these:
    9860963 (98.16%) aligned 0 times
    14678 (0.15%) aligned exactly 1 time
    170153 (1.69%) aligned >1 times
1.84% overall alignment rate

real	9m18.101s
user	90m15.170s
sys	1m17.688s
```

Extra assignment! Aor those who finishes fast. Consider running the non-filtered file for later comparison, using the same command. Or changing the `bowtie2` option settings


# filtering and refining alignment output using bam-filter and the metaDMG compressbam function

In this step we clean out the header of the alignment file, to only contain references that have

```
time /shared/mikkelpedersen/acad_test/euks_programmes/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam
```
The stdout should look similar to this once `compressbam` is done. (This tool is valuable when handling databases with large amounts of reference genomes. Especially when the header of the bam files exceeds 2Gb sizes then samtools cannot handle the header as a bam file format, which is the reason we developed this tool)
```
  /shared/mikkelpedersen/acad_test/euks_programmes/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam
	-> compressbam: (compressbam.cpp;Feb 20 2024;09:58:30): '/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam'
	-> input: ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam; output: ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam; out format: wb; ref: (null); nthreads: 4
	-> Header has now been read. Will now start list of refIDs to use
	-> Now at read:     47800001
	-> Done reading list of refids: to keep 57260
	-> Number of alignments parsed: 47825234
	-> Done writing new header as text
	-> header info before: nref: 204143 bytesize: 6209816
	-> header info  after: nref: 57260 bytesize: 1674655
	-> header info reduction nref: 0.280490 bytesize: 0.269679
	-> Done writing new header as binary
	->  Now at read:     47800001 	-> Done writing file: 'ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam'
	-> [ALL done] cpu-time used =  161.29 sec walltime used =  84.00 sec
```

The alignment file now only have references in the header that also received an alignment, and we now turn our focus to the alignments in the "header compressed" bam file. As we have specified the -k 1000 option in bowtie2, reads can potentially have up to 1000 alignments, lets double check that we did not saturate the number of alignments allowed and thereby not allowing the read an equal chance to map against all alignments. This is just simple text gymnastics. We cut out the readID in each alignments sort these and count how many times the read occurs. if = 1000 it is saturated.

```
samtools view ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam | cut -f1 | sort | uniq -c | sort -k1 -n
```

As ngsLCA (also the version embedded in metaDMG), does not look specifically at each alignment to determine which is closest, we have developed a small algorithm to run over the bam and keep alignments that are closest to reference inside `bam-filter`. This discard alignments that have poorer alignment stats than the set values in options.
```
time filterBAM reassign --bam ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam -t 12 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.reassigned.bam -n 3

for file in *.comp.bam; do time filterBAM reassign --bam $file -t 12 -i 0 -A 92 -m 8G -o $file.reassigned.bam -n 3; done
```


```
filterBAM filter -e 0.6 -m 8G -t 12 -n 3 -A 92 -a 95 -N --bam  MED-2021-20-ver15-2LFQY-210811_S18.reassigned.bam --stats MED-2021-20-ver15-2LFQY-210811_S18.stats.tsv.gz --stats-filtered MED-2021-20-ver15-2LFQY-210811_S18.stats-filtered.tsv.gz --bam-filtered MED-2021-20-ver15-2LFQY-210811_S18.filtered.bam
for file in *.reassigned.bam; do time filterBAM filter -e 0.6 -m 8G -t 12 -n 3 -A 92 -a 95 -N --bam $file --stats $file.stats.tsv.gz --stats-filtered $file.stats-filtered.tsv.gz --bam-filtered $file.filtered.bam ; done
```

## Taxonomic classification (ngsLCA) and DNA damage estimation (metaDMG)
**While `metaDMG lca` is running through the alignment file, it also collects mismatch information for all reads and their alignments**


```
#acc2tax=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz
nodes=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/nodes.dmp
names=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/names.dmp
acc2tax=/projects/lundbeck/people/npl206/TMP/small_accession2taxid.txt.gz

for file in ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.sort.bam
do
time /projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp lca --names $names --nodes $nodes --acc2tax $acc2tax --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 12 --bam $file --out_prefix $file
done

```

Now familiarize yourself with the log output that `metaDMG lca`, importantly if there are accession numbers

hint https://www.ncbi.nlm.nih.gov/
another important look-up tool https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

are there other ways we can go back and find out what is going on?

hint the database and the taxonomy used.

Now
https://bioinf.shenwei.me/csvtk/
https://github.com/TabViewer/tabview

Tabview dependencies was not compatible with the acad-euks_1 environment, but was in the acad-euks_2. Lets activate this
```
conda activate acad-euks_2
```


quick and dirty version (for the impatient)
```
for file in *.sort.bam *filtered.bam
do
time /projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp dfit $file.bdamage.gz --names $names --nodes $nodes --showfits 2  --lib ds --out $file
done
```

Do full stats (for the patient)
```
for file in *.sort.bam *filtered.bam
do
time /projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp dfit $file.bdamage.gz --names $names --nodes $nodes --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234 --lib ds --out $file
done
```

Aggregate (sum) all results up the taxonomic tree from nodes to root.
```
for file in *.sort.bam *filtered.bam
do
time /projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp aggregate $file.bdamage.gz --lcastat $file.stat.gz --names $names --nodes $nodes
done
```

Paste the resulting lca stats and damage stats for a final metaDMG out.
```
zcat $file.bdamage.gz.stat.gz | paste - <(zcat $file.dfit.gz) > $file.combined_metaDMG_output.txt

```

Now lets explore the output, perhaps using tabview
```
tabview $file.combined_metaDMG_output.txt
```
