# ACAD aeDNA Workshop Commands Eukaryotes

Login into the cloud server as the other days using ssh

```
ssh user@acadworkshop.uoa.cloud
```
Check that there are available resources up and running 
```
sinfo
```

As we will be going into detail on each command/programme lets open an interactive bash shell on the slurm
```
srun --export=ALL --ntasks-per-node 4 --nodes 1 --mem 40G  -t 02:00:00 --pty bash
```

For checking your ressources and jobs running use `htop` / `top`
```
htop
```
press q for quit to exit the screen

First download the data, taxonomy, database and basic scripts for this tutorial and unzip all files (except the accc2taxid file which needs to stay zipped.

```
wget -O 2m_tutorial_data.tar.gz https://sid.erda.dk/share_redirect/cM9JHwqbPz
tar -xzvf 2m_tutorial_data.tar.gz 
cd 2m_tutorial_data

for file in *fq.gz; do gunzip $file & done
for file in database/refseq211_small_dedup.fa.*; do gunzip $file & done
for file in small_taxonomy/n*; do gunzip $file & done
```

## Make github repo and clone it to your home directory in the cloud
GitHub is a web-based platform used for version control and collaboration on software development projects. It provides a variety of features for developers and teams to work together efficiently on coding projects.

Log in to GitHub: Go to https://github.com and log in to your GitHub account. If you don't have one, you'll need to sign up first.

**Consider whether you want it to be private or public (It is possible to apply for an educational licence [for free] and then you have the private option), you can change it from private to public, when private only you and people you have invited can see the repo**

To work with the repository locally on your computer, you need to clone it. To do this, click on the green "Code" button on your forked repository's page, copy the HTTPS or SSH URL, then use Git to clone the repository to your local machine.

```
git clone https://github.com/your-username/repository.git
```

Add a Remote (Optional): By default, Git will add a remote named "origin" that points to your forked repository on GitHub. If you want to keep track of the original repository that you forked from, you can add a remote with a different name.

Make Changes and Push: Now that you have a local copy of your repository, you can start writing the code, commit them, and push them back to your online repository on GitHub.

Add a header to the first line (Header) in the README.md file using a text editor (vim, vi, nano). Now add the changes, commit writing a message and push it to your copy of the repo.
```
git status
git add .
git commit -m "Name added to header"
git push
```
You can now navigate to your repo in your browser and check that the change has been made. (Refresh the page if necessary)

Now let us list the commits
```
git log
git log --oneline --decorate --graph
```
That should print something similar to this:
* 8425734 (HEAD -> main, origin/main, origin/HEAD) added github tutorial and fix metaDMG commands
* f8b9f70 Update README.md
* 106148d Update README.md
* 2081ee6 Initial commit

You can roll back to previous versions using 'git checkout', maybe if you made a mistake or just need to run a previous version. 

```
git checkout xxXXXxx 
```

and you can roll back to the main (most recent version)

```
git checkout main  
```

These are the most common commands and would be useful perhaps for each project. A bit of good advice is to keep your repo updated (git push) but only when you have made working / significant changes! Make one for each of your projects, this is good practice and you need to have this for your eventual publication anyway. Trust me you will save time by doing this on the fly instead of waiting to publication is close.

The last cool thing about git is that it allows for multiple users to collaborate on the same code, while still keeping track of all versions. And as a Centre or group, you can have an organization in GitHub like https://github.com/GeoGenetics where you have public and private repos that are shared.

You can also fork a repo e.g. make your copy of the repo, and develop it for yourself or help on the repo development

To fork a Git repository means to create a copy of it in your own GitHub account (or any other Git hosting service). Here are the general steps to fork a Git repository:

Find the Repository to Fork: Navigate to the repository you want to fork. You can do this by searching for the repository in the GitHub search bar or by accessing it through a direct link or finding it on the list of repositories on Mikkels GitHub profile.

Fork the Repository: Once you're on the repository's page, you'll see a button labelled "Fork" at the top right corner of the page. Click on this button. GitHub will then create a copy of the repository under your GitHub account.

```
git remote add upstream https://github.com/original-owner/original-repository.git
```
Syncing with the Original Repository (Optional): If you want to keep your forked repository up-to-date with changes made to the original repository, you can fetch the changes from the original repository and merge them into your local repository.

```
git fetch upstream
```
From this point you handle it as your repo, using git add, git commit and git push. Now if you forked it to fix an issue for the original owner and users, you can create/open a pull request (or merge request) on Git Hub. This pull request has to include information about the changes made, the purpose of the changes, and any other relevant details.

The original developers, collaborators, or team members then review the pull request. They can leave comments, ask questions, suggest changes, or approve the changes. Once the pull request has been reviewed and approved, the changes can be merged into the main branch of the repository. Alternatively, if the changes are not ready or if there are issues, the pull request can be closed without merging.

Pull requests are a fundamental part of the collaborative development process in Git-based projects. They enable code review, collaboration, and integration of changes in a controlled and organized manner.

## Setting up your conda environments and other dependencies
Create two (acad-euks_1 and acad-euks_2) environment files and hereafter the corresponding conda environment.
First, use a text editor in the terminal, it could be vi/vim/nano, copy and paste the content in the box below and save the file as acad-euks_1.yaml.

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
  - r-dplyr
  - r-tidyr
  - r-readr
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
- r-base
- r-ggplot2
- r-dplyr
- r-tidyr
- r-readr
```

Lets create both environments (NOTE! that each will take a bit of time ~5 min)!
```
conda env create -f acad-euks_1.yaml
conda env create -f acad-euks_2.yaml

wget https://raw.githubusercontent.com/genomewalker/bam-filter/master/environment.yml
conda env create -f environment.yml
conda activate bam-filter
pip install bam-filter
conda install r-base r-tidyverse r-ggplot2 r-readr
conda deactivate bam-filter
```

For your information, we have also installed the following programmes that are not part of the `conda` packages:
```
metaDMG-cpp https://github.com/metaDMG-dev/metaDMG-cpp
ngsLCA https://github.com/miwipe/ngsLCA
```

The dataset we are going to play with is downscaled from raw data, and 10 mill. reads were extracted from a QC'ed, dereplicated, low complexity trimmed fastq file. It has further been spiced up with a few surprises. Raw data can be downloaded here. But please don
```
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/077/ERR10493277/ERR10493277_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/077/ERR10493277/ERR10493277_2.fastq.gz
```

This was just for your information, and if you wanted to play around with the full-size dataset later. So do not download. And as James has been going through adaptor trimming and QC I will not go into much detail about this here. How I did this can be found here https://github.com/miwipe/KapCopenhagen

For this tutorial, I created not just a downscaled sample, described above, but also a small database and taxonomy for a rapid walk through each step. It should allow you to run through all in half a day, including the installation of environment tools, and plottings. The data, database and corresponding taxonomy can be downloaded here.
```

```
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
time sga preprocess -m 30 --dust-threshold=1 ERR10493277_small-FINAL.vs.fq  -o ERR10493277_small-FINAL.vs.ds1.fq
```

It is good practice to double-check that your output is good, and what consequences the filters have had for the data. Let us check what sequences were filtered out.
```
# prints the readIDs of each sequence that parsed filters
grep 'M_A00706' ERR10493277_small-FINAL.vs.ds1.fq

# we can make this into a text file with all readIDs that are parsed
grep 'M_A00706' ERR10493277_small-FINAL.vs.ds1.fq | cut -f2 -d@ > ERR10493277_small-FINAL.vs.ds1.readID.tmp

# we can then inverse grep these readIDs in the original fastq file, to get the readIDs of the filtered sequences
grep 'M_A00706' ERR10493277_small-FINAL.vs.fq | grep -f ERR10493277_small-FINAL.vs.ds1.readID.tmp -v | cut -f2 -d@ > ERR10493277_small-FINAL.vs.ds1.readID_lowcom.tmp

# how many reads did we filter out, and does it correspond to the stdout sga printed?
wc -l ERR10493277_small-FINAL.vs.ds1.readID_lowcom.tmp

# we can use seqtk to grep out all sequences that were categorized as low complex, what do you see?
seqtk subseq ERR10493277_small-FINAL.vs.fq ERR10493277_small-FINAL.vs.ds1.readID_lowcom.tmp
```

First validation step extract and plot the read lengths of your QCed fasta file.


# Readlength distribution from fastq files
```
cat ERR10493277_small-FINAL.vs.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ERR10493277_small-FINAL.vs.read_length.txt &

cat ERR10493277_small-FINAL.vs.ds1.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ERR10493277_small-FINAL.vs.ds1.read_length.txt &
```

Lets install the missing libraries in your R installation. By first opening the R console 

```
R
```

Define a variable for package names you miss, add or remove what you need. 
```
packages_to_install <- c("tidyverse", "readr", "dplyr")
```
Install packages
```
install.packages(packages_to_install)
```

To exit the R console 
```
q()
```

Now let us run a script that takes one or more read_length.txt files in the folder you are standing in and plots them using ggplot. Make sure that the R conda environment is activated or that there is a system-wide R installation

```
/data/scripts/readlengthPLOT_fastq.sh
```

Let us download the pdf and have a look at it. Are there differences between the fqs read length distributions?
Use scp or any other ways you normally get data downloaded to your local machine if you are not running it on a desktop machine or like. For example.
```
scp -r user@name.of.server:/path/2/files/filename.pdf . 
```

How does the read-length distribution look? Remember this is before taxonomic classification and ancient DNA filtering, so this can be a composite of ancient and living organisms. Typically it would be. 

## The Mapping 
Next, what can `bowtie2` do and what options have we set? Look at the command below and check the options and settings.

```
bowtie2 --help
```

Let us map the reads to our database (the database was generated by indexing the multi-reference fasta file refseq211_small_dedup.fa with `bowtie2-build -@ 4 refseq211_small_dedup.fa`). Depending on the file size, it will run for a short or longer time. We skip this step, due to limited time. A bowtie2 database consists of 6 subfiles ending .bt2l. Let us first check they are there.

```
ls -lh  database/refseq211_small_dedup.fa*
```

Map all QC'ed reads against a database, let's create a slurm bash script. The script essentially aligns sequencing reads (fastq files) to a reference genome database using `Bowtie2` and then converts the output to BAM format using `samtools`. OBS if you are not on a slurm system, you can simply copy the commands below and run them individually or make a loop in a bash script. save the below as bowtie2.sh, to submit the jobs to the queue. 

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --job-name=bowtie2_mapping
#SBATCH --mem=30G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
#SBATCH --output=bowtie2_mapping_%A_%a.out

time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.ds1.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.bam
time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.ds1.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.bam 
time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 1 -L 21 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.ds1.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.bam 

time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.bam 
time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.bam 
time bowtie2 --threads 12 -k 5000 -x database/refseq211_small_dedup.fa -D 15 -R 2 -N 1 -L 21 -i S,1,1.15 -U ERR10493277_small-FINAL.vs.fq --no-unal | samtools view -bS - > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.bam
```
Then run this and your jobs are submitted (go get coffee, maybe lunch and hopefully your jobs are done by then, however, it depends on how busy the system is and your slurm setup and more.) 

```
sbatch bowtie2.sh
```

### Bowtie options explained 
#!/bin/bash: This line is called a `shebang` and indicates that the script should be interpreted by the `Bash shell`.
#SBATCH directives: These lines are specific to Slurm and are used to specify parameters for the job submission to the HPC cluster. Here's what each directive means:
--nodes=1: Requests 1 compute node.
--time=2:00:00: Requests a maximum of 2 hours of runtime for the job.
--job-name=bowtie2_Bianca: Specifies the name of the job.
--mem=30G: Requests 30 gigabytes of memory.
--output=bowtie2.out: Redirects the standard output of the job to the file bowtie2.out.
bowtie2: This is a command-line tool used for aligning sequencing reads to a reference genome using the Burrows-Wheeler Transform algorithm. By default the alignments 
--threads 4: Specifies the number of CPU threads to use for the alignment process.
k 1000: Reports up to 1000 valid alignments for each read. OBS the number of alignment extensions bowtie2 performs will significantly increase cpu time, but might not yield significantly more alignments! The optimal setting will depend on your database and data. 
-x /shared/database/refseq211_small_dedup.fa: Specifies the path to the Bowtie2 index of the reference genome.
-D 15: is used to set the "maximum edit distance" e.g. max number of mismatches to the reference sequence
-R 2: sets the maximum number of times Bowtie2 will attempt to extend an alignment in a row for each read. 
-N 0: maximum number of mismatches to the seed
-L 22: Length of seed. 
-i S,1,1.15: Sets various alignment parameters for Bowtie2.

The -i parameter in Bowtie2 sets the interval between which the quality of the alignment and the quality of the read sequence are rescaled for alignment score computation. This parameter influences how Bowtie2 computes alignment scores based on the quality values of the read sequences.
-i <func>: This parameter sets the interval between which the quality of the alignment and the quality of the read sequence are rescaled for alignment score computation. The default value is "5,50,0.50", which means that the quality values of the read are rescaled such that quality value 5 maps to 50 with a mapping function of 0.50.
The -i parameter allows for adjusting how Bowtie2 interprets the quality values of the read sequences during alignment. By rescaling the quality values, Bowtie2 can more accurately compute the alignment scores and determine the quality of the alignments.

So, with -i S,1,1.15, quality values of the read sequence are rescaled linearly such that quality value 1 is mapped to 1.15 for alignment score computation.
S: Indicates a "scale" function, which means the quality values of the read sequence are linearly scaled.
1: Indicates that quality value 1 is mapped.
1.15: Indicates that quality value 1 is mapped to 1.15.
-U ERR10493277_small-FINAL.vs.ds1.fq: Specifies the input file containing unpaired reads.
--no-unal: Suppresses SAM records for reads that failed to align.
|: This is a pipe, used to redirect the output of bowtie2 to the input of `samtools`.
samtools view -bS - > ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.bam: Converts SAM output from `bowtie2` into BAM format and writes it to the file ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.bam.

When bowtie2 is done mapping it prints to stdout and should look similar to this

```
Mapping ERR10493277_small-FINAL.vs.ds1.fq against /shared/database/refseq211_small_dedup.fa
10045794 reads; of these:
  10045794 (100.00%) were unpaired; of these:
    9860963 (98.16%) aligned 0 times
    14678 (0.15%) aligned exactly 1 time
    170153 (1.69%) aligned >1 times
1.84% overall alignment rate

real	13m20.957s
user	68m48.312s
sys	1m7.282s
```

Extra assignment! For those who finish fast. Consider running the non-filtered file for later comparison, using the same command. Or changing the `bowtie2` option settings could be the seed length and or the mismatches allowed on the seed. 

# Filtering and refining alignment output using bam-filter and the metaDMG compressbam function

If you are on a slurm system, you can likely ask for a bash shell using the below command (depending on the system setup, contact your local admin). If not, please neglect the following command!
```
srun --export=ALL --ntasks-per-node 8 --nodes 1 --mem 10G  -t 02:00:00 --pty bash
```
You can also just copy the commands below. 

In this step we clean out the header of the alignment file, to only contain references that have reads aligning to. please double check that the path to the metaDMG-cpp/misc/compressbam is correct, it depends on where it is installed. The misc folder is inside the metaDMG-cpp misc folder

``` 
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=compressbam
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
#SBATCH --output=compress_header_out_%A_%a.out

/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.bam --output ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.log.txt
/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.bam --output ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.log.txt
/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.bam --output ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.log.txt
/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.bam --output ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.bam &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.log.txt 
/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.bam --output ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.bam &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.log.txt
/metaDMG-cpp/misc/compressbam --threads 12 --input ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.bam --output ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.bam &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.log.txt
```

The stdout should look similar to this once `compressbam` is done. (This tool is valuable when handling databases with large amounts of reference genomes. Especially when the header of the bam files exceeds 2Gb sizes then samtools cannot handle the header as a bam file format, which is the reason we developed this tool). The stdout I saved into the log files ending *.comp.log.txt  
```
/apps/software/metaDMG-cpp/misc/compressbam  --threads 4 --input ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.bam --output ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.bam
	-> compressbam: (compressbam.cpp;Feb 20 2024;09:58:30): '/apps/software/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.comp.bam'
	-> input: ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam; output: ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.fa.comp.bam; out format: wb; ref: (null); nthreads: 4
	-> Header has now been read. Will now start list of refIDs to use
	-> Now at read:     47800001
	-> Done reading list of refids: to keep 57260
	-> Number of alignments parsed: 47825234
	-> Done writing new header as text
	-> header info before: nref: 204143 bytesize: 6209816
	-> header info  after: nref: 57260 bytesize: 1674655
	-> header info reduction nref: 0.280490 bytesize: 0.269679
	-> Done writing new header as binary
	->  Now at read:     47800001 	-> Done writing file: 'ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.bam'
	-> [ALL done] cpu-time used =  126.91 sec walltime used =  48.00 sec
```

The alignment file now only have references in the header that also received an alignment, and we now turn our focus to the alignments in the "header compressed" bam file. As we have specified the -k 1000 option in bowtie2, reads can potentially have up to 1000 alignments, lets double check that we did not saturate the number of alignments allowed and thereby not allowing the read an equal chance to map against all alignments. This is just simple text gymnastics. We cut out the readID in each alignments sort these and count how many times the read occurs. if = 1000 it is saturated.

```
samtools view ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.fa.comp.bam | cut -f1 | sort | uniq -c | sort -k1 -n > counts_of_alignments.txt

tail -500 counts_of_alignments.txt 
```

**Free assignment:** Perhaps grep for one of the reads with 5000 hits, and blast the sequence, what is it hitting maybe a conserved region? Or why does it attract so many alignments? 

As `ngsLCA` (also the version embedded in `metaDMG`), does not look specifically at each alignment to determine which is closest, we have developed a small algorithm to run over the bam and keep alignments that are closest to reference inside `bam-filter`. This discards alignments that have poorer alignment stats than the set values in options.

Now as the conda version of `bam-filter` that we installed was an old version we just need to install the latest, but in a new environment. `pip` has the latest version. 

```
conda activate bam-filter
```

Now let us run the reassignment of each read hits to the reference.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=reassign
#SBATCH --mem=30G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=filterbam_reaasign_out_%A_%a.out

filterBAM reassign --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.log.txt 
filterBAM reassign --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.log.txt 
filterBAM reassign --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.bam.log.txt

filterBAM reassign --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.log.txt
filterBAM reassign --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.log.txt
filterBAM reassign --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.bam -t 4 -i 0 -A 92 -m 8G -o ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam -n 10 &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam.log.txt
```

Now we can do the same check as above, how does it look? 
```
samtools view ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam | cut -f1 | sort | uniq -c | sort -k1 -n
```

Now let us use bamfilter on the reassigned alignments to make the final filtering, bam-filter now creates statistics for both the non-filtered and the filtered file. Go explore! 


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --job-name=filterBAM
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=filterBAM_filter_out_%A_%a.out

filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam.log.txt 
filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam.log.txt 
filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam.log.txt 

filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam.log.txt 
filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam.log.txt 
filterBAM filter -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam --stats ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam.stats.tsv.gz --stats-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam.stats-filtered.tsv.gz --bam-filtered ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam &> ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam.log.txt 
```

Have a look at the log files and familiarize yourself with this general output.

```
cat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam.log.txt 
```
Okay, now go explore the stats, first, let us start with the non-filtered file. Start R, and install packages needed (we did this for the first environment but not in acad-euks_2). You can also simply download the filtered and non-filtered stats files *.stats-filtered.tsv.gz and *.stats.tsv.gz 

```
R
install.packages("readr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr)

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the table
data <- readr::read_tsv("ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.bam.stats-filtered.tsv.gz")

# Histogram of a read_length_mean 
plot <- ggplot(data, aes(x = read_length_mean)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "read_length_mean", x = "Length in Bp", y = "Frequency")

# Save the plot as a PDF file
ggsave("read_length_mean_histogram_plot.pdf", plot, width = 8, height = 6)
```
You can try and change the x variable too, look into the https://github.com/genomewalker/bam-filter documentation for variable explanations. 

```
# Create a point plot
point_plot <- ggplot(data, aes(x = breadth, y = breadth_exp_ratio)) +
  geom_point() + # Add points
  labs(title = "Breadth vs. Breadth Exp Ratio", x = "Breadth", y = "Breadth Exp Ratio")

# Save the plot as a PDF file
ggsave("point_plot.pdf", point_plot, width = 8, height = 6)
```
or lastly, make a correlogram using GGally

```
install.packages("GGally")
library(GGally)

data_filtered <- data %>% filter(n_reads >= 500) %>% select(n_reads, read_ani_median, coverage_mean_trunc, exp_breadth, cov_evenness, coverage_covered_mean, norm_gini, breadth, gini, c_v, norm_entropy, edit_distances)

# Create the pairwise scatterplot matrix
 cor_plot <- ggpairs(data_filtered, 
        progress = FALSE,
        aes(alpha = 0.65))

ggsave("correlogram_plot.pdf", cor_plot, width = 10, height = 10)
```

Download plot into a folder on your local machine, here I download it to my Download folder
```
noglob scp -r $USER@path2.cloud:/path2project/*.pdf .
```

Try and play around with the different variables and their relationship with each other, select two plots with variables that we can discuss in plenum. Consider which statistic would be good to use for filtering, what values and why?

```
# Create another plot
point_plot2 <- ggplot(data_filtered , aes(x = cov_evenness, y = norm_entropy, color = n_reads > 700, size = n_reads)) +
  geom_point() + # Add points
  geom_text(data = subset(data, n_reads > 700), aes(label = reference), vjust = -0.5) + # Add labels for n_reads > 700
  scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "black")) + # Set color for points
  scale_y_continuous(limits = c(0, 1.01)) + 
  scale_x_continuous(limits = c(0, 1.01)) + 
  labs(title = "cov_evenness vs. norm_entropy", x = "cov_evenness", y = "norm_entropy") +
  theme(text = element_text(size = 17))

# Save the plot as a PDF file
ggsave("cov_evennessVSnorm_entropy_plot.pdf", point_plot2, width = 8, height = 6)

```

## Taxonomic classification (ngsLCA) and DNA damage estimation (metaDMG)
**While `metaDMG-cpp lca` is running through the alignment file, it also collects mismatch information for all reads and their alignments**

We start by running the lca analysis, OBS please note that while going through the alignment file, `metaDMG-ccp lca` also spits out look-up files (such as mismatch matrices, read lengths, gc content and more) for the later dfit.


Now we move to the taxonomic classification of the reads (ngsLCA embedded in metaDMG-cpp), if on a slurm ask for a bash shell. 
```
srun --export=ALL --ntasks-per-node 3 --nodes 1 --mem 8G  -t 02:00:00 --pty bash
```

Run the LCA. 

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=metadmgLCA
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=metadmgLCA_out_%A_%a.out

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp lca --names small_taxonomy/names.dmp --nodes small_taxonomy/nodes.dmp --acc2tax small_taxonomy/small_accession2taxid.txt.gz --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 4 --bam ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bam --out_prefix ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered
```

Now familiarize yourself with the log output that `metaDMG lca`, importantly if there are accession numbers

hint https://www.ncbi.nlm.nih.gov/
another important look-up tool https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

are there other ways to go back and find out what is happening?

hint the database and the taxonomy used.

Cool tools to know! 
https://bioinf.shenwei.me/csvtk/
https://github.com/TabViewer/tabview

Tabview dependencies were incompatible with the acad-euks_1 environment but were in the acad-euks_2. Let us activate this

`Tabview`
Navigation:
**Arrow keys**: Move the cursor up, down, left, or right.
**Page Up and Page Down**: Scroll the view up or down by one page.
**Home and End**: Move to the beginning or end of the data.
Sorting:
**S or s** for numbers >1 also works for text strings. 
**A or a** for decimal numbers 1-0 
Searching:
**/**: Start searching and anter your query.
**Enter**: Find the next occurrence of the search term.
**Escape**: Exit search mode.

```
conda activate acad-euks_2
```

Quick and dirty damage calculations (for the impatient)
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=metadmgDfit
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=metadmgDfit_out_%A_%a.out

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered

```

For full stats (for the patient).
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=metadmgDfit
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=metadmgDfit_out_%A_%a.out

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered

metaDMG-cpp/metaDMG-cpp dfit ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234  --lib ds --out ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered

```

Aggregate (sum) all results up the taxonomic tree from nodes to root.
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=metadmgAgg
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --output=metadmgAgg_out_%A_%a.out

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

metaDMG-cpp/metaDMG-cpp aggregate ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz -lcastat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.stat.gz --names small_taxonomy/names.dmp --nodes  small_taxonomy/nodes.dmp

```

Paste the resulting lca stats and damage stats for a final metaDMG out.

```
zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.combined_metaDMG_output.txt

zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.combined_metaDMG_output.txt

zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.combined_metaDMG_output.txt

zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.combined_metaDMG_output.txt

zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.combined_metaDMG_output.txt

zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.bdamage.gz.stat.gz | paste - <(zcat ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.dfit.gz) > ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.combined_metaDMG_output.txt
```
Lastly, we concatenate the files together to one longtable while printing the corresponding filenames to each line. 

```
header_file="ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.combined_metaDMG_output.txt"

# Get the header
header=$(head -n 1 "$header_file")

# Add the header to the concatenated file
echo -e "filename\t$header" > concatenated_metaDMGdata.txt

for file in ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.combined_metaDMG_output.txt \
            ERR10493277_small-FINAL.vs.ds1.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.combined_metaDMG_output.txt \
            ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.comp.reassign.filtered.combined_metaDMG_output.txt \
            ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L22.N1.comp.reassign.filtered.combined_metaDMG_output.txt \
            ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.L21.N1.comp.reassign.filtered.combined_metaDMG_output.txt
do
    tail -n +2 "$file" | while read -r line; do
        echo -e "$file\t$line" >> concatenated_metaDMGdata.txt
    done
done
```

Now let us explore the output, perhaps using tabview or less -S or ??
```
tabview concatenated_data.txt
```

Make plots perhaps initially with one sample
```

# Read in the table
data <- read.table("ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.L22.comp.reassign.filtered.combined_metaDMG_output.txt", 
                   header = TRUE, # If the file has a header row
                   sep = "\t",    # Specify the delimiter used in the file (tab-separated)
                   stringsAsFactors = FALSE)  # Ensure character columns are read as strings

# Display the first few rows of the data
head(data)

# Set variable for filtering, maybe start without and then set variable depending on your data?
# Minimum amount of damage (filtered for)
A_Min = 0.20
#Lambda Likelihood Ratio
Significance = 3.0
# Minimum reads for parsing taxa
MinRead = 200
# Minimum mean read length
MeanLength = 35
###### subsetting the table using grepl and filter, parameters you need to set and possibly add more but here we look at plants only
df1 <- data %>% filter(A > A_Min, nreads > MinRead, mean_rlen >= MeanLength, Zfit > Significance,  grepl("Viridiplantae",taxa_path), grepl("\\bgenus\\b", rank))


# Create the ggplot bar plot
bar_plot <- ggplot(df1, aes(x = name, y = nreads)) +
  geom_bar(stat = "identity", fill = "lightblue") +  # Create bars with nreads values
  labs(title = "Bar Plot of nreads by genus", x = "genus", y = "number of reads") +  # Add labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees

# Save the plot as a PDF file
ggsave("bar_plot.pdf", plot = bar_plot, width = 8, height = 6)



all_data <- read.table("concatenated_data.txt", 
                   header = TRUE,
                   sep = "\t", 
                   stringsAsFactors = FALSE) 


# Create the ggplot tile plot
tile_plot <- ggplot(all_data, aes(x = nreads, y = name)) +
  geom_tile(aes(fill = nreads)) +  # Add tiles with nreads as fill color
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Set gradient fill colors
  labs(title = "", x = "Number of reads", y = "Taxa (genus)")  # Add labels

# Save the plot as a PDF file
ggsave("tile_plot.pdf", plot = tile_plot, width = 8, height = 6)




```



