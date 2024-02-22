# ACAD aeDNA Workshop Commands Eukaryotes

# 
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
srun --export=ALL --ntasks-per-node 15  --nodes 1 --mem 120G  -t 05:00:00 --pty bash
```

For checking your ressources and jobs running use htop / top
```
htop
```
press q for quit to exit the screen


**OBS remember to deactivate your node when done!!!!!**
```
bash /apps/scripts/cluster_allocation.sh 0
```

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
```

Lets create both environments
```
conda env create -f acad-euks_1.yaml
conda env create -f acad-euks_2.yaml
```
and activate the first environment
```
conda activate acad-euks_1.yaml
```

Lastly, we will install bam-filter in the acad-euks_1 using pip
```
pip install bam-filter
```

For your information we have also installed the following programmes that are not part of the conda packages:
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

This was just for you information, and if you wanted to play around with the full size dataset later. So do not dowload. And as James have been going through adaptor trimming and QC I will not go into much detail about this here. How I did this can be found here https://github.com/miwipe/KapCopenhagen


# we start by


```
conda activate acad-euks_1
```


Remove duplicates (100% identical)
```
time vsearch --fastx_uniques ERR10493277_small-FINAL.fq.gz --fastqout ERR10493277_small-FINAL.vs.fq --minseqlength 30 --strand both
```

# removes low complexity reads dust ranges from 1-4 where 1 being the most stringent.
```
time sga preprocess -m 30 --dust-threshold=1 ERR10493277_small-FINAL.vs.fq  -o ERR10493277_small-FINAL.vs.d1.fq
```


grep 'M_A00706' ERR10493277_small-FINAL.vs.d1.fq

for file in ERR10493277_small-FINAL.vs.d1.fq; do grep 'M_A00706' $file | cut -f2 -d@ > $file.readID.tmp ; done

grep 'M_A00706' ERR10493277_small-FINAL.vs.fq | grep -f readID.tmp -v | cut -f2 -d@ > $file.readID_lowcom.tmp

wc -l ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp

seqtk subseq ERR10493277_small-FINAL.vs.fq ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp



for file in ERR10493277_small-FINAL.vs.d1.fq
do
DB=/shared/mikkelpedersen/acad_test/euks_database/refseq211_small_dedup.fa
echo Mapping $file against $DB
time bowtie2 --threads 15 -k 1000 -x $DB -U $file --no-unal | samtools view -bS - > $file.$(basename $DB).bam
done &> refseq211_small_dedup_mapping.log.txt

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

Extra assignment, for those who finishes first you could consider running the non low complexity filtered file for later comparison, using the same command.

```
Mapping ERR10493277_small-FINAL.vs.fq against /projects/lundbeck/people/npl206/databases/organelles_04042022/refseq211_small_dedup.fa
10046690 reads; of these:
  10046690 (100.00%) were unpaired; of these:
    9861326 (98.15%) aligned 0 times
    14720 (0.15%) aligned exactly 1 time
    170644 (1.70%) aligned >1 times
1.85% overall alignment rate

real	7m1.370s
user	89m31.775s
sys	2m25.342s
```

# filtering and refining alignment output using bam-filter and the metaDMG compressbam function

In this step we clean out the header of the alignment file, to only contain references that have

```
/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam
```
```
  /projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam
	-> compressbam: (compressbam.cpp;Feb 20 2024;09:58:30): '/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/misc/compressbam --threads 4 --input ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam --output ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam'
	-> input: ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.bam; output: ERR10493277_small-FINAL.vs.d1.fq.refseq211_small_dedup.fa.comp.bam; out format: wb; ref: (null); nthreads: 4
	-> Header has now been read. Will now start list of refIDs to use
	->  Now at read:     47800001
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

The alignment file now only have references in the header that also received an alignment, and we now turn our focus

```
time filterBAM reassign --bam MED-2021-20-ver15-2LFQY-210811_S18.bam -t 12 -i 0 -A 92 -m 8G -o MED-2021-20-ver15-2LFQY-210811_S18.reassigned.bam -n 3
for file in *.comp.bam; do time filterBAM reassign --bam $file -t 12 -i 0 -A 92 -m 8G -o $file.reassigned.bam -n 3; done
```


```
filterBAM filter -e 0.6 -m 8G -t 12 -n 3 -A 92 -a 95 -N --bam  MED-2021-20-ver15-2LFQY-210811_S18.reassigned.bam --stats MED-2021-20-ver15-2LFQY-210811_S18.stats.tsv.gz --stats-filtered MED-2021-20-ver15-2LFQY-210811_S18.stats-filtered.tsv.gz --bam-filtered MED-2021-20-ver15-2LFQY-210811_S18.filtered.bam
for file in *.reassigned.bam; do time filterBAM filter -e 0.6 -m 8G -t 12 -n 3 -A 92 -a 95 -N --bam $file --stats $file.stats.tsv.gz --stats-filtered $file.stats-filtered.tsv.gz --bam-filtered $file.filtered.bam ; done


acc2tax=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz
nodes=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/nodes.dmp
names=/projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/names.dmp

for file in *.bam.filtered.bam
do
/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp lca --names $names --nodes $nodes --acc2tax $acc2tax --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 --fix_ncbi 0 --threads 12 --bam $file --out_prefix $file
done 



for file in ERR10493277_small-FINAL.vs.fq.refseq211_small_dedup.fa.fastlocal.bam.comp.bam.reassigned.bam.filtered.bam
do
/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp getdamage --run_mode 1 -n 8 -l 30 -p 30 -o $file.getdamage $file 
done





/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp dfit $file.getdamage.bdamage.gz --names $names--nodes $nodes --threads 12 --lcastat $file.getdamage.stat.gz --showfits 2 --nopt 10 --nbootstrap 20 --doboot 1 --seed 1234 --lib ds --out $file.dfit

