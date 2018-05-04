
<p align="center">
  <a>
    <img height="350" src="logo/PmR.png">
  </a>
</p>

Pimp my RAD (PmR) is a pipeline providing different tools to guide the treatment of ddRAD sequencing data. In addition to automatically treat the sequences using Stacks functions (Catchen et al. 2013), the pipeline allow for data cleaning, data visualization and help users to choose pertinent thresholds during loci reconstruction and genetic dataset export. 


### 1 . Installation
--------------------

PmR is tested on Unix environment and requires:
+ bbduk (http://jgi.doe.gov/data-and-tools/bbtools/)
+ R (http://www.r-project.org/)
+ Stacks (version 1.46, http://catchenlab.life.illinois.edu/stacks/)
+ Needs awk, java, gcc and g++


### 2. Input files
------------------


PmR requires the paired-end ddRAD sequences files (R1 and R2), a file of barcodes/individuals matching, a file of parameters and a file of adapters sequences. All these files have to be directly placed in the **DATA/** directory.


#### 2.1 - Parameters files (params.txt)


The following section describes the params files required for PmR. Users have to modify the params.txt file provided before running PmR. Default values are set for filtering and assembly steps. As the step 4 of PmR performs a range of parameter setting to optimise the assembly of loci, some values have to be modified after the step 4 (see section 3.). 

`less params.txt`

```
# Global parameters :
TOOLS=~/Path_2_TOOLS/tools.sh                  ## path to file with tools aliases
DATA=~/Path_2_DATA/                            ## path to directory containing the input library reads and the barcode file (see next section to specifications)
RES=~/Path_2_RES/                              ## path to directory to write output
THREADS=15                                     ## Number of threads which will be used
VERBOSE=1                                      ## Set verbose to TRUE (1) or FASE (0)

# Step 1 
ADAP=~/Path_2_PmR/ressources/adapters.fa       ## list of adapters for filtering reads
MIN_LENGTH=60                                  ## minimum length of reads kept after adapters trimming

# Step 2
ENZYME_1=pstI                                  ## first restriction enzyme used
ENZYME_2=mspI                                  ## second restriction enzyme used
MAX_LENGTH=90                                  ## truncate final read length to this value 

# Step 3

# Step 4 
COVERAGE_LOC_MIN=5                             ## Minimum depth of coverage to create a stack (m)
MISMATCH_LOC_IND_START=0                       ## Starting value of the Maximum distance (in nucleotides) allowed between stacks (M)
MISMATCH_LOC_IND_END=7                         ## Ending value of the Maximum distance (in nucleotides) allowed between stacks (M)
NB_INDIV_M=10                                  ## Number of individuals used for the game


# Step 5
M_CHOSEN=6                                     ## Number of mismatches allowed between sample loci when building the catalog (M)

# Step 6 
MISMATCH_CATALOG_MAX=8                         ## Number of mismatches allowed between a read and a loci in the catalog (n) 

# Step 7 

# Step 8 
POP_INFOS=${DATA}/Populations_table.txt        ## population file
NB_POP=3                                       ## the locus must be present in at least X populations
PROP_POP=0.10                                  ## minimum percentage of individuals in a population required to process a locus for that population
MAF=0.005                                      ## minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5)
```

#### 2.2 - Dependencies (tools.txt)

The access path of all dependencies required by PmR must be supplied in the tools.sh file, using following command: 

`cat tools.sh`

```
#!/bin/bash

BBDUK=/mypathToBbduk/bbduk.sh
PROCESS_RADTAG=/mypathToprocess/process_radtags
USTACKS=/mypathToustacks/ustacks
CSTACKS=/mypathTocstacks/cstacks
SSTACKS=/mypathTosstacks/sstacks
POPULATIONS=/mypathTopopulation/populations
```

#### 2.3 - Raw data 


All Raw Data must be in the **DATA/** repertory.
this repertory must contain :
+ two sequences files
+ a barcode file
+ a population file (optional)

The sequences files are in .fastq.gz format and should be named :
+ Basename_Lib_R1.fastq.gz
+ Basename_Lib_R2.fastq.gz

where <Basename_Lib> is the name of the library.

`ls DATA/`

```
Library_test_R1.fastq.gz
Library_test_R2.fastq.gz
Library_test.barcode
```

`less DATA/Library_test_R1.fastq.gz`

```
@K00268:48:H7YLJBBXX:2:1101:3285:998 1:N:0
ATTTATATATAAATCATCGAGCGAGGAGAGCTAGAGAGC
+
#AAFFA77<F7FFJJJJJJJJJJJJJJJJJJJJJJJJFJJAFJJJJJJ
```

The barcodes file is a table with barcodes and samples matches. Each line should have one name and one barcode, separated by a tab, as in the following example:

This file have to be called Basename_Lib.barcode, where <Basename_Lib> is the name of the library.

`head DATA/Library_test.barcode`

```
CGATA <tab> Sample1
CTGTT <tab> Sample2
```

The population map file links individuals to different population groups for summary statistics computing such as π, FIS and FST. Each line should have one name and one barcode, separated by a tab, as in the following example:

`head Populations_table.txt`
```
Sample1Name<tab>pop01
Sample2Name<tab>pop01
Sample3Name<tab>pop02
```

### 3. Pipeline description
---------------------------

PmR mostly used the Stacks functions to treat the ddRAD sequencing data. We strongly recommend to look at the Stacks documentation before using the pipeline. 

#### 3.1 - Adapters trimming (Step 1)

If sequencing adapters are present in the sequences, they are removed by bbduk. The list of illumina adapters (containing all indexes) is provided with PmR in the **ressources/** directory (adapters.fa) and could be directly be used. All reads that were present in the input file will also be present in the output file, even reads that were trimmed entirely (because the adapter was found in the very beginning).

A directory named **STEP_1_FITLER_READS/**, is created at this step containing the clean data which will be used by PmR.

#### 3.2 - Data demultiplexing and reads filtering (Step 2)

The process_radtags function from Stacks is run on each library to demultiplex and filter the reads on their quality with a sliding window of 15% of the read length and a Phred score limit to 10. By default the quality score of reads is encoded in phred33. Reads with low quality or uncalled bases were removed. This step requires to specify a barcode file linking barcodes to samples and the couple of enzymes used during the lab protocol.  

Currently supported enzymes by Stacks include:
'aciI', 'ageI', 'aluI', 'apeKI', 'apoI', 'aseI', 'bamHI', 'bfaI', 'bgIII', 'bspDI', 'bstYI', 'claI', 'ddeI', 'dpnII', 'eaeI', 'ecoRI', 'ecoRV', 'ecoT22I', 'hindIII', 'kpnI', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nheI', 'nlaIII', 'notI', 'nsiI', 'pstI', 'rsaI', 'sacI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'speI', 'sphI', 'taqI', 'xbaI', or 'xhoI'.

This function output four files per barcode in a directory **STEP_2_DEMULTIPLEX/**, two for the single-end reads and two for the paired-end reads.

`ls STEP_2_DEMULTIPLEX/`

```
sample_ACTCG.1.fq
sample_ACTCG.rem.1.fq
sample_ACTCG.2.fq
sample_ACTCG.rem.2.fq
```

PmR will create a `inds.tsv` file containing all demultiplexed samples after the step 2.


#### 3.3 - Samples files preparation for loci reconstruction (Step 3)

The results of the cleaning and the demultiplexing steps have to be checked before the reconstruction of the loci. For example, according to the distribution graph of the reads number, samples with too few reads should be identified.

PmR will perform a de novo assembly according to a list of samples, allowing the analysis of multiple libraries into a single run. PmR will create a inds.tsv file containing all demultiplexed samples after the step 2.
If the same sample is sequenced on different runs (i.e. differents libraries), PmR will concatenate sequences of the sample to be analyses as a single sample. Samples identified with too few reads must be removed from the list (by removing the line), and will not be assembled further. Users have to modify the name of the sample in the last file column according the following example:

`head inds.tsv`
```
160817_SNK268_B_L002_GWM-817-4<tab>AOS_7_Cla<tab> AOS_7_Cla
160817_SNK268_B_L002_GWM-817-4<tab>AOS_7_Ens<tab>AOS_7_Ens
160817_SNK268_B_L002_GWM-817-4<tab>AOS_8 <tab>AOS_8
160817_SNK268_B_L002_GWM-817-6<tab>AOS_8<tab>AOS_8
160817_SNK268_B_L002_GWM-817-6<tab>AOS_12<tab>AOS_12
```

In this case, the samples AOS_7_Cla and AOS_7_Ens are a same AOS_7 individual from the same sequencing run (160817_SNK268_B_L002_GWM-817-4) and AOS_8 had been sequenced twice in 160817_SNK268_B_L002_GWM-817-4 and 160817_SNK268_B_L002_GWM-817-6 libraries. AOS_12 is not processed because of too few reads (identify at the step 1). To perform assembling step, we thus keep only the AOS_7 and AOS_8 names and we remove AOS_12 from the list.

`nano inds.tsv`
```
160817_SNK268_B_L002_GWM-817-4<tab>AOS_7_Cla<tab> AOS_7
160817_SNK268_B_L002_GWM-817-4<tab>AOS_7_Ens<tab>AOS_7
160817_SNK268_B_L002_GWM-817-4<tab>AOS_8<tab>AOS_8
160817_SNK268_B_L002_GWM-817-6<tab>AOS_8<tab>AOS_8
```

All the output files will be in the **STEP_3_PREPARE_SAMPLE/** directory.

#### 3.4 - Range of M threshold for loci reconstruction (Step 4)

A de novo assembly is performed on aligned reads for each sample from the provided list using the ustacks function. The minimum coverage to create a stack of identical reads (-m) was fixed (params file) and a range of values of different nucleotides was tested to merge two different stacks into one polymorphic locus (-M). Highly-repetitive stacks and over merged tags were dropped using both the “Removal algorithm” (-r) and the “Deleveraging algorithm” (-d). PmR will create a STEP_4_ASSEMBLE_LOCI directory containing a directory for each M value supplied in the range. 

Each of M_<Value> directory contains three files for each sample provided in the assembly file: 
+ sampleXXX.tags.tsv: containing each loci assembled
+ sampleXXX.snps.tsv: containing each loci assembled  
+ sampleXXX.alleles.tsv: Haplotypes/alleles recorded from each locus

The choice of the optimal M value depends on the distribution of the number of polymorphic stacks according to the range of M values (given in output) and is arbitrary fixed by the user. The retained value should correspond for example to the threshold M value for which Stacks will not merge anymore loci in spite of the increase of the distance allowed between loci. 

#### 3.5 - Loci reconstruction for each individuals (Step 5)

Considering the chosen M value, a reconstruction of loci is made for each individuals separatle by the ustacks function.
Output files will be in the **STEP_5_7_RUN_STACKS/** directory. 


#### 3.6 - Build the catalog of loci (Step 6)

A catalog of the loci from all the individuals’ loci data sets is built using the cstacks function by producing consensus alleles for each locus. A maximum value of mismatches allowed for considering two individual tags as the same locus and to merge them in the catalogue (-n) must be supplied in the params.txt  file.
Output files will be in the **STEP_5_7_RUN_STACKS/** directory. 

#### 3.7 - Match of individual loci to the catalog (Step 7)

The sets of loci constructed within each individual reads pool (step 5) were then searched against the newly produced catalog (step 6) using the sstacks function. Two files are created batch_X.catalog.tags.tsv corresponding to the catalog and batch_X.matches.tsv containing the matches of the catalog. 
Output files will be in the **STEP_5_7_RUN_STACKS/** directory. 

#### 3.8 - Genetic dataset export (Step 8)

Populations genetics outputs and statistic are performed with the populations function of Stacks. PmR writes output file in FASTA, VCF, GENEPOP, STRUCTURE and PHYLIP format to produce SNP datasets with individuals genotype which could be used for further analysis. Fasta format contains full sequence for each unique haplotype, from each sample locus in FASTA format, regardless of plausibility. The PHYLIP format contains nucleotides that are fixed-within, and variant among populations in PHYLIP format for phylogenetic tree construction. This alignment could be direct use for phylogenetic inference. The populations function will also calculate population genetics statistics such as expected/observed heterozygosity, π, and FIS at each nucleotide position. 

A population file is required, and options must be edited in the params file (see input files section). NB_POP and PROP_POP will use the population file provided to select loci in the final dataset according the given value. For example, considering NB_POP=3 and PROP_POP=0.75, populations will create a final dataset with loci recovered in at least three populations and 75% of individuals in each population. 

For further phylogenetic analysis, if a different population for each sample is supplied (in case of each sample correspond to unique species, see input files section), NB_POP=3 will fix the condition that a locus must be present in three species at least. In this case, PROP_POP will not be considered in the populations function. 
Output files will be in the STEP_8_FINAL_OUTPUT/ directory. 
4. Running process


### 4. Running PmR
------------------

PmR could be called step by step or all the steps in one time. Recommendations about steps are given in the previous description (section 3). After edition of the tools.sh and params.txt files, the first step of quality analysis and demultiplexing should be called as following: 

`./PmR.sh -s 1 -p params.txt`

To run the complete pipeline:

`./PmR.sh -s all -p params.txt`


### 5. References
--------------------


Catchen, Julian, Paul A. Hohenlohe, Susan Bassham, Angel Amores, and William A. Cresko. 2013. “Stacks: An Analysis Tool Set for Population Genomics.” Molecular Ecology 22 (11): 3124–40.
