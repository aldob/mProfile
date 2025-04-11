# mProfile
## Installation
    pip install mprofile-mut
<br>
## Overview
A package for processing targeted sequencing (amplicon) data for high-resolution mutation profiling. 

The callMUT tool converts mutation data from the Samtools mpileup format, which contains raw mutation calls for every read, into a lightweight, analysis ready mprofile: a table cotaining the percentage rate of every mutation type at every nucleotide. 

The TransloCapture tool maps all possible crossover events between different targets in a targeted sequencing experiment. 

The StructureMap tool maps structural variants in whole genome data using paired-end alignment files.



<br>
<br>
<br>
<br>
<br>



## callMUT
The default output of Samtools Mpileup is very large and cumbersome to parse, a problem that increases with alignment file size. 

Mprofiles are comparitively very lightweight and do not expand with increasing alignment file size. They provide per-nucloetide mutation rates (% of reads) in a table that makes plotting and further analysis simple and fast! 

Also able to calculate the differential between a treated(-i) and a control(-c) sample (either in mpileup or mprofile format).

#### Arguments
    callMUT -i input.mpileup -o output.mprofile
    
    Required arguments:
      --input INPUT, -i INPUT
                            Input mpileup or mprofile to process.
      --output OUTPUT, -o OUTPUT
                            Output mprofile.

    Additional arguments:
      --indelcut [INDELCUT], -ic [INDELCUT]
                            Minimum rate for an indel's sequence to be reported,
                            default=1, If 'NA', no cutoff will be applied.
      --control CONTROL, -c CONTROL
                            mpileup/mprofile to normalise to (e.g. untreated).
      --preproc, -pp        Specifies input files are mprofiles, not mpileups
                            (requires --control to be set).
      --quiet QUIET, -q QUIET
                            Removes all messages.
      --help -h HELP
                        show this help message and exit

    Example: 
      callMUT -i treated.mpileup -c untreated.mpileup -o treated.mprofile -ic 0.001

<br>

#### Example .mprofile table
IMPORTANT: callMUT requires that samtools mpileup is run using the -aa command and with the same bed file for every condition being compared! This ensures that all mpileups have the same length and report the same nucleotides.

callMUT accepts the standard output format of samtools mpileup and the default output is a tab seperated table that reports every nucleotide in the mpileup input, the readcount and all the mutation rates at that nucleotide. 

With -ic NA set, all indels found in the input file will be reported in the Common.Indels column. <br>
These are reported as a comma seperated list of all indels for that nucleotide position in the format "indel sequence:indel rate". <br>
Example: -5ATTAT:0.012 means the following five nucleotides ATTAT are deleted at a rate of 0.012% 

Chromosome|Coordinate|Ref.Base|Readcount|A.Mutations|T.Mutations|G.Mutations|C.Mutations|Transitions|Transversions|Total.SNVs|Insertions|Deletions|Common.Indels
----------|----------|--------|---------|------|------|------|------|-----------|-------------|----------|----------|---------|-------------
chr1	|88992916	|A	|1486103	|0	|0.012355217	|0.02932642	|0.008889362	|0.02932642	|0.02124458	|0.050571	|0.020751696	|0.020646069|+1T:0.02032100326111,-9TCGCGGAAT:6.74002509985e-05,
chr1	|88992917	|T	|1471955	|0.002030475	|0	|0.010923032	|0.021511631	|0.021511631	|0.012953507	|0.034465138	|0.014427304	|0.018301666|
chr1	|88992918	|C	|1480451	|0.003216938	|0.000509182	|0.028957037	|0	|0.000509182	|0.032173975	|0.032683158	|0.000578539	|0.005032715|-3CTG:0.005010032611



<br>
<br>
<br>
<br>
<br>



## TransloCapture
In the multiplex PCR used for targeted sequencing, crossover events (translocations) between the targets will also be amplified and sequenced. TransloCapture takes unprocessed fastq files and identifies reads that are crossover events. The tool works by identifying the primers used for amplification of the targets at either end of the reads.

Requires a csv file(--primers, -p) containing three columns "target name", "forward primer sequence", "reverse primer sequence" for each amplicon target. The output is a .csv file containing a matrix of all possible crossovers and their frequencies (% of reads) and can also output the crossover reads to a new fastq file.

Accepts either single-read(SR) or paired-end(PE) data, however PE is highly reccommended for this analysis! SR will significantly reduce the accruacy.

#### Arguments
    TransloCapture -1 input_read1.fastq -2 input_read2.fastq -o output.csv

    Required arguments:
      --input INPUT, -i INPUT
                              Input fastq file for SR sequencing
      --read1 READ1, -1 READ1
                              Fastq read 1 from PE sequencing
      --read2 READ2, -2 READ2
                              Fastq read 2 from PE sequencing
      --output OUTPUT, -o OUTPUT
                              Output file to write to, format is csv
      --primers PRIMERS, -p PRIMERS
                              A 3 column .csv file of the name, foward primer
                              sequence and reverse primer sequence (reverse
                              complement) for each site to be analysed.

    Additional arguments:
      --control CONTROL, -c CONTROL
                              The fastq you want to normalise to (e.g. untreated).
                              If unspecified, will not normalise.
      --control1 CONTROL1, -c1 CONTROL1
                              Read 1 of the fastq you want to normalise to (e.g.
                              untreated). If unspecified, will not normalise.
      --control2 CONTROL2, -c2 CONTROL2
                              Read 2 of the fastq you want to normalise to (e.g.
                              untreated). If unspecified, will not normalise.
      --preproc, -pp          If specified, --input (-i) and --control (-c) must be
                              already quantified TransloCapture matrices. Output
                              will be a new matrix that is the differential of
                              input-control.
      --fastqout FASTQOUT, -fo FASTQOUT
                              Fastq file to write translocated sequences to. If
                              unspecified, will not write
      --fastqout1 FASTQOUT1, -fo1 FASTQOUT1
                              Fastq file to write read1 of translocated sequences
                              to. If unspecified, will not write.
      --fastqout2 FASTQOUT2, -fo2 FASTQOUT2
                              Fastq file to write read2 of translocated sequences
                              to. If unspecified, will not write.
      --quiet, -q             Removes all messages.

      --help -h HELP
                              show this help message and exit

    Example: 
      TransloCapture -1 treated_R1.fastq -2 treated_R2.fastq -o treated_translomap.csv -p target_primers.csv
      
<br>

#### Example Primer table
First, here is an example primer table that is used to create the example translocation table.<br>
The table should have no header, but the columns are target name, forward primer sequence, reverse primer sequence.

.          |.                         |.                           |
-----------|--------------------------|----------------------------|
CTRL_Target|TACGCTACGACTAGCAGCTATCGACT|ATCGCATCTAGACTGATCACGATCTACG|
Target1|TACGACACTACGACTATCGA|ATCGCTACGACGATAGCTTCGT|
Target2|GCTACGCTACAGCGACTACTA|GAGCTACGACATCACCGCTGCATCA|
Target3|GACTACGACTACGACTACG|CGCGACTACGACGCTACGCGC|
Target4|CGCGACTACGCGATCACGACTACTAGC|CGATCACGACTACGACGCATCG|

<br>

#### Example Crossover table
The x-axis are the sites recognised in Read1 (or at the start if SR) and the y-axis are the sites recognised in Read2 (or at the end if SR).<br>
Cells where x=y are the un-translocated targets and are therefore represented by NA.

Site|CTRL_Target|Target1|Target2|Target3|Target4
----|-----|----|----|----|----|
CTRL_Target|NA|0|0|0|0|
Target1|0|NA|0.01673|0.0002336|0.23451|
Target2|0|0.010564|NA|0.01568|0.006168|
Target3|0|0.000005169|0.000168168|NA|0.000079841|
Target4|0|0.0016546|0.078965|0.0165169|NA|



<br>
<br>
<br>
<br>
<br>



## StructureMap
By identifying the position of each mate in a read pair, StructureMap determines if this fragment has been translocated (different chromosomes), had a large deletion (insert size greater than threshold) or other variants.

It requires a sorted alignment file that can be SAM, BAM or CRAM. 

There is an option to output a file that contains all the insert sizes for every succesfully aligned read (excludes translocations etc).

DEPENDENCY: pysam, unlike other mprofile tools StructureMap has a dependency on pysam to read the alignment files

#### Arguments
    StructureMap_v3.py [-h] [--input INPUT] [--output OUTPUT]
                              [--insert_max [INSERT_MAX]]
                              [--fragment_sizes FRAGMENT_SIZES]

    StructureMap -i input.sam -o output.sprofile

    optional arguments:
      -h, --help            show this help message and exit

    Required arguments:
      --input INPUT, -i INPUT
                            Input alignment file of SAM/BAM/CRAM format MUST be
                            paired end and sorted
      --output OUTPUT, -o OUTPUT
                            Output structure profile, referred to as .sprofile
                            files, is a tab delimited table.

    Additional arguments:
      --insert_max [INSERT_MAX], -im [INSERT_MAX]
                            Maximum size of input DNA fragments, over which
                            deletions are annotated, default=2000
      --fragment_sizes FRAGMENT_SIZES, -fs FRAGMENT_SIZES
                            Optional output file for the fragment size of
                            correctly mapped reads



<br>
<br>
<br>
<br>
<br>



## Performance
The mProfile tools are very lightweight, they require only 10-20MB of system memory and are generally single-threaded.

It is recommened to locally store all files to be processed, i.e. on fast storage directly attached to the processing computer, not on network attached storage. Files are read line by line, therefore slow storage can result in significant slow downs.

#### Threading
When a control file is specified for normalisation, mprofile tools run two threaded, simultaneously processing both samples.<br>
There is no other way to multi-thread these tools.

Since the tools are lightweight and single runs are relatively fast (see below), it's recommended that to improve speed to simultaneously process multiple runs via command line.

#### Expected speed
callMUT processing time is mostly altered by the number of nucleotides analysed, however the number of reads originally analysed can also have an effect.

Including a control sample for normalisation will cause increased processing time, though not double.

Processing time for a 4500 nucleotide, 18 million read file was 3min 45s.<br>
The above example done again but with a control file took 6min 4s

An important note is that the processing time for indels can cause amplicons with high levels of indels compared to the reference genome (e.g. a common alteration in the cell line) to result in substantial delays. This normally requires the indels to be in a majority of reads and so will not impact most cases. 
<br><br><br>
TransloCapture is not as effected by having a control file as callMUT.<br>
Processing time for 15.9M reads with 21 primer pairs with a control file was 2min 56s, ~5 million reads/minute. <br>
Without a control file, this process took 2min 53s.

Both readcount and number of primer pairs almost linearly alter processing time.<br>
The above example done with 42 primer pairs took 5min.
