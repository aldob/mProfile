# mProfile
## Overview
A package for processing targeted sequencing (amplicon) data for high-resolution mutation profiling. 

The callMUT tool converts mutation data from the Samtools mpileup format, which contains raw mutation calls for every read, into a lightweight, analysis ready mprofile: a table cotaining the percentage rate of every mutation type at every nucleotide. 

The TransloCapture tool maps all possible crossover events between different targets in a targeted sequencing experiment. 



<br>
<br>
<br>
<br>
<br>



## callMUT
The default output of Samtools Mpileup is very large and cumbersome to parse, a problem that increases with alignment file size. 

Mprofiles are comparitively very lightweight and do not expand with increasing alignment file size. They provide per-nucloetide mutation rates (% of reads) in a table that makes plotting and further analysis simple and fast! 

Also able to calculate the differential between a treated(-i) and a control(-c) sample (either in mpileup or mprofile format). NOTE: To calculate the differential the mpileups must be created using the argument -aa to ensure they all report the exact same coordinates.

#### Arguments
    mprofile callMUT -i <input> -o <output>
    
    Required arguments:
      --input -i INPUT
                        Input mpileup/mprofile to process   
      --output -o OUTPUT
                        Output mprofile
    Additional arguments:
      --indelcut -ic [INDELCUT]
                        Minimum rate (% of reads) for an indel to be reported,
                        default=1
                        If 'NA', no cutoff will be applied
      --control -c CONTROL
                        mpileup/mprofile to normalise to (e.g. untreated)
      --preproc -pp PREPROC
                        Specifies input files are mprofiles, not mpileups
                        (requires --control to be set)
      --help -h HELP
                        show this help message and exit


    Example: 
      mprofile callMUT -i treated.mpileup -c untreated.mpileup -o treated.mprofile -ic 0.001

<br>

#### Example .mprofile table
With -ic NA set, all indels found in the input file will be reported in the Common.Indels column.

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
In the multiplex PCR used for targeted sequencing, crossover events (translocations) between the targets will also be sequenced. TransloCapture takes unprocessed fastq files and identifies reads that are crossover events. The tool works by identifying the primers used for amplification of the targets at either end of the reads.

Requires a csv file(-p) containing three columns "target name", "forward primer sequence", "reverse primer sequence" for each target. The output is a .csv file containing a matrix of all possible crossovers and their frequencies (% of reads) and can also output the crossover reads to a new fastq file.
NOTE: The feature of writing the translocated reads to a new fastq is not currently available.

Similar to callMUT, it can also calculate the differential of the input to a control and can accept pre-generated output .csv files instead of fastq to generate these differentials seperately (-pp). 

Accepts either single-read(SR) or paired-end(PE) data, however PE is highly reccommended for this analysis! SR will significantly reduce the accruacy of the tool.

#### Arguments
    mprofile TransloCapture -1 <input_read1> -2 <input_read2> -o <output>
    
    Required arguments:
      --input -i INPUT
                        Fastq from SR data
      --read1 -1 READ1
                        Fastq read 1 from PE data
      --read2 -2 READ2
                        Fastq read 2 from PE data
      --output -o OUTPUT
                        .csv file to write result matrix
      --primers -p PRIMERS
                        A .csv file containing the name of the target and the foward and reverse
                        primer (reverse complement) sequences for each site to be analysed 
    Additional arguments:
      --control -c CONTROL
                        The fastq you want to normalise to (e.g. untreated)
      --control1 -c1 CONTROL1
                        Read 1 of the fastq you want to normalise to (e.g.
                        untreated)
      --control2 -c2 CONTROL2
                        Read 1 of the fastq you want to normalise to (e.g.
                        untreated)
      --preproc -pp PREPROC
                        Specifies input files are pre-made .csv translocation output files
                        (requires --control to be set)
      --fastq -fq FASTQ
                        Fastq file to write translocated reads to. Default will not write the reads.
      --help -h HELP
                        show this help message and exit


    Example: 
      mprofile TransloCapture -1 treated_R1.fastq -2 treated_R2.fastq -o treated_translomap.csv -p target_primers.csv
      
      
       -h, --help            show this help message and exit

<br>

#### Example Primer table
First, here is an example primer table that we will use to create the example translocation table.<br>
Columns are target name, forward primer sequence, reverse primer sequence.

.          |.                         |.                           |
-----------|--------------------------|----------------------------|
CTRL_Target|TACGCTACGACTAGCAGCTATCGACT|ATCGCATCTAGACTGATCACGATCTACG|
Target1|TACGACACTACGACTATCGA|ATCGCTACGACGATAGCTTCGT|
Target2|GCTACGCTACAGCGACTACTA|GAGCTACGACATCACCGCTGCATCA|
Target3|GACTACGACTACGACTACG|CGCGACTACGACGCTACGCGC|
Target4|CGCGACTACGCGATCACGACTACTAGC|CGATCACGACTACGACGCATCG|

<br>

#### Example Crossover table
The x-axis are the sites recognised in Read1(or at the start if SR) and the y-axis are the sites recognised in Read2(or at the end if SR).<br>
Cells where x=y are the un-translocated targets and are represented by NA as the counts are very high. 

Site|CTRL_Target|Target1|Target2|Target3|Target4
----|-----|----|----|----|----|
CTRL_Target|NA|0|0|0|0|
Target1|0|NA|0.01673|0.0002336|0.23451|
Target2|0|0.010564|NA|0.01568|0.006168|
Target3|0|0.000005169|0.000168168|NA|0.000079841|
Target4|0|0.0016546|0.078965|0.0165169|NA|
