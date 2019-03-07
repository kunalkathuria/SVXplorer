# SVXplorer
A comprehensive-approach structural variant caller making use of PE, SR and read-depth information.

### SUMMARY

SVXplorer accepts a BAM file of target as input and outputs a BEDPE file (and an equivalent VCF file) containing deletions (DEL), tandem duplications (TD), inversions (INV), non-tandem-duplications, translocations (see below), novel sequence insertions (DN_INS) and undetermined types tagged as "BND." The variant tags listed in parentheses pertain to the BEDPE file, whereas the VCF file follows VCF 4.3 specifications. The insertions contain 3 breakpoints (1 = source location 1, 2 = source location 2, 3 = paste location) and are tagged as INS (copy-paste insertion/duplication), INS_C_P (cut-paste insertion or translocation that is fully identified), INS_C (cut-paste insertion with ambiguity in 2 of its breakpoints, i.e. source vs paste location) or INS_I (inverted copy-paste insertion) or in the BEDPE file. The BND events contains all variants that are not fully identified, like INS_C, partially supported inversions, non-deletion "FR" clusters, deletions or duplications not supported by local depth of coverage, interchromosomal insertions only containing 2 breakpoints etc. and are indicated as such in the comment/INFO field respectively.

SVXplorer addresses many of the standard limitations in SV detection by its comprehensive 3-tier approach of sequentially using discordant paired-end (PE) alignment, split-read (SR) alignment and read-depth information to capture as many SVs as possible and simultaneously or progressively weeding out poor candidates. Significant attention is given to categorizing alignments, grouping alignments correctly into respective clusters, eliminating cluster “conflict,” consolidating clusters meticulously into variants, integrating PE and SR calls precisely, dynamically calculating PE and SR SV-support thresholds, retaining all clusters and choosing variants based on final support, corroborating SVs using streamlined local read-depth information etc.

### METHODOLOGY

SVXplorer first forms discordant clusters from paired-end reads via formation of maximal cliques in a weight-thresholded bidirectional graph and consolidates them further into PE-supported variants. It then integrates split reads and read-depth information to call putative variants, enhancing/filtering out existing variants or identifying new ones along the way. 

### REQUIREMENTS

Unix-based OS with bash, python2 with basic and some other libraries (including pysam, pybedtools, pandas, networkx, bitarray, interlap, scikit-learn), "bedtools" command on user path. 

### INSTALLATION

Download SVXplorer from GitHub following GitHub instructions. Follow usage instructions below. 

### USAGE

Type "make" from the installation path. Now you are ready to run SVXplorer.

./bin/SVXplorer [options]

A good check to see if the tool is working properly is to provide as input the test alignment files in the SVXplorer/testCases folder and check if the resulting vcf file found in the newly created "test" folder is identical with the one contained in SVXplorer/testFiles. Simply type:

make test

and a message will be printed notifying whether the test was successful.

The input BAM file should be generated and then indexed with BWA (or potentially any aligner). For example, using standard formatting options,

```bash
$BWA mem -R '@RG\tID:foo\tSM:bar' -Y -t $threads $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> input.bam
```

In addition, a split-alignment file and a file containing all fragments that align discordantly should be provided as input. If unsure, scripts are provided in the scripts/ folder to create an discordant alignment file as well as a split-read alignment file (the latter uses a LUMPY script). To create split reads, use the script as follows: ./createSplitReads.sh BAMFILE OUTPUT_DIR OUTPUT_FILE_NAME PATH_TO_LUMPY N_THREADS. One may alternatively use SAMBLASTER: samtools view -h samp.bam | samblaster -a -e -s samp.split.sam -o /dev/null. To create discordants, run: ./createDiscordants.sh BAMFILE OUTPUT_DISC_FILE_PATH N_THREADS.

All SVXplorer command line options are accessed via ./SVXplorer -h. A file to ignore alignments in certain chromosome/genomic units for B37 (-c) and a file to exclude certain regions of alignment for B37 (-i) are included in the data folder. 

A typical (and recommended) call on sequenced data might be, with appropriate path replacement:

path_to_SVXplorer/bin/SVXplorer discordant.bam splitters.bam sample.bam reference.fa -i exclude.bed -c ignore_CHR.txt -m non_repeat_regions.bed -w pathToWorkingDirectory

Option -m expects a file listing regions not containing frequently repeated sequences, for use in assessing coverage, and SVC provides one in the data folder in zipped form.

### RESULTS

The final BEDPE and VCF files are written in the results folder. Intermediate BEDPE files using only PE/PE and SR mappings are stored in the "workspace" directory.

