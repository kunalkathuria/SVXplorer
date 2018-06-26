# SVXplorer
A comprehensive-approach structural variant caller making use of PE, SR and read-depth information.

### SUMMARY

SVXplorer accepts a BAM file of target as input and outputs 6 bedpe files (and 1 equivalent vcf file) containing deletions, non-tandem duplications (insertions.bedpe), inversions, tandem duplications, novel sequence insertions (insertions_dn.bedpe) and BNDs in standard bedpe format. The insertions are tagged as INS (copy-paste duplication), INS_C_P (cut-paste or translocation), INS_I (inverted copy-paste), INS_C_P, INS_C_I_P. The BND file ("unknowns.bedpe") contains all unidentified variants, like partially supported inversions, small non-deletion "FR" clusters, translocations (INS_C) whose source and insert locations have not been ascertained etc.

SVXplorer addresses many of the standard limitations in SV detection by its comprehensive 3-tier approach of sequentially using discordant paired-end (PE) alignment, split-read (SR) alignment and read-depth information to capture as many SVs as possible and simultaneously or progressively weeding out poor candidates. Significant attention is given to categorizing alignments, grouping alignments correctly into respective clusters, eliminating cluster “conflict,” consolidating clusters meticulously into variants, integrating PE and SR calls precisely, dynamically calculating PE and SR SV-support thresholds, retaining all clusters and choosing variants based on final support, corroborating SVs using streamlined local read-depth information etc.

### METHODOLOGY

SVXplorer first forms discordant clusters from paired-end reads via formation of maximal cliques in a weight-thresholded bidirectional graph and consolidates them further into PE-supported variants. It then integrates split reads and read-depth information to call putative variants, enhancing/filtering out existing variants or identifying new ones along the way. 

### REQUIREMENTS

Unix-based OS with bash, python with pysam and other common libraries, samtools, git. It is best to have basic executables like "python", "samtools" , "bedtools" etc. on user path. 

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

In addition, a split-alignment file and a file containing all fragments that align discordantly should be provided as input (e.g. samtools view -F 3586 -b -o DISC_OUTPATH BAM). If unsure, scripts are provided in the scripts/ folder to create discordants as well as split reads (latter uses a LUMPY script).

All SVXplorer command line options are accessed via ./SVXplorer -h. A file to ignore alignments in certain chromosome/genomic units for B37 (-c) and a file to exclude certain regions of alignment for B37 (-i) are included. 

A typical (and recommended) call on sequenced data might be, with appropriate path replacement:

path_to_SVXplorer/bin/SVXplorer discordant.bam splitters.bam sample.bam reference.fa -i exclude.bed -c ignore_CHR.txt -m non_repeat_regions.bed -f -w pathToworkDir -s 100

Option -m expects a file listing regions not containing frequently repeated sequences, for use in assessing coverage, and SVC provides one.

### RESULTS

Six bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (chr_source, start, stop, chr_insert, start, end) non-tandem duplication (NTD) format. The kind of NTD (cut,copy,inverted cut etc.) is identified in the last column as indicated above. 1 vcf file (All_SVs.vcf) is also created containing all variants together.

These files are written in the working directory. Intermediate bedpe files using only PE (PE and SR) mappings are stored in the above directory as well in a folder titled 'pe_results' (sr_results). These results may give better results against "truth sets" if the truth set's breakpoint precision is questionable.

### NOTES

1. One can set SPLIT_INS to True for *good-quality* diploid data which allows for better performance with this particular filter (reevaluates non-tandem duplications/translocations and break them into deletions and tandem duplications if indicated by local read depth).
2. One can set LIB_INV to FALSE if wish to call only those inversions that are supported by 2 opposing PE clusters, and not 1 PE and 1 SR cluster only (see manuscript for more details). 
