# VARSECer
A comprehensive-approach structural variant caller making use of PE, SR and read-depth information.

### SUMMARY

VARSECer accepts a BAM file of target as input and outputs 6 bedpe files (and 1 equivalent vcf file) containing deletions, non-tandem duplications (insertions.bedpe), inversions, tandem duplications, novel sequence insertions (insertions_dn.bedpe) and BNDs repectively in standard bedpe format. The insertions are tagged as INS (copy-paste duplication), INS_C_P (cut-paste or translocation), INS_I (inverted copy-paste), INS_C_P, INS_C_I_P. The BND file ("unknowns.bedpe") contains all unidentified variants, like partially supported inversions, small non-deletion "FR" clusters etc. It also contains translocations (INS_C) whose source and insert locations have not been ascertained. In the last column, the partially identified variant type is also provided.

### METHODOLOGY

VARSECer first forms discordant clusters from paired-end reads via formation of maximal cliques in a weight-thresholded bidirectional graph and matches them up with each other to form PE variants. It then integrates split reads and read-depth information from samtools(c) to call putative variants. Thus, it iteratively uses paired-end mappings, split reads and read-depth-based pile-up filters to enhance existing variants and identiy new ones. It also filters out variants based on unique-fragment-support criteria: if a variant has sufficient support from fragments that uniquely support that variant, it is likely to be true. 

### REQUIREMENTS

Unix-based OS with bash, python with pysam and other common libraries, samtools, git. It is best to have basic executables like "python", "samtools" etc. on user path. 

### INSTALLATION

Download VARSECer from GitHub following GitHub instructions. Follow usage instructions below. 

### USAGE

Simply type "make" from the installation path. Now you are ready to run VARSECer.

./varsecer [options]

A good check to see if the tool is working properly is to provide as input the test alignment files in the VARSECER/data/bams/test folder and check if the resulting bedpe files found in the working directory are identical with those contained in the VARSECER/results/test folder. So, first run the code thus (in this example the current directory is the VARSECER parent directory):

./bin/varsecer -a ../data/bams/test_1k.bam -r ../data/bams/test_1k_splitters.ns.bam -z ../data/bams/discordants.bam -A 1 -m 0 -E 0 -W ../results/text

Then, ascertain that the bedpe files in VARSECER/results/test match the respective ones just created in the working directory (there should be 3 non-empty ones: deletions.bedpe, tandemDuplications.bedpe and unknowns.bedpe). If this is not the case, please recheck all the paths and checkout the master branch from GitHub again if necessary.

The input BAM file should be generated and then indexed with BWA (or potentially any aligner). For example, using standard formatting options,

```bash
$BWA mem -R '@RG\tID:foo\tSM:bar' -Y -t $threads $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> input.bam
```

In addition, a name-sorted file containing fragments that align discordantly (if unsure, simply uncomment the lines in ./varsecer and provide a dummy input to -z dummy_disc_file) and a name-sorted split-reads BAM file (uncomment the lines in run_RD_SR.sh if unavailable and have the indicated LUMPY SR script handy) should also be available for required input. If using CNV (read-depth) signal as input, set the RDS flag (-m) to 1 and provide a space or tab-delimited file in the format (chr start stop copy-number) as input (-x). The script run_cnv.filters.sh is provided for this purpose assuming CNVNATOR is installed. Only the first 4 columns are required. If one does not wish to use third party CNV calls to call DE NOVO variants, but only to support existing ones, simply use -F 0 when running VARSECer. Using 3rd party read-depth calls is not highly recommended due to the unreliability of reported breakpoints.

All VARSECer command line options are accessed via (./varsecer or ./varsecer -h). A file to ignore alignments in certain chromosome/genomic units for B37 (-C), a file to exclude certain regions of alignment for B37 (-D) and a file listing regions not containing high repeats for use in assessing coverage (-I) are included in the VARSECer/data/bed and VARSECER/data/txt.

A typical (and recommended) call on sequenced data might be, with appropriate path replacement:

./varsecer -A 0 -B 32 -z disc.bam -a sample.bam -r splitters.ns.bam -C exclude.bed -D ignore_CHR.txt -I good_regions.bed -W ../results/text

### RESULTS

Six bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (chr_source, start, stop, chr_insert, start, end) non-tandem duplication (NTD) format. The kind of NTD (cut,copy,inverted cut etc.) is identified in the last column as indicated above. 1 vcf file (All_SVs.vcf) is also created containing all variants together.

These files are written in the working directory. Intermediate bedpe files using only PE (PE and SR) mappings are stored in the above directory as well in a folder titled 'pe_results' (sr_results). These results may give better results against "truth sets" if the truth set's breakpoint precision is questionable.

### NOTES

1. To enhance sensitivity  and precision of variant calls for good data sets, one can set MQ_THRESH to 1. By default it is set to 10, but setting it to 1 would work well for high quality data (low spread in insert length, good coverage, low bias etc.). 
2. One can also set MQ_THRESH to 0 and use the option of secondary alignments by assigning a value less than 1 to MATCHRATIO (-d;see command line help) for higher sensitivity, but which may increase run-time considerably.
3. One can set SPLIT_INS to True for good-quality data which allows for better performance with this particular filter (reevaluates non-tandem duplications/translocations and break them into deletions and tandem duplications if indicated by local read depth). 
