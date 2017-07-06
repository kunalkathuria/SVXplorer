# VARSECer
A comprehensive-approach structural variant caller making use of all secondary discordant PE alignments.

### SUMMARY

VARSECer accepts a BAM file of target as input and outputs 6 bedpe files (and 1 equivalent vcf file) containing deletions, non-tandem duplications (insertions.bedpe), inversions, tandem duplications, novel sequence insertions (insertions_dn.bedpe) and BNDs repectively in standard bedpe format. The insertions are tagged as INS (copy-paste duplication), INS_C_P (cut-paste or translocation), INS_I (inverted copy-paste), INS_C_P, INS_C_I_P. The BND file ("unknowns.bedpe") contains all unidentified variants, like partially supported inversions, deletions etc. It also contains translocations (INS_C) whose source and insert locations have not been ascertained (Breakpt 1 and Breakpt 3 are the source and insert locations, but it could not be determined which they were respectively). In the last column, the partially identified variant type is also provided.

### METHODOLOGY

One unique feature of VARSECer is that it uses all secondary PE alignments above a certain threshold to form variants. VARSECer forms discordant clusters with specific inclusion/filtering criteria and classifies them into all basic variant types after grouping them with all possible cluster matches. It then uses an approach based on unique-fragment-support to call its final variants-- the basic idea being that if a variant has sufficient unique support, it is likely to be true. It uses an iterative approach of using paired-end mappings, split reads and read-depth-based INDEL calls respectively to enhance existing variants and identiy new ones. For the latter, it integrates third-party read-depth-based INDEL calls with existing variants. For split reads, a name-sorted split read file (no secondary alignments) should be provided as mentioned below.

### REQUIREMENTS

Unix-based OS with bash, python with pysam and other common libraries, samtools, git. It is best to have basic executables like "python", "samtools" etc. on user path. 

### INSTALLATION

Download VARSECer from GitHub following GitHub instructions. Follow usage instructions below. 

### USAGE

Example call (from sv_caller/src):

./run_sv_caller.sh [options]

A good check to see if the tool is working properly is to give it as input the test BAM files in the sv_caller/data/bams/test folder and check if the resulting bedpe files found in the sv_caller/results/text folder are identical with those contained in the sv_caller/results/test folder. So, first run the code with all default parameter values thus from the sv_caller/src folder:

./run_sv_caller.sh -a ../data/bams/test/test_1k.bam -b ../data/bams/test/test_1k.ns.bam -r ../data/bams/test/test_1k_splitters.ns.bam -z ../data/bams/test/discordants.bam -A 1 -E 0

Then, ascertain that the 4 bedpe files in sv_caller/results/test match the respective ones just created in sv_caller/results/text (there should be 2 non-empty ones: deletions.bedpe and tandemDuplications.bedpe). If this is not the case, please recheck all the paths and checkout the master branch from GitHub again if necessary.

The input BAM file should be generated (and then indexed) with BWA using the following options:

```bash
$BWA mem -R '@RG\tID:foo\tSM:bar' -a -Y -t 1 $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> input.bam
```

In addition, a name-sorted file containing *all* fragments that align discordantly (if unsure, simply uncomment the lines in ./run_sv_caller.sh and provide a dummy input to -z dummy_disc_file) and a name-sorted split-reads BAM file (uncomment the lines in run_RD_SR.sh if unavailable and have the indicated LUMPY SR script handy) should also be available for required input. If using CNV (read-depth) signal as input, set the RDS flag (-m) to 1 and provide a space or tab-delimited file in the format (chr start stop copy-number) as input (-x). The script run_cnv.filters.sh is provided for this purpose assuming CNVNATOR installed. Only the first 4 columns are required.

All VARSECer command line options are accessed via (./run_sv_caller.sh or ./run_sv_caller.sh -h)

### RESULTS

Six bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (chr_source, start, stop, chr_insert, start, end) non-tandem duplication (NTD) format. The kind of NTD (cut,copy,inverted cut etc.) is identified in the last column as indicated above. 1 vcf file is also created containing all variants together.

These files are stored in the sv_caller/results/text directory. Intermediate bedpe files using only PE (PE and SR) mappings are stored in the above directory as well in a folder titled 'pe_results' (sr_results). These results may give better results against "truth sets" if the truth set's breakpoint precision is questionable.

### NOTES

1. To enhance sensitivity of variant calls, one can set MQ_THRESH to 0 in the main script. By default it is set to 1, which provides a faster run-time than most commonly used tools.
