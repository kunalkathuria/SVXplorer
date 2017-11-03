# VARSECer
A comprehensive-approach structural variant caller making use of all secondary discordant PE alignments.

### SUMMARY

VARSECer accepts a BAM file of target as input and outputs 6 bedpe files (and 1 equivalent vcf file) containing deletions, non-tandem duplications (insertions.bedpe), inversions, tandem duplications, novel sequence insertions (insertions_dn.bedpe) and BNDs repectively in standard bedpe format. The insertions are tagged as INS (copy-paste duplication), INS_C_P (cut-paste or translocation), INS_I (inverted copy-paste), INS_C_P, INS_C_I_P. The BND file ("unknowns.bedpe") contains all unidentified variants, like partially supported inversions, deletions etc. It also contains translocations (INS_C) whose source and insert locations have not been ascertained (Breakpt 1 and Breakpt 3 are the source and insert locations, but it could not be determined which they were respectively). In the last column, the partially identified variant type is also provided.

### METHODOLOGY

VARSECer forms discordant clusters via formation of maximal cliques in a weight-thresholded bidirectional graph and classifies them into all basic variant types after grouping them with all possible cluster matches. It then uses an approach based on unique-fragment-support to call its final variants-- the basic idea being that if a variant has sufficient unique support, it is likely to be true. It uses an iterative approach of using paired-end mappings, split reads and read-depth-based pile-up filters respectively to enhance existing variants and identiy new ones. For split reads, a name-sorted split read file (no secondary alignments) should be provided as mentioned below.

### REQUIREMENTS

Unix-based OS with bash, python with pysam and other common libraries, samtools, git. It is best to have basic executables like "python", "samtools" etc. on user path. 

### INSTALLATION

Download VARSECer from GitHub following GitHub instructions. Follow usage instructions below. 

### USAGE

Example call (from VARSECER/src):

./varsecer [options]

A good check to see if the tool is working properly is to give it as input the test BAM files in the VARSECER/data/bams/test folder and check if the resulting bedpe files found in the VARSECER/results/text folder are identical with those contained in the VARSECER/results/test folder. So, first run the code with all default parameter values thus from the VARSECER/src folder:

./varsecer -a ../data/bams/test/test_1k.bam -b ../data/bams/test/test_1k.ns.bam -r ../data/bams/test/test_1k_splitters.ns.bam -z ../data/bams/test/discordants.bam -A 1 -m 0 -E 0

Then, ascertain that the 4 bedpe files in VARSECER/results/test match the respective ones just created in VARSECER/results/text (there should be 2 non-empty ones: deletions.bedpe and tandemDuplications.bedpe). If this is not the case, please recheck all the paths and checkout the master branch from GitHub again if necessary.

The input BAM file should be generated (and then indexed) with BWA (or potentially any aligner listing all alignments of each read together)  using the following options:

```bash
$BWA mem -R '@RG\tID:foo\tSM:bar' -a -Y -t 1 $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> input.bam
```

In addition, a name-sorted file containing *all* fragments that align discordantly (if unsure, simply uncomment the lines in ./varsecer and provide a dummy input to -z dummy_disc_file) and a name-sorted split-reads BAM file (uncomment the lines in run_RD_SR.sh if unavailable and have the indicated LUMPY SR script handy) should also be available for required input. If using CNV (read-depth) signal as input, set the RDS flag (-m) to 1 and provide a space or tab-delimited file in the format (chr start stop copy-number) as input (-x). The script run_cnv.filters.sh is provided for this purpose assuming CNVNATOR is installed. Only the first 4 columns are required. If one does not wish to use third party CNV calls to call DE NOVO variants, but only to support existing ones, simply use -F 0 when running VARSECer.

All VARSECer command line options are accessed via (./varsecer or ./varsecer -h)

A typical call on real sequenced data might be, with appropriate path replacement:

./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 1

### RESULTS

Six bedpe files are created in the standard 6-column format (chr1, start, stop, chr2, start, stop) except for insertions.bedpe which is in the standard (chr_source, start, stop, chr_insert, start, end) non-tandem duplication (NTD) format. The kind of NTD (cut,copy,inverted cut etc.) is identified in the last column as indicated above. 1 vcf file (All_SVs.vcf) is also created containing all variants together.

These files are stored in the VARSECer/results/text directory. Intermediate bedpe files using only PE (PE and SR) mappings are stored in the above directory as well in a folder titled 'pe_results' (sr_results). These results may give better results against "truth sets" if the truth set's breakpoint precision is questionable.

### NOTES

1. To enhance sensitivity  and precision of variant calls, one can set MQ_THRESH to 0 in the main script. By default it is set to 1, which provides a faster run-time than most commonly used tools with  coverage > 50x. But for lower coverages especially, it will be faster and more beneficial.
