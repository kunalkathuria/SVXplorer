# SVXplorer
A structural variant caller that uses discordant read-pairs (PE), split-reads (SR) and read-depth (RD) information.

### SUMMARY

SVXplorer accepts a BAM file of target as input and outputs a BEDPE file (and an equivalent VCF file) containing deletions (DEL), tandem duplications (TD), inversions (INV), non-tandem-duplications, translocations (see below), novel sequence insertions (DN_INS) and undetermined types tagged as "BND." The variant tags listed in parentheses pertain to the BEDPE file, whereas the VCF file follows VCF 4.3 specifications. The insertions contain 3 breakpoints (1 = source location 1, 2 = source location 2, 3 = paste location) and are tagged as INS (copy-paste insertion/duplication), INS_C_P (cut-paste insertion or translocation that is fully identified), INS_C (cut-paste insertion with ambiguity in 2 of its breakpoints, i.e. source vs paste location) or INS_I (inverted copy-paste insertion) in the BEDPE file. The BND events contains all variants that are not fully identified, like INS_C, partially supported inversions, non-deletion "FR" clusters, deletions or duplications not supported by local depth of coverage, interchromosomal insertions only containing 2 breakpoints etc. and are indicated as such in the comment/INFO field in VCF/BEDPE respectively.

SVXplorer addresses many of the common limitations in SV detection by its comprehensive 3-tier approach of sequentially using discordant paired-end (PE) alignment, split-read (SR) alignment and read-depth information to capture as many SVs as possible and simultaneously or progressively weeding out poor candidates. Significant attention is given to categorizing alignments, grouping alignments correctly into respective clusters, eliminating cluster “conflict,” consolidating clusters meticulously into variants, integrating PE and SR calls precisely, dynamically calculating PE and SR SV-support thresholds, retaining all clusters and choosing variants based on final support, and corroborating SVs using streamlined local read-depth information.

### METHODOLOGY

SVXplorer first forms discordant clusters from paired-end reads via formation of maximal cliques in a weight-thresholded bidirectional graph and consolidates them further into PE-supported variants. It then integrates split reads and read-depth information to call putative variants, enhancing/filtering out existing variants or identifying new ones along the way. 

### READ-GROUPS

SVXplorer should be run on a BAM file with a single read-group. For datasets with multiple read-groups, we recommend running SVXplorer on every read-group separately. We will add support for consolidation of these results soon. Please follow issue #81 for further details.

### REQUIREMENTS

SVXplorer should run on any Unix-based OS with bash, python (version 2.6 or 2.7) and libraries as specified in the requirements file. "bedtools" and "samtools" executables should be on the user PATH. In addition, if a split read file is not available in the typical splitters format (2 entries per query name with 2 distinct, split queries) a script is provided to extract this from the alignment file using LUMPY's extractBwaMem_reads script (https://raw.githubusercontent.com/arq5x/lumpy-sv/master/scripts/extractSplitReads_BwaMem).

### INSTALLATION

Download latest SVXplorer release from GitHub and unzip the directory. Alternatively, clone the repository. Then ensure that "samtools" and "bedtools" are on the PATH by running 

```
which samtools
which bedtools
``` 

Then install all the python libraries and install SVXplorer using

```
cd SVXplorer*
pip2 install -r requirements.txt
make
```

### USAGE

SVXplorer can now be run from the "bin" sub-directory of the installation (which can be added to user path to run it from any directory).

```
./bin/SVXplorer -h
```

will show all options available to the user. Then run

```
make test
```

to run a simple test-case which will ensure that SVXplorer is running as expected. SVXplorer will run using the test alignment files in the SVXplorer/testCases folder and check if the resulting vcf file found in the newly created "test" folder is identical with the one contained in SVXplorer/testFiles. A message will be printed notifying whether the test was successful.

The input BAM file should be generated and then indexed with BWA (or potentially any aligner). For example, using standard formatting options,

```bash
$BWA mem -R '@RG\tID:foo\tSM:bar' -Y -t $threads $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> input.bam
```

In addition, a split-alignment file and a file containing all fragments that align discordantly (using "samtools view -F 3842") should be provided as input. If unsure, scripts are provided in the scripts/ folder to create a discordant alignment file as well as a split-read alignment file (the latter needs a LUMPY script found in PATH_TO_LUMPY/scripts). To create split reads, use the script as follows: ./createSplitReads.sh BAMFILE OUTPUT_FILE_PATH PATH_TO_LUMPY N_THREADS (please make sure not to have soft-clipped "secondary" alignments as generated by -M option in BWA present in BAMFILE). To create discordants, run: ./createDiscordants.sh BAMFILE OUTPUT_DISC_FILE_PATH N_THREADS.

All SVXplorer command line options are accessed via ./SVXplorer -h. A file to ignore alignments in certain chromosome/genomic units for Human Genome Reference Build b37 (-c) and a file to exclude certain regions of alignment for b37 (-i) are included in the data folder. The reference file should be indexed using "samtools faidx." 

A typical (and recommended) call on sequenced data might be, with appropriate path replacement:

path_to_SVXplorer/bin/SVXplorer discordant.bam splitters.bam sample.bam reference.fa -i exclude.bed -c ignore_CHR.txt -m non_repeat_regions.bed -w pathToWorkingDirectory

Option -m expects a file listing regions not containing frequently repeated sequences, for use in assessing coverage, and SVC provides one for b37 in the data folder in zipped form, which will need to be unzipped before use.

### RESULTS

The final BEDPE and VCF files are written in the results folder. Intermediate BEDPE files using only PE/PE and SR mappings are stored in the "workspace" directory.

### TUNABLE PARAMETERS

For the advanced user, there are certain parameters in SVXplorer that can be manipulated and tailored to the nature of the input data. These are not accessible from the command line. Unless noted otherwise, they are found in the sequential global variable list in bin/SVXplorer.

1. MAP_THRESH: minimum mapping quality for a PE alignment to enter the SVXplorer pipeline. By default, this is set to 1, which implies that reads where the given alignment is more likely than the second best alignment are included. Since mapping quality is the phred-scaled posterior probability that the mapping position of this read is incorrect, the user might want to modify this to improve precision, based on their needs.

2. PE_THRESH_MIN, PE_THRESH_MAX: minimum and maximum cutoffs for the dynamic support threshold calculation via linear model (please see manuscript for details) for PE-only variant calls. We require at least 3 paired-end reads to support the same variant, before it is considered. As the coverage increases, the minimum reuired coverage is dynamically selected between [PE_THRESH_MIN, PE_THRESH_MAX]. Increasing these should increase precision, though it may lead to decrease in sensitivity.

3. SR_THRESH_MIN, SR_THRESH_MAX: same as above for variant calls that are only supported by split-reads. By default these are set to 3 and 6 respectively.

4. MQ_SR: minimum mapping quality for split alignments to be considered. By default this is set to 10.

5. SLOP_SR: maximum distance between reference locations of split alignments to become part of the same variant. Split read alignments can be heavily influenced by the sorrounding bases. If two split alignments are within SLOP_SR bases, then they are considered to be supporting the same variant. If the user runs a local indel realignment, then this can be reduced without a loss in sensitivity. 

6. DEL_CN_SUPP_THRESH (DUP_CN_SUPP_THRESH): variant-region coverage threshold to call a deletion (duplication) as a fraction of median chromosome coverage. Local coverage lower (higher) than this value will cause a deletion (duplication) to be called. By default this is set to 0.8(1.15). The user should maniplulate these based on the coverage distribution of their dataset.

7. MIN_PILEUP_THRESH in bin/covPUFilter.py: specifies the minimum number of mappable bases that should be used to calculate local coverage in a variant region for the thresholds in 6) to be effective. ("MIN_PILEUP_THRESH_NH" in the same file should be set to the same value.) This value is set to 80 by default.

8. pe_low, pe_high, sr_low, sr_high, covg_low, covg_high in bin/uniqueSupportFilter.py: (covg_low,pe_low) and (covg_high, pe_high) define the (x,y) points to be joined by a line and interpolated between to set the dynamic support threshold given the coverage of the data in question (please see manuscript for details on the model). Same for sr_high, sr_low. Changing these parameters alters the slope/position of the linear model.

Users familiar with their data and SVXplorer scripts can alter these parameters to suit the nature of their input data. For example, if one is processing very high-quality alignment data or is processing simulated data for specific analysis, one could lower the above mapping quality thresholds (in addition to choosing a low threshold via the "-q" option in command-line) to any value greater than 0. One can also use lower read-support values corresponding to "covg_low" and "covg_high" to alter the dynamic support model calculated by SVXplorer. The fixed local coverage thresholds to call deletions and duplications can also be made more strict or relaxed depending upon the needs of the user and quality of input data. MIN_PILEUP_THRESH can be manipulated based on the scale precision of the user's mappable region file (if using one different from SVXplorer's). It can be increased to be conservative or decreased to be more liberal in filtering variants. 

We plan to expose more of these variables as we document their effects on the final results.
