# SVXplorer
A structural variant caller that uses discordant read-pairs (PE), split-reads (SR) and read-depth (RD) information.

### SUMMARY

SVXplorer accepts a BAM file of target as input and outputs a BEDPE file (and an equivalent VCF file) containing deletions (DEL), tandem duplications (TD), inversions (INV), non-tandem-duplications, translocations (see below), novel sequence insertions (DN_INS) and undetermined types tagged as "BND." The variant tags listed in parentheses pertain to the BEDPE file, whereas the VCF file follows VCF 4.3 specifications. The insertions contain 3 breakpoints (1 = source location 1, 2 = source location 2, 3 = paste location) and are tagged as INS (copy-paste insertion/duplication), INS_C_P (cut-paste insertion or translocation that is fully identified), INS_C (cut-paste insertion with ambiguity in 2 of its breakpoints, i.e. source vs paste location) or INS_I (inverted copy-paste insertion) or in the BEDPE file. The BND events contains all variants that are not fully identified, like INS_C, partially supported inversions, non-deletion "FR" clusters, deletions or duplications not supported by local depth of coverage, interchromosomal insertions only containing 2 breakpoints etc. and are indicated as such in the comment/INFO field respectively.

SVXplorer addresses many of the common limitations in SV detection by its comprehensive 3-tier approach of sequentially using discordant paired-end (PE) alignment, split-read (SR) alignment and read-depth information to capture as many SVs as possible and simultaneously or progressively weeding out poor candidates. Significant attention is given to categorizing alignments, grouping alignments correctly into respective clusters, eliminating cluster “conflict,” consolidating clusters meticulously into variants, integrating PE and SR calls precisely, dynamically calculating PE and SR SV-support thresholds, retaining all clusters and choosing variants based on final support, and corroborating SVs using streamlined local read-depth information.

### METHODOLOGY

SVXplorer first forms discordant clusters from paired-end reads via formation of maximal cliques in a weight-thresholded bidirectional graph and consolidates them further into PE-supported variants. It then integrates split reads and read-depth information to call putative variants, enhancing/filtering out existing variants or identifying new ones along the way. 

### READ-GROUPS

SVXplorer should be run on a BAM file with a single read-group. For datasets with multiple read-groups, we recommend running SVXplorer on every read-group separately. We will add support for consolidation of these results soon. Please follow issue #81 for further details.

### REQUIREMENTS

SVXplorer should run on any Unix-based OS with bash, python > 2.6 and libraries as specified in the requirements file. "bedtools" and "samtools" executables should be on the user PATH. In addition, if a split read file is not available in the typical splitters format (2 entries per query name with 2 distinct, split queries) a script is provided to extract this from the alignment file using LUMPY's extractBwaMem_reads script (https://raw.githubusercontent.com/arq5x/lumpy-sv/master/scripts/extractSplitReads_BwaMem).

### LOCAL INSTALLATION

Download latest SVXplorer release from GitHub and unzip the directory. Alternatively, clone the repository. Then ensure that "samtools" and "bedtools" are on the PATH by running 

```
which samtools
which bedtools
``` 

Then install all the python libraries and install SVXplorer using

```
cd SVXplorer*
pip install -r requirements.txt
python setup.py bdist_wheel
pip install -e .
```

### USAGE

After installation SVXplorer and all its helper scripts are added to path and can be run as a command-line applications.

```
SVXplorer -h
```

will show all options available to the user. Then run

```
./run_test.sh
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

SVXplorer discordant.bam splitters.bam sample.bam reference.fa -i exclude.bed -c ignore_CHR.txt -m non_repeat_regions.bed -w pathToWorkingDirectory

Option -m expects a file listing regions not containing frequently repeated sequences, for use in assessing coverage, and SVC provides one for b37 in the data folder in zipped form, which will need to be unzipped before use.

### RESULTS

The final BEDPE and VCF files are written in the results folder. Intermediate BEDPE files using only PE/PE and SR mappings are stored in the "workspace" directory.

