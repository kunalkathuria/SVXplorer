CNVNATOR=/m/cphg-RLscratch/cphg-RLscratch/share/CNVnator_v0.3/src/cnvnator
REF=/home/kk7t/scratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37_decoy.fasta #/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta
BAM=/m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam #/home/kk7t/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/sample.na12878.bam
dir=../results/text/
rootfile=cnv_out.na12878.srr505885.root

#csplit --digits=2  --quiet --prefix=chr $REF "/>/" "{*}"
$CNVNATOR -root $dir/$rootfile -tree $BAM
$CNVNATOR -root $dir/$rootfile -his 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -stat 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -partition 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -call 100 -d ~/store/REFs > $dir/cnvs.txt

