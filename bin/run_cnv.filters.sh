CNVNATOR=/m/cphg-RLscratch/cphg-RLscratch/share/CNVnator_v0.3/src/cnvnator
REF=/home/kk7t/scratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37_decoy.fasta #/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta
BAM=$2 #/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim3/sim100.all.sorted.bam
#BAM=/m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/${1}/alignments/sample.bam #/home/kk7t/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/sample.na12878.bam
dir=../results/text/
rootfile=cnv.s2.${1}.root
sex=$4
GAPFILE=~/scratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/gaps.bed

#csplit --digits=2  --quiet --prefix=chr $REF "/>/" "{*}"
$CNVNATOR -genome $3 -root $dir/$rootfile -tree $BAM -unique
$CNVNATOR -root $dir/$rootfile -his 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -stat 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -partition 100 -d ~/store/REFs
$CNVNATOR -root $dir/$rootfile -call 100 -d ~/store/REFs > $dir/${1}.cnv

cat $dir/${1}.cnv \
        | awk -v OFS="\t" '$NF < 0.5 {split($2,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2],$1,$3,$4,$5,$6,$7,$8,$9}' \
        > $dir/${1}.qcnv

if [ ${sex} == "2" ]; then
    bedtools intersect \
	-a <(cat $dir/${1}.qcnv | awk '$1 != "Y"') \
	-b $GAPFILE -f 0.95 -v > $dir/${1}.nogaps.cnv
else
    bedtools intersect -a $dir/${1}.qcnv \
	-b $GAPFILE -f 0.95 -v > $dir/${1}.nogaps.cnv
fi

cat $dir/${1}.nogaps.cnv | cut -f1,2,3,6 > $dir/${1}.final.cnv
