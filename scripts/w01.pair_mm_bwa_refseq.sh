reference=/data/Analysis/fanxiaoying/database/mm9/03.ucsc_ercc/ucsc.ercc.mRNA.fa
fq1=$2
fq2=$3
outdir=$1
bam_prefix=$4
stat=$5
summarize=$6
insert=$7

mkdir $outdir
bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq1 > $outdir/$insert.1.sai
bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq2 > $outdir/$insert.2.sai
bwa sampe $reference $outdir/$insert.1.sai  $outdir/$insert.2.sai  $fq1 $fq2 >$outdir/$insert.sam
samtools view  -u -b -S -t $reference.fai $outdir/$insert.sam | samtools sort -m 200000000 - $bam_prefix
samtools flagstat $bam_prefix* >$stat
samtools view -F 4 $bam_prefix* | awk '{print $3}' | perl -ne 'chomp;$hash{$_}++;END{foreach(keys %hash){print "$_\t$hash{$_}\n"}}' >$summarize

rm $outdir/$insert.1.sai
rm $outdir/$insert.2.sai
rm $outdir/$insert.sam



