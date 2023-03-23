# ============================================================
# prepare GTF file that include novel transcribed region
# ============================================================

mkdir -p /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/log
mkdir /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/code

# randomly sample fastq from final.samples.to.begin.analysis
# seprate into conditions (age + sex combination)
# from each fastq, sample 0.5M reads
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF
mkdir tmp.baseline.f
mkdir tmp.baseline.m
mkdir tmp.3wk.f
mkdir tmp.3wk.m

while read line
do
  sample=`echo $line | awk '{print $1}'`
  fastq=`echo $line | awk '{print $10}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=$fastq,sample=$sample,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.f,mem=16G,head=500000 --out=log/sampleFastq.baseline.f.$sample.out --error=log/sampleFastq.baseline.f.$sample.err code/sampleFastq.sbatch
done < <(awk '$2=="baseline" && $4=="F"' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis)


while read line
do
  sample=`echo $line | awk '{print $1}'`
  fastq=`echo $line | awk '{print $10}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=$fastq,sample=$sample,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.m,mem=16G,head=500000 --out=log/sampleFastq.baseline.m.$sample.out --error=log/sampleFastq.baseline.m.$sample.err code/sampleFastq.sbatch
done < <(awk '$2=="baseline" && $4=="M"' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis)


while read line
do
  sample=`echo $line | awk '{print $1}'`
  fastq=`echo $line | awk '{print $10}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=$fastq,sample=$sample,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.f,mem=16G,head=500000 --out=log/sampleFastq.3wk.f.$sample.out --error=log/sampleFastq.3wk.f.$sample.err code/sampleFastq.sbatch
done < <(awk '$2=="3wk" && $4=="F"' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis)


while read line
do
  sample=`echo $line | awk '{print $1}'`
  fastq=`echo $line | awk '{print $10}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=$fastq,sample=$sample,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.m,mem=16G,head=500000 --out=log/sampleFastq.3wk.m.$sample.out --error=log/sampleFastq.3wk.m.$sample.err code/sampleFastq.sbatch
done < <(awk '$2=="3wk" && $4=="M"' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis)


cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.f && ls |wc -l #474
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.m && ls |wc -l #472
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.f && ls |wc -l #404
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.m && ls |wc -l #408

cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.f && cat * > /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.baseF.fastq
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.baseline.m && cat * > /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.baseM.fastq
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.f && cat * > /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.3F.fastq
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp.3wk.m && cat * > /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.3M.fastq

rm -r /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/tmp*

# hisat2 alignment
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF
mkdir -p baseF/log
mkdir -p baseM/log
mkdir -p 3F/log
mkdir -p 3M/log

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/baseF,fastq=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.baseF.fastq --out=log/sample.baseF.hisat2.out --error=log/sample.baseF.hisat2.err code/sample.hisat2.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/baseM,fastq=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.baseM.fastq --out=log/sample.baseM.hisat2.out --error=log/sample.baseM.hisat2.err code/sample.hisat2.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/3F,fastq=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.3F.fastq --out=log/sample.3F.hisat2.out --error=log/sample.3F.hisat2.err code/sample.hisat2.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/3M,fastq=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/pooled.sample.3M.fastq --out=log/sample.3M.hisat2.out --error=log/sample.3M.hisat2.err code/sample.hisat2.sbatch

# stringtie
# ============================================================
sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/baseF,gtf=/mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all.gtf --out=log/sample.baseF.stringtie.out --error=log/sample.baseF.stringtie.err code/sample.stringtie.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/baseM,gtf=/mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all.gtf --out=log/sample.baseM.stringtie.out --error=log/sample.baseM.stringtie.err code/sample.stringtie.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/3F,gtf=/mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all.gtf --out=log/sample.3F.stringtie.out --error=log/sample.3F.stringtie.err code/sample.stringtie.sbatch

sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,dir=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/3M,gtf=/mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all.gtf --out=log/sample.3M.stringtie.out --error=log/sample.3M.stringtie.err code/sample.stringtie.sbatch

# stringtie merge
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF
mkdir stringtie.merge
/mnt/research/qgg/software/stringtie-1.3.6.Linux_x86_64/stringtie --merge -p 8 -G dmel-all.gtf -o stringtie.merge/stringtie.merged.gtf stringtie.merge/merge.txt > log/stringtie.merge.log 2>&1

# stringtie estimate expression
# ============================================================

/mnt/research/qgg/software/stringtie-1.3.6.Linux_x86_64/stringtie -e -B -p 8 -G stringtie.merge/stringtie.merged.gtf -o baseF/baseF.eB.gtf -A baseF/gene_abund.tab baseF/hisat2.map.bam > log/stringtie.eB.baseF.log 2>&1 &

/mnt/research/qgg/software/stringtie-1.3.6.Linux_x86_64/stringtie -e -B -p 8 -G stringtie.merge/stringtie.merged.gtf -o baseM/baseM.eB.gtf -A baseM/gene_abund.tab baseM/hisat2.map.bam > log/stringtie.eB.baseM.log 2>&1 &

/mnt/research/qgg/software/stringtie-1.3.6.Linux_x86_64/stringtie -e -B -p 8 -G stringtie.merge/stringtie.merged.gtf -o 3F/3F.eB.gtf -A 3F/gene_abund.tab 3F/hisat2.map.bam > log/stringtie.eB.3F.log 2>&1 &

/mnt/research/qgg/software/stringtie-1.3.6.Linux_x86_64/stringtie -e -B -p 8 -G stringtie.merge/stringtie.merged.gtf -o 3M/3M.eB.gtf -A 3M/gene_abund.tab 3M/hisat2.map.bam > log/stringtie.eB.3M.log 2>&1 &

# get all transcripts expression in each condition to find cutoff
# their were some transcripts split by stringtie
# e.g. FBtr0084079
# but this is OK, later on the transcripts are pulled
# from the original annotaitons
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF
awk '$3=="transcript" {print $12"\t"$18}' baseF/baseF.eB.gtf | tr -d '"' | tr -d ';' | sort -k1,1 | join -t $'\t' - <(awk '$3=="transcript" {print $12"\t"$18}' baseM/baseM.eB.gtf | tr -d '"' | tr -d ';' | sort -k1,1) | join -t $'\t' - <(awk '$3=="transcript" {print $12"\t"$18}' 3F/3F.eB.gtf | tr -d '"' | tr -d ';' | sort -k1,1) | join -t $'\t' - <(awk '$3=="transcript" {print $12"\t"$18}' 3M/3M.eB.gtf | tr -d '"' | tr -d ';' | sort -k1,1) > /mnt/research/qgg/dgrp-age-exp/figureData/rnaseq.stringtie.tx.exp

# get transcrips with expression > 1TPM in at least one condition
# ============================================================
awk '$3=="transcript" {print $18"\t"$12}' baseF/baseF.eB.gtf | tr -d '"' | tr -d ';' | awk '$1>1 {print $2}' > tpm1.transcript
awk '$3=="transcript" {print $18"\t"$12}' baseM/baseM.eB.gtf | tr -d '"' | tr -d ';' | awk '$1>1 {print $2}' >> tpm1.transcript
awk '$3=="transcript" {print $18"\t"$12}' 3F/3F.eB.gtf | tr -d '"' | tr -d ';' | awk '$1>1 {print $2}' >> tpm1.transcript
awk '$3=="transcript" {print $18"\t"$12}' 3M/3M.eB.gtf | tr -d '"' | tr -d ';' | awk '$1>1 {print $2}' >> tpm1.transcript

sort tpm1.transcript | uniq > tpm1.transcript.uniq

# Run cuffcompare
# ============================================================
/mnt/research/qgg/software/cufflinks-2.2.1.Linux_x86_64/cuffcompare -o stringtie.merge/comp -r dmel-all.gtf -s /mnt/research/qgg/resource/flybase/r6.36/genome/genome.fa stringtie.merge/stringtie.merged.gtf > log/compare.stringtie.merge.log 2>&1

# Extract candidate transcripts
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF
awk '$3 == "exon" {print $0" class_code \"=\";"}' dmel-all.gtf > stringtie.merge/r6.36.plus.candidate.gtf
perl extract_ju.pl stringtie.merge/comp.combined.gtf > extract_ju.gtf
awk '/"u"/ {print $16}' extract_ju.gtf | tr -d '"' | tr -d ';' | paste - <(awk '/"u"/' extract_ju.gtf) | sort -k1,1 | join -1 1 -2 1 -t $'\t' - tpm1.transcript.uniq | cut -f 2- >> stringtie.merge/r6.36.plus.candidate.gtf
awk '/"j"/ {print $18}' extract_ju.gtf | tr -d '"' | tr -d ';' | paste - <(awk '/"j"/' extract_ju.gtf) | sort -k1,1 | join -1 1 -2 1 -t $'\t' - tpm1.transcript.uniq | cut -f 2- >> stringtie.merge/r6.36.plus.candidate.gtf

# summarize transcripts
# ============================================================

# number of genes + NTRs
cat /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.gtf | perl -wne 'chomp $_; if ($_ =~ m/gene_id \"(.*?)\";/) { print $1, "\n"; }' | uniq | sort | uniq | wc -l
# transcripts
cat /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.gtf | perl -wne 'chomp $_; if ($_ =~ m/transcript_id \"(.*?)\";/) { print $1, "\n"; }' | uniq | sort | uniq | wc -l
# NTR
cat /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.gtf | perl -wne 'chomp $_; if ($_ =~ m/gene_id \"(.*?)\";/) { print $1, "\n"; }' | uniq | sort | uniq | grep -v FBgn | wc -l
# NTR transcripts
cat /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.gtf | perl -wne 'chomp $_; if ($_ =~ m/gene_id \"(.*?)\";.*transcript_id \"(.*?)\";/) { print $1, "\t",$2, "\n"; }' | awk '!($1 ~ /FBgn/) {print $2}' | uniq | sort | uniq | wc -l

# also prepare a GTF file for estimating gene expression
# ============================================================

# remove rRNA gene; 2. remove 7SLRNA genes 
awk '$2 ~ /rRNA/ || $2 ~ /7SLRNA/ {print $1}' /mnt/research/qgg/dgrp-age-exp/figureData/gene.info | sort | uniq > /mnt/research/qgg/dgrp-age-exp/gtf/rRNA.7SLRNA.gene.id

awk -F "\"" '{print $2"\t"$0}' /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.gtf | sort -k1,1 | join -v 1 -t $'\t' - /mnt/research/qgg/dgrp-age-exp/gtf/rRNA.7SLRNA.gene.id | cut -f 2- | sort -k1,1 -k2,2n > /mnt/research/qgg/dgrp-age-exp/gtf/r6.36.plus.candidate.no.rRNA.no.7SLRNA.gtf





