# lift to generate line specific genome and gtf
# ============================================================
mkdir -p /mnt/gs18/scratch/users/tansuxu/dgrp/lift/log
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift

while read LINE
do 
  line=`echo $LINE | awk '{print $2}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,vcf=/mnt/gs18/scratch/users/tansuxu/dgrp/jgil/all.jgil.hpp.vcf,gtf=/mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/stringtie.merge/r6.36.plus.candidate.gtf,line=$line --output=log/$line.lift.out --error=log/$line.lift.err /mnt/research/qgg/dgrp-age-exp/sbatch/lift.sbatch
done < /mnt/gs18/scratch/users/tansuxu/dgrp/bam.file.list

# hisat2 build
# ============================================================
cut -d " " -f 2 /mnt/gs18/scratch/users/tansuxu/dgrp/bam.file.list | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 1 -j 50 "
bash /mnt/gs18/scratch/users/tansuxu/dgrp/lift/hisat2.build.bash --resource /mnt/research/qgg/dgrp-age-exp --wdir /mnt/gs18/scratch/users/tansuxu/dgrp/lift --line {} > /mnt/gs18/scratch/users/tansuxu/dgrp/lift/{}/{}.hisat2.build.out 2>&1" > log/hisat2.build.parallel.02222021.log &

# get info of RNA-seq samples, hisat2 mapping
# ============================================================
mkdir -p /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/log
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2

awk '$2=="baseline" {print $1"\t"$3"\t"$3"_"$4 $5"\t"$6"\t"$7"\t"$10}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 6 -j 50 "
bash /mnt/gs18/scratch/users/tansuxu/dgrp/lift/hisat2.mapping.bash --resource /mnt/research/qgg/dgrp-age-exp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --wdir /mnt/gs18/scratch/users/tansuxu/dgrp/lift --out /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2 --sample {1} --line {2} --individual {3} --fc {4} --lane {5} --fastq {6} > {1}.out 2>&1" > log/hisat2.mapping.parallel.02222021.log &
# sacct: error: slurm_persist_conn_open_without_init: failed to open persistent connection to slurm-db-00:6819: Connection refused
# then use sbatch to get jobID

# merge bam files for each individual with more than one fastq.gz files, such as line 100 female 1 (100_F1)
# ============================================================
awk '$2=="baseline" {print $3"_"$4 $5"\t"$0}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis > final.samples.to.begin.analysis.baseline

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -i final.samples.to.begin.analysis.baseline|awk '/,/' > BaseLine.file.d.list

awk '$2=="baseline" {print $3"_"$4 $5}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | sort | uniq -u > BaseLine.file.u.list
ls /mnt/gs18/scratch/users/tansuxu/tmp > tmpList

cat BaseLine.file.d.list | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 50 "
bash /mnt/gs18/scratch/users/tansuxu/dgrp/lift/mergeBamInd.bash --resource /mnt/research/qgg/dgrp-age-exp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --out /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2 --list /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/tmpList --individual {1} --bamInd {2} > {1}.out 2>&1" > log/mergeBamInd.parallel.02222021.log &

# move bam files for individual with one fastq.gz file
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2
 
grep -w -f BaseLine.file.u.list final.samples.to.begin.analysis.baseline | awk '{print $1"\t"$2}' > BaseLine.uniq

cat BaseLine.uniq | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 20 " mkdir {1} && mv /mnt/gs18/scratch/users/tansuxu/tmp/{2}/map.bam {1}/sort.bam && rm -r /mnt/gs18/scratch/users/tansuxu/tmp/{2} && /mnt/research/qgg/software/samtools-1.10/samtools index {1}/sort.bam  > {1}.out 2>&1" > log/mv.BaseLine.parallel.02222021.log &

# use stringtie to obtain gene expression
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2
awk '$2=="baseline" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | sort | uniq | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 50 " bash /mnt/research/qgg/dgrp-age-exp/bash/stringtie.bash --resource /mnt/research/qgg/dgrp-age-exp --wdir /mnt/gs18/scratch/users/tansuxu/dgrp/lift/ --bam /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/{1}/sort.bam --out /mnt/research/qgg/dgrp-age-exp/exp/baseline --individual {1} --line {2} > /mnt/research/qgg/dgrp-age-exp/exp/baseline/log/{1}.stringtie.out 2>&1" > /mnt/research/qgg/dgrp-age-exp/exp/baseline/log/stringtie.baseline.parallel.05172021.log &

# 789 individuals in baseline
# check if all successfully complete
for sample in `awk '$2=="baseline" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | cut -f 1 | sort | uniq`
do
  if [ `wc -l /mnt/research/qgg/dgrp-age-exp/exp/baseline/$sample/stringtie.err | awk '{print $1}'` -ne 0 ]
  then
    echo $sample
  fi
done


# do similar things after hisat2 build for 3WK fastq 
# after making sure BaseLine pipeline is ok

# get info of RNA-seq samples, hisat2 mapping
# ============================================================
mkdir -p /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/log
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2

awk '$2=="3wk" {print $1"\t"$3"\t"$3"_"$4 $5"\t"$6"\t"$7"\t"$10}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 6 -j 50 "bash /mnt/gs18/scratch/users/tansuxu/dgrp/lift/hisat2.mapping.bash --resource /mnt/research/qgg/dgrp-age-exp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --wdir /mnt/gs18/scratch/users/tansuxu/dgrp/lift --out /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2 --sample {1} --line {2} --individual {3} --fc {4} --lane {5} --fastq {6} > {1}.out 2>&1" > log/hisat2.mapping.3w.parallel.02222021.log &

# merge bam files for each individual with more than one fastq.gz 
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2

awk '$2=="3wk" {print $3"_"$4 $5"\t"$0}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis > final.samples.to.begin.analysis.3wk

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -i final.samples.to.begin.analysis.3wk | awk '/,/' > 3wk.file.d.list
awk '$2=="3wk" {print $3"_"$4 $5}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | sort | uniq -u > 3wk.file.u.list
ls /mnt/gs18/scratch/users/tansuxu/tmp > tmpList

cat 3wk.file.d.list | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 50 "
bash /mnt/gs18/scratch/users/tansuxu/dgrp/lift/mergeBamInd.bash --resource /mnt/research/qgg/dgrp-age-exp --tmp /mnt/gs18/scratch/users/tansuxu/tmp --out /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2 --list /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/tmpList --individual {1} --bamInd {2} > {1}.out 2>&1" > log/mergeBamInd.3wk.parallel.02222021.log &

# move bam files for individual with one fastq.gz file
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2
 
grep -w -f 3wk.file.u.list final.samples.to.begin.analysis.3wk | awk '{print $1"\t"$2}' > 3wk.uniq

cat 3wk.uniq | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 20 " mkdir {1} && mv /mnt/gs18/scratch/users/tansuxu/tmp/{2}/map.bam {1}/sort.bam && rm -r /mnt/gs18/scratch/users/tansuxu/tmp/{2} && /mnt/research/qgg/software/samtools-1.10/samtools index {1}/sort.bam > {1}.out 2>&1" > log/mv.3wk.parallel.02222021.log &

# use stringtie to obtain gene expression
# ============================================================
cd /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2
awk '$2=="3wk" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | sort | uniq | tr '\t' '\n' | /mnt/research/qgg/software/parallel-20200722/bin/parallel -N 2 -j 50 " bash /mnt/research/qgg/dgrp-age-exp/bash/stringtie.bash --resource /mnt/research/qgg/dgrp-age-exp --wdir /mnt/gs18/scratch/users/tansuxu/dgrp/lift/ --bam /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/{1}/sort.bam --out /mnt/research/qgg/dgrp-age-exp/exp/3wk --individual {1} --line {2} > /mnt/research/qgg/dgrp-age-exp/exp/3wk/log/{1}.stringtie.out 2>&1" > /mnt/research/qgg/dgrp-age-exp/exp/3wk/log/stringtie.baseline.parallel.05192021.log &

# 807 individuals in 3wk
# check if complete
for sample in `awk '$2=="3wk" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | cut -f 1 | sort | uniq`
do
  if [ `wc -l /mnt/research/qgg/dgrp-age-exp/exp/3wk/$sample/stringtie.err | awk '{print $1}'` -ne 0 ]
  then
    echo $sample
  fi
done


# combine gene count
# ============================================================

# prepare sorted gene list
cut -f 1 /mnt/research/qgg/dgrp-age-exp/exp/*/*/gene_abund.tab | grep -v Gene | sort | uniq > /mnt/research/qgg/dgrp-age-exp/exp/uniq.gene.id

# get tpm
mkdir /mnt/research/qgg/dgrp-age-exp/exp/join/
for sample in `awk '$2=="baseline" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | cut -f 1 | sort | uniq`
do
  tail -n+2 /mnt/research/qgg/dgrp-age-exp/exp/baseline/$sample/gene_abund.tab | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 9 -o sum | join -a 1 -o '2.2' -e 'NA' /mnt/research/qgg/dgrp-age-exp/exp/uniq.gene.id - | cat <(echo "baseline_$sample") - > /mnt/research/qgg/dgrp-age-exp/exp/join/baseline.$sample.tpm
done &

for sample in `awk '$2=="3wk" {print $3"_"$4 $5"\t"$3}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | cut -f 1 | sort | uniq`
do
  tail -n+2 /mnt/research/qgg/dgrp-age-exp/exp/3wk/$sample/gene_abund.tab | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 9 -o sum | join -a 1 -o '2.2' -e 'NA' /mnt/research/qgg/dgrp-age-exp/exp/uniq.gene.id - | cat <(echo "3wk_$sample") - > /mnt/research/qgg/dgrp-age-exp/exp/join/3wk.$sample.tpm
done &

# join expression, separate by female and male
cat <(echo "gene") /mnt/research/qgg/dgrp-age-exp/exp/uniq.gene.id | paste - /mnt/research/qgg/dgrp-age-exp/exp/join/*.*F*.tpm > /mnt/research/qgg/dgrp-age-exp/joinExp/female.exp.tpm

cat <(echo "gene") /mnt/research/qgg/dgrp-age-exp/exp/uniq.gene.id | paste - /mnt/research/qgg/dgrp-age-exp/exp/join/*.*M*.tpm > /mnt/research/qgg/dgrp-age-exp/joinExp/male.exp.tpm

# summarize alignment 
# ============================================================
while read line
do
  age=`echo $line | awk '{print $2}'`
  dir=`echo $line | awk '{print $1}'`
  if [[ $age = "baseline" ]]
  then
    cat /mnt/gs18/scratch/users/tansuxu/dgrp/lift/BaseHisat2/$dir/log/map.summary | awk '{print $1}' | paste - - - - - - | cut -f 1,3,4,5 | awk -v ll="$line" '{print $0"\t"ll}'
  else
    cat /mnt/gs18/scratch/users/tansuxu/dgrp/lift/3wkHisat2/$dir/log/map.summary | awk '{print $1}' | paste - - - - - - | cut -f 1,3,4,5 | awk -v ll="$line" '{print $0"\t"ll}'
  fi
done < /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis > /mnt/research/qgg/dgrp-age-exp/figureData/sample.map.summary &
# first four columns are 1. total number of reads; 2. reads mapped 0 times; 3. reads
# mapped 1 time and 4. reads mapped > 1 times (non-unique)

