# prepare file information to filter samples based on distance
# between RNA-Seq and DNA
# also assess among the RNA-Seq samples, whether
# all samples deviate from the DNA-Seq or only a few
# get information on min, median, max
# ============================================================

# prepare info file for baseline samples
# ============================================================

for file in `find /mnt/research/qgg/dgrp-age-exp/checkID -type f -name "*.dist.txt" | grep -v 3WKDGRP`
do
  location=`echo $file | awk -F '/' '{print "/mnt/ufs18/scratch/huangw53/dgrp-rna/DGRP_BaseLine/"$NF}' | sed 's/\.dist\.txt/\.fastq\.gz/'`
  sample=`echo $file | sed 's/\/mnt\/research\/qgg\/dgrp-age-exp\/checkID\///' | sed 's/\.dist\.txt//' | sed 's/R1_//' | awk -F '[-_]' '{print $1"_"$2"_"$4"_"$5}'`
  line=`echo $sample | awk -F '_' '{print $1}'`
  sex=`echo $sample | awk -F '_' '{print $2}'| awk -F '[0-9]' '{print $1}'`
  number=`echo $sample | awk -F '_' '{print $2}' | sed -E -e 's/F|M//' ` 
  fc=`echo $sample | awk -F '_' '{print $4}'`
  lane=`echo $sample | awk -F '_' '{print $3}'`
  leastID=`awk '{print $2"\t"$1}' $file | sort -k1,1n | head -1 | awk '{print $2}'`
  leastDist=`awk '{print $2"\t"$1}' $file | sort -k1,1n | head -1 | awk '{print $1}'`
  echo $sample "baseline" $line $sex $number $fc $lane $leastID $leastDist $location
done > /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.info

# prepare info file for 3WK samples
# ============================================================

for file in `find /mnt/research/qgg/dgrp-age-exp/checkID -type f -name "*.dist.txt" | grep 3WKDGRP`
do
  location=`echo $file | awk -F '/' '{print "/mnt/ufs18/scratch/huangw53/DGRP_3WK/"$NF}' | sed 's/\.dist\.txt/\.fastq\.gz/'`
  sample=`echo $file | sed 's/\/mnt\/research\/qgg\/dgrp-age-exp\/checkID\///' | sed 's/\.dist\.txt//' | sed 's/F_/F/' | sed 's/M_/M/' | awk -F '[-_]' '{print $2"_"$3"_"$5"_"$7}' | tr [:lower:] [:upper:]`
  line=`echo $sample | awk -F '_' '{print $1}'`
  sex=`echo $sample | awk -F '_' '{print $2}'| awk -F '[0-9]' '{print $1}'`
  number=`echo $sample | awk -F '_' '{print $2}' | sed -E -e 's/F|M//' ` 
  fc=`echo $sample | awk -F '_' '{print $4}'`
  lane=`echo $sample | awk -F '_' '{print $3}'`
  leastID=`awk '{print $2"\t"$1}' $file | sort -k1,1n | head -1 | awk '{print $2}'`
  leastDist=`awk '{print $2"\t"$1}' $file | sort -k1,1n | head -1 | awk '{print $1}'`
  echo $sample "3wk" $line $sex $number $fc $lane $leastID $leastDist $location
done >> /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.info

# summarize for each line
# get rid of the flow cell C7TRFANXX, all of them have
# distance above 0.35
# even if the RNA-Seq sample's minimum distance agrees with 
# DNA ID, they cannot be retained.
# ============================================================

awk '$6 != "C7TRFANXX"' /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.info | sed 's/ /\t/g' | awk '{print $8"\t"$0}' | sort -k1,1 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,10,10,10,10 -o collapse,collapse,min,max,median > /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.summary.stat

# take a few steps to identify samples to filtered out
# ============================================================

# 1. get rid of lines completely when the median distance is > 0.12, why 0.12?
#    let's assume 0.02 error rate + 0.10 difference to be tolorated.
#    this removes 4 lines, > 0.10 removes 10 lines
#    this difference is a risk that we are willing to take
#    it is arbitrary, indeed.

awk '$6 > 0.12 {print $2}' /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.summary.stat | sed 's/,/\n/g' > /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.samples.to.remove

# 2. get rid of samples that are > 0.10 && more than 0.05 above the median within
#    that line. For example, most of the 313 samples
#    have dist < 0.05 for the baseline, but < 0.10 for the 3wk samples
#    3wk samples are consistently more distant than the baseline samples
#    but there is one 3wek sample that is 0.123 and > 0.2 from the median

cat /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.summary.stat | perl -wne 'chomp $_; @line = split /\t/, $_; @samples = split /,/, $line[1]; @dist = split /,/, $line[2]; $min = $line[3]; $max = $line[4]; $median = $line[5]; for (my $i = 0; $i <= $#samples; $i++) {print $samples[$i], "\t", $dist[$i], "\t", $min, "\t", $max, "\t", $median, "\n"; }' | awk '$2 - $5 > 0.05 && $2 > 0.10 {print $1}' >> /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.samples.to.remove

# 3. finally, get rid of these samples and also check
#    if the least dist ID is the same from the claimed ID
#    in file names
# ============================================================

awk '$6 != "C7TRFANXX"' /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.dist.info | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' -v 1 - <(sort /mnt/research/qgg/dgrp-age-exp/checkID/rnaseq.samples.to.remove | uniq) | awk '$3 == $8' > /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis
