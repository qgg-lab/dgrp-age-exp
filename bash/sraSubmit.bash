# ==================================================
# = upload DGRP 3WK data from MSU HPCC to NCBI SRA =
# ==================================================

# Bioproject description, may need changes later

# Aging is an inevitable biological process that could lead to physiological senescence and diseases. The underlying mechanisms are still being elucidated. Drosophila is an excellent model for aging studies with many distinct advantages. Here, we obtained the transcriptomic data of 200 sequenced Drosophila Genetic Reference Panel (DGRP) inbred lines at 3-week old and combined them with those from 3-5 day old flies from the same DGRP lines (Everett et al., 2020). We aim to evaluate the age effect on gene expression and uncover the genetic mechanisms of aging by performing a series of analyses, including quantitative genetics, genotype by age interaction, expression QTL mapping, mediation analysis, and network construction. This study would provide new insight into the aging process and pave the road for future functional analyses.

# design_description in SRA_metadata.xlsx, may need changes later

# Flies were maintained at the standard 25 °C, 60-75% relative humidity, and 12-hour light-dark cycle condition on cornmeal-molasses-agar medium. At 3 weeks post eclosion, two biological replicates of 25 females or 30 males were collected for each DGRP line and stored at -80 °C. Total RNA was extracted and subjected to ribosomal RNA (rRNA) depletion before being converted to a strand-specific cDNA libraries using a dUTP based protocol (Everett et al., 2020). The cDNA libraries were barcoded and randomly multiplexed into 50 pools of 16 libraries each. Each pool was sequenced on an Illumina HiSeq 2500 with 125bp single end reads. 

# prepare file for Bioproject and Biosample, then edit and upload using Excel

awk '$2=="3wk" {print "DGRP_Line_"$3"_3wk_"$4"_Rep_"$5" DGRP_"$3" "$4" "$10}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis > ~/sraSubmit/dgrp.sra.txt

grep '3wk' /mnt/ufs18/rs-015/qgg/dgrp-age-exp/normalize/bad.sample.id | sed 's/F/F_/g' | sed 's/M/M_/g' | awk -F _ '{print "DGRP_Line_"$2"_3wk_"$3"_Rep_"$4}' > ~/sraSubmit/bad.id

grep -f ~/sraSubmit/bad.id -v ~/sraSubmit/dgrp.sra.txt > ~/sraSubmit/dgrp.sra.submit #wc -l 796

sed 's/ /\t/g' ~/sraSubmit/dgrp.sra.submit | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 4 -o collapse -i - > ~/sraSubmit/dgrp.sra.submit.collapse #wc -l 794

# creat softlink and transfer data to sra

awk '$2=="3wk" {print $10}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis| sed 's/\/mnt\/ufs18\/scratch/\/mnt\/gs18\/scratch\/users/g'|xargs -I{} ln -s {} /mnt/home/tansuxu/sra.submision

module load Aspera-Connect/3.9.6
ascp -i /mnt/home/tansuxu/aspera.openssh -QT -l100m -k1 -d /mnt/home/tansuxu/sra.submision subasp@upload.ncbi.nlm.nih.gov:uploads/tansuxu_msu.edu_VUq6CPJ8

 