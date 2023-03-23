# ================
# = prepare data =
# ================

# fastq's downloaded from hyperion 
# baseline: /home2/hiseq2500/FastqFiles/DGRP_BaseLine
# 3wk: /home2/hiseq2500/FastqFiles/DGRP_3WK
# ============================================================

# get SampleSheet.csv information from hiseq2500
# this helps find information on the runs, including dates, flow cells, lanes, etc.
# they are also available in the fastq's themselves
# ============================================================

for file in `find /home2/hiseq2500/ -maxdepth 3 -name "SampleSheet.csv"`
do
  target=`echo $file | sed 's/\//ddd/g'`
  cp $file $target
done &
  
