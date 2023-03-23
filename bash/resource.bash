# ==============================
# = prepare resources for dgrp =
# ==============================

# 1. make hisat2 index in flybase
# ============================================================

/mnt/research/qgg/software/hisat2-2.2.1/hisat2_extract_exons.py $dir/lift.gtf > $dir/lift.exons
$HISAT2SS $dir/lift.gtf > $dir/lift.splice_sites

$HISAT2BD -p 8 $dir/lift.fa --ss $dir/lift.splice_sites --exon $dir/lift.exons $dir/lift.genome
