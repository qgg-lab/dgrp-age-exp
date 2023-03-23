# ======================================================
# = analysis of the dgrp baseline and 3wk RNA-Seq data =
# ======================================================

# 1. check RNA-Seq data against the DGRP sequences
# make sure sample IDs are right
# run steps in checkID.bash
# ============================================================

# 2. after check IDs
# make filters, prepare file locations, etc.
# run steps in prepare.info.bash
# ============================================================

# 3. prepare final GTF for analysis
# run steps in prepare.GTF.bash
# ============================================================

# 4. perform line-specific alignments, obtain gene expression data and merge together
# run steps in align.bash
# ============================================================

# 5. normalize, filter and adjust expression
# run steps in normalize.bash
# ============================================================


# 5. perform PCA and sva analysis
# run steps in PCA.bash
# ============================================================

# 6. perform qutatitative genetics analysis
# run steps in qg.bash
# ============================================================

# 7. perform GSEA analysis
# run steps in gsea.bash
# ============================================================



