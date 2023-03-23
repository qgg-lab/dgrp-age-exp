# ==================================
# = table for network preservation =
# ==================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.wgcna.RData", "../figureData/male.wgcna.RData")

library("openxlsx")

# read and process data
# ============================================================

load(args[1])
female.mod.pres <- mod.pres

load(args[2])
male.mod.pres <- mod.pres

# final.table
# ============================================================

final.table <- rbind(cbind("female", 1:nrow(female.mod.pres$observed), female.mod.pres$observed, female.mod.pres$p.values),
										 cbind("male", 1:nrow(male.mod.pres$observed), male.mod.pres$observed, male.mod.pres$p.values))

# write excel file
# ============================================================

write.xlsx(list('table' = final.table), file = args[3], colNames = TRUE, rowNames = FALSE)
