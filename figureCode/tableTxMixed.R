# ==================================== 
# = table for transcript mixed model =
# ====================================

args = commandArgs(TRUE) # args <- c("../figureData/female.qtt.mixed.RData", "../figureData/male.qtt.mixed.RData")
library("openxlsx")

# female
load(args[1])
female.table <- rbind(c("Female", "Young", young.mix.h2),
											c("Female", "Old", old.mix.h2),
											c("Female", "Old - Young", diff.mix.h2))

load(args[2])
male.table <- rbind(c("Male", "Young", young.mix.h2),
										c("Male", "Old", old.mix.h2),
 										c("Male", "Old - Young", diff.mix.h2))
											
											
# table

final.table <- rbind(female.table, male.table)

# save table
# ============================================================

write.xlsx(list(table = final.table), file = args[3], colNames = FALSE, rowNames = FALSE)
