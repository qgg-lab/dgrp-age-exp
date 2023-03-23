# ===================================
# = table for eigengene association =
# ===================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.phototaxis.decline.qtt.eigen.assoc.RData", "../figureData/female.fecundity.decline.qtt.eigen.assoc.RData", "../figureData/female.lifespan.qtt.eigen.assoc.RData", "../figureData/male.phototaxis.decline.qtt.eigen.assoc.RData", "../figureData/male.speed.decline.qtt.eigen.assoc.RData", "../figureData/male.endurance.decline.qtt.eigen.assoc.RData", "../figureData/male.lifespan.qtt.eigen.assoc.RData")

library("openxlsx")

# read and process data
# ============================================================

load(args[1])
female.table <- rbind(cbind("Female", "Young", "Decline in phototaxis", module.assoc.young),
                      cbind("Female", "Aged", "Decline in phototaxis", module.assoc.old),
                      cbind("Female", "Aged - Young", "Decline in phototaxis", module.assoc.diff))


load(args[2])
female.table <- rbind(female.table,
											cbind("Female", "Young", "Decline in fecundity", module.assoc.young),
											cbind("Female", "Aged", "Decline in fecundity", module.assoc.old),
											cbind("Female", "Aged - Young", "Decline in fecundity", module.assoc.diff))


load(args[3])
female.table <- rbind(female.table,
											cbind("Female", "Young", "Lifespan", module.assoc.young),
											cbind("Female", "Aged", "Lifespan", module.assoc.old),
											cbind("Female", "Aged - Young", "Lifespan", module.assoc.diff))



load(args[4])
male.table <- rbind(cbind("Male", "Young", "Decline in phototaxis", module.assoc.young),
										cbind("Male", "Aged", "Decline in phototaxis", module.assoc.old),
										cbind("Male", "Aged - Young", "Decline in phototaxis", module.assoc.diff))


load(args[5])
male.table <- rbind(male.table,
										cbind("Male", "Young", "Decline in climbing speed", module.assoc.young),
										cbind("Male", "Aged", "Decline in climbing speed", module.assoc.old),
										cbind("Male", "Aged - Young", "Decline in climbing speed", module.assoc.diff))

load(args[6])
male.table <- rbind(male.table,
										cbind("Male", "Young", "Decline in climbing endurance", module.assoc.young),
										cbind("Male", "Aged", "Decline in climbing endurance", module.assoc.old),
										cbind("Male", "Aged - Young", "Decline in climbing endurance", module.assoc.diff))


load(args[7])
male.table <- rbind(male.table,
										cbind("Male", "Young", "Lifespan", module.assoc.young),
										cbind("Male", "Aged", "Lifespan", module.assoc.old),
										cbind("Male", "Aged - Young", "Lifespan", module.assoc.diff))

# final.table
# ============================================================

final.table <- rbind(female.table, male.table)

# write excel file
# ============================================================

write.xlsx(list('table' = final.table), file = args[8], colNames = TRUE, rowNames = FALSE)
