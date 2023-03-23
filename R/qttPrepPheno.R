# ======================================
# = prepare phenotype for QTT analysis =
# ======================================

# 1. Durham paper, fecundity and lifespan

fecundity.wk1 <- read.csv("durhamFecundityWeek1.csv", header = F, as.is = T, na.strings = ".")[1:189, ]
fecundity.wk3 <- read.csv("durhamFecundityWeek3.csv", header = F, as.is = T, na.strings = ".")[1:189, ]
durham.lifespan <- read.csv("durhamLifespan.csv", header = F, as.is = T, na.strings = ".")[1:189, ]

fecundity.wk1and3 <- merge(merge(fecundity.wk1, fecundity.wk3, by.x = 1, by.y = 1), durham.lifespan, by.x = 1, by.y = 1)
fecundity.wk1and3 <- within(fecundity.wk1and3, decline <- V2.x - V2.y)
rownames(fecundity.wk1and3) <- fecundity.wk1and3[, 1]

# 2. Carbone paper, visual senescence

phototaxis <- read.csv("carbonePhototaxis.csv", header = F, as.is = T)
phototaxis[, 1] <- gsub("RAL", "line", phototaxis[, 1])

phototaxis.female.wk1 <- with(phototaxis[phototaxis[, 3] == "F" & phototaxis[, 4] == 4, ], sapply(split(V5, V1), mean))
phototaxis.female.wk1 <- data.frame(line = names(phototaxis.female.wk1), score = phototaxis.female.wk1)

phototaxis.female.wk4 <- with(phototaxis[phototaxis[, 3] == "F" & phototaxis[, 4] == 28, ], sapply(split(V5, V1), mean))
phototaxis.female.wk4 <- data.frame(line = names(phototaxis.female.wk4), score = phototaxis.female.wk4)

phototaxis.female.wk1and4 <- merge(phototaxis.female.wk1, phototaxis.female.wk4, by.x = 1, by.y = 1)
phototaxis.female.wk1and4 <- within(phototaxis.female.wk1and4, decline <- score.x - score.y)
rownames(phototaxis.female.wk1and4) <- phototaxis.female.wk1and4[, 1]

phototaxis.male.wk1 <- with(phototaxis[phototaxis[, 3] == "M" & phototaxis[, 4] == 4, ], sapply(split(V5, V1), mean))
phototaxis.male.wk1 <- data.frame(line = names(phototaxis.male.wk1), score = phototaxis.male.wk1)

phototaxis.male.wk4 <- with(phototaxis[phototaxis[, 3] == "M" & phototaxis[, 4] == 28, ], sapply(split(V5, V1), mean))
phototaxis.male.wk4 <- data.frame(line = names(phototaxis.male.wk4), score = phototaxis.male.wk4)

phototaxis.male.wk1and4 <- merge(phototaxis.male.wk1, phototaxis.male.wk4, by.x = 1, by.y = 1)
phototaxis.male.wk1and4 <- within(phototaxis.male.wk1and4, decline <- score.x - score.y)
rownames(phototaxis.male.wk1and4) <- phototaxis.male.wk1and4[, 1]

# 3. Gabrawy paper, activity

climbing.speed.wk1 <-  read.csv("gabrawyClimbingSpeedWeek1.csv", header = F, as.is = T)
climbing.speed.wk5 <-  read.csv("gabrawyClimbingSpeedWeek5.csv", header = F, as.is = T)

climbing.speed.wk1and5 <- merge(climbing.speed.wk1, climbing.speed.wk5, by.x = 1, by.y = 1)
climbing.speed.wk1and5 <- within(climbing.speed.wk1and5, decline <- V2.x - V2.y)
climbing.speed.wk1and5[, 1] <- paste("line_", climbing.speed.wk1and5[, 1], sep = "")
rownames(climbing.speed.wk1and5) <- climbing.speed.wk1and5[, 1]

endurance.wk1 <-  read.csv("gabrawyEnduranceWeek1.csv", header = F, as.is = T)
endurance.wk5 <-  read.csv("gabrawyEnduranceWeek5.csv", header = F, as.is = T)

endurance.wk1and5 <- merge(endurance.wk1, endurance.wk5, by.x = 1, by.y = 1)
endurance.wk1and5 <- within(endurance.wk1and5, decline <- V2.x - V2.y)
endurance.wk1and5[, 1] <- paste("line_", endurance.wk1and5[, 1], sep = "")
rownames(endurance.wk1and5) <- endurance.wk1and5[, 1]


# 4. lifespan

female.lifespan <- read.table("female.25c.mean.pheno", header = FALSE, as.is = TRUE)
female.lifespan[, 1] <- paste("line_", female.lifespan[, 1], sep = "")
male.lifespan <- read.table("male.25c.mean.pheno", header = FALSE, as.is = TRUE)
male.lifespan[, 1] <- paste("line_", male.lifespan[, 1], sep = "")
rownames(female.lifespan) <- female.lifespan[, 1]
rownames(male.lifespan) <- male.lifespan[, 1]

# merge and save file for plotting
# ============================================================

# female
female.pheno <- merge(merge(fecundity.wk1and3[, c(1, 5)], phototaxis.female.wk1and4[, c(1, 4)], by.x = 1, by.y = 1, all = TRUE), female.lifespan, by.x = 1, by.y = 1, all = TRUE)
colnames(female.pheno) <- c("line", "fecundity.decline", "phototaxis.decline", "lifespan")

# male
male.pheno <- merge(merge(merge(phototaxis.male.wk1and4[, c(1, 4)], climbing.speed.wk1and5[, c(1, 4)], by.x = 1, by.y = 1, all = TRUE), endurance.wk1and5[, c(1, 4)], by.x = 1, by.y = 1, all = TRUE), male.lifespan, by.x = 1, by.y = 1, all = TRUE)
colnames(male.pheno) <- c("line", "phototaxis.decline", "speed.decline", "endurance.decline", "lifespan")

save(female.pheno, male.pheno, file = "../figureData/pheno.merge.RData")

# make adjustment
# ============================================================

load("../figureData/adjustData.RData")
adjustPheno <- function(raw.pheno) {

# raw.pheno is a two column data.frame where the first is
# line id and the second is raw phenotype

  common.line <- intersect(rownames(wolba), raw.pheno[, 1])
  rownames(raw.pheno) <- raw.pheno[, 1]
  # get data
  pheno.data <- cbind(raw.pheno[common.line, 2], wolba[common.line, 1], inv[common.line, c("In_2L_t", "In_2R_NS", "In_3R_P", "In_3R_K", "In_3R_Mo")], pcs[common.line, ])
  colnames(pheno.data)[2] <- "wolba"
  fit.form <- "pheno.data[, 1] ~ 1"
  if (length(unique(pheno.data[, "wolba"])) > 1) {
    fit.form <- paste(fit.form, " + factor(wolba)")
  }
  if (length(unique(pheno.data[, "In_2L_t"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2L_t)")
  }
  if (length(unique(pheno.data[, "In_2R_NS"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_2R_NS)")
  }
  if (length(unique(pheno.data[, "In_3R_P"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_P)")
  }
  if (length(unique(pheno.data[, "In_3R_K"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_K)")
  }
  if (length(unique(pheno.data[, "In_3R_Mo"])) > 1) {
    fit.form <- paste(fit.form, " + factor(In_3R_Mo)")
  }
  fit.form <- paste(fit.form, " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11", sep = "")

  lm.fit <- lm(formula(fit.form), data = pheno.data)
  
  covar.coef <- summary(lm.fit)$coefficients
  
  pheno.adjust <- residuals(lm.fit) + covar.coef[1]
  
  return(cbind(names(pheno.adjust), raw.pheno[names(pheno.adjust), 2], pheno.adjust))
  
}

# output phenotype
# adjust phenotype leads to spurious h2, also do raw line means
female.fecundity.decline.pheno <- adjustPheno(female.pheno[, c(1, 2)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.fecundity.decline.pheno[, 1], female.fecundity.decline.pheno[, 1], female.fecundity.decline.pheno[, 2])), file = "female.fecundity.decline.raw.pheno", sep = " ", col.names = F, row.names = F, quote = F)
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.fecundity.decline.pheno[, 1], female.fecundity.decline.pheno[, 1], female.fecundity.decline.pheno[, 3])), file = "female.fecundity.decline.pheno", sep = " ", col.names = F, row.names = F, quote = F)

female.phototaxis.decline.pheno <- adjustPheno(female.pheno[, c(1, 3)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.phototaxis.decline.pheno[, 1], female.phototaxis.decline.pheno[, 1], female.phototaxis.decline.pheno[, 2])), file = "female.phototaxis.decline.raw.pheno", sep = " ", col.names = F, row.names = F, quote = F)
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.phototaxis.decline.pheno[, 1], female.phototaxis.decline.pheno[, 1], female.phototaxis.decline.pheno[, 3])), file = "female.phototaxis.decline.pheno", sep = " ", col.names = F, row.names = F, quote = F)

female.lifespan.pheno <- adjustPheno(female.pheno[, c(1, 4)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.lifespan.pheno[, 1], female.lifespan.pheno[, 1], female.lifespan.pheno[, 2])), file = "female.lifespan.raw.pheno", sep = " ", col.names = F, row.names = F, quote = F)
write.table(rbind(c("FID", "IID", "pheno"), cbind(female.lifespan.pheno[, 1], female.lifespan.pheno[, 1], female.lifespan.pheno[, 3])), file = "female.lifespan.pheno", sep = " ", col.names = F, row.names = F, quote = F)


male.phototaxis.decline.pheno <- adjustPheno(male.pheno[, c(1, 2)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(names(male.phototaxis.decline.pheno), names(male.phototaxis.decline.pheno), male.phototaxis.decline.pheno)), file = "male.phototaxis.decline.pheno", sep = " ", col.names = F, row.names = F, quote = F)

male.speed.decline.pheno <- adjustPheno(male.pheno[, c(1, 3)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(names(male.speed.decline.pheno), names(male.speed.decline.pheno), male.speed.decline.pheno)), file = "male.speed.decline.pheno", sep = " ", col.names = F, row.names = F, quote = F)

male.endurance.decline.pheno <- adjustPheno(male.pheno[, c(1, 4)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(names(male.endurance.decline.pheno), names(male.endurance.decline.pheno), male.endurance.decline.pheno)), file = "male.endurance.decline.pheno", sep = " ", col.names = F, row.names = F, quote = F)

male.lifespan.pheno <- adjustPheno(male.pheno[, c(1, 5)])
write.table(rbind(c("FID", "IID", "pheno"), cbind(names(male.lifespan.pheno), names(male.lifespan.pheno), male.lifespan.pheno)), file = "male.lifespan.pheno", sep = " ", col.names = F, row.names = F, quote = F)


