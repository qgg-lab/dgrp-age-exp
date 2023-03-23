# =======================================================
# = make figure of quantiles of standardized expression =
# =======================================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.gene.exp.scaled.quantile.RData", "../figureData/male.gene.exp.scaled.quantile.RData", "../figure/figureScaledExpOutlier.pdf", "Myriad Pro")

library("RColorBrewer")

load(args[1])
female.baseline.quantile <- baseline.quantile
female.3wk.quantile <- old.quantile
female.age <- age
female.flow.cell <- flow.cell
female.sample.name <- sample.name
female.flow.cell[, 2] <- gsub(",.*", "", female.flow.cell[, 2])

load(args[2])
male.baseline.quantile <- baseline.quantile
male.3wk.quantile <- old.quantile
male.age <- age
male.flow.cell <- flow.cell
male.sample.name <- sample.name
male.flow.cell[, 2] <- gsub(",.*", "", male.flow.cell[, 2])

# group samples by flow cell
# ============================================================

female.baseline.quantile <- female.baseline.quantile[, order(female.flow.cell[colnames(female.baseline.quantile), 2])]
male.baseline.quantile <- male.baseline.quantile[, order(male.flow.cell[colnames(male.baseline.quantile), 2])]

female.3wk.quantile <- female.3wk.quantile[, order(female.flow.cell[gsub("^X", "", colnames(female.3wk.quantile)), 2])]
male.3wk.quantile <- male.3wk.quantile[, order(male.flow.cell[gsub("^X", "", colnames(male.3wk.quantile)), 2])]

# prepare file to plot
# ============================================================

plot.col = c(brewer.pal(n = 9, name = "Blues")[c(9, 7, 4)], "grey80", brewer.pal(n = 9, name = "Reds")[c(4, 7, 9)])

file.width = 120
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width * 1.6/25.4, family = args[4])
#cairo_pdf(file = "Figure_PreliminaryNormalizationOutlier.pdf", width = file.width/25.4, height = file.width * 1.6/25.4, family = "Myriad Pro")

par(las = 1, tcl = -0.2, mar = c(2, 2, 1, 0.9), ps = 7, lwd = 0.5, mfrow = c(4, 1))

plot(c(1, 405), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(female.baseline.quantile[i, ], col = plot.col[i])
}

bad.array <- which(female.baseline.quantile[1, ] < -5 | female.baseline.quantile[7, ] > 5)

if (length(bad.array) > 0) { 
  
  for (i in 1:length(bad.array)) { 
    
    k <- ifelse(abs(female.baseline.quantile[1, bad.array[i]]) > abs(female.baseline.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], female.baseline.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(names(bad.array[i]), split = "_"))[2:3]
    text(bad.array[i], female.baseline.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))
  
  }

}
  
legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Young females")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw flow cell

female.baseline.fc <- female.flow.cell[gsub("^X", "", colnames(female.baseline.quantile)), 2]
female.baseline.fc.unique <- unique(female.flow.cell[gsub("^X", "", colnames(female.baseline.quantile)) ,2])

for (i in 1:length(female.baseline.fc.unique)) {

  segments(min(which(female.baseline.fc == female.baseline.fc.unique[i])), -15, max(which(female.baseline.fc == female.baseline.fc.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(female.baseline.fc == female.baseline.fc.unique[i])), ifelse(i %% 2, -14, -12.5), female.baseline.fc.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  
}


box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 405), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(female.3wk.quantile[i, ], col = plot.col[i])
}

bad.array <- which(female.3wk.quantile[1, ] < -5 | female.3wk.quantile[7, ] > 5)

if (length(bad.array) > 0) {

  for (i in 1:length(bad.array)) {

    k <- ifelse(abs(female.3wk.quantile[1, bad.array[i]]) > abs(female.3wk.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], female.3wk.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(names(bad.array[i]), split = "_"))[2:3]
    text(bad.array[i], female.3wk.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))

  }

}

legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Aged females")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw flow cell

female.3wk.fc <- female.flow.cell[gsub("^X", "", colnames(female.3wk.quantile)), 2]
female.3wk.fc.unique <- unique(female.flow.cell[gsub("^X", "", colnames(female.3wk.quantile)), 2])

for (i in 1:length(female.3wk.fc.unique)) {

  segments(min(which(female.3wk.fc == female.3wk.fc.unique[i])), -15, max(which(female.3wk.fc == female.3wk.fc.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(female.3wk.fc == female.3wk.fc.unique[i])), ifelse(i %% 2, -14, -12.5), female.3wk.fc.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])

}


box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 405), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(male.baseline.quantile[i, ], col = plot.col[i])
}

bad.array <- which(male.baseline.quantile[1, ] < -6 | male.baseline.quantile[7, ] > 6)

if (length(bad.array) > 0) {

  for (i in 1:length(bad.array)) {

    k <- ifelse(abs(male.baseline.quantile[1, bad.array[i]]) > abs(male.baseline.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], male.baseline.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(names(bad.array[i]), split = "_"))[2:3]
    text(bad.array[i], male.baseline.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))

  }

}

legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Young males")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw flow cell

male.baseline.fc <- male.flow.cell[gsub("^X", "", colnames(male.baseline.quantile)), 2]
male.baseline.fc.unique <- unique(male.flow.cell[gsub("^X", "", colnames(male.baseline.quantile)), 2])

for (i in 1:length(male.baseline.fc.unique)) {

  segments(min(which(male.baseline.fc == male.baseline.fc.unique[i])), -15, max(which(male.baseline.fc == male.baseline.fc.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(male.baseline.fc == male.baseline.fc.unique[i])), ifelse(i %% 2, -14, -12.5), male.baseline.fc.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])

}


box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)





plot(c(1, 405), c(-15, 15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
for (i in 1:7) {
  lines(male.3wk.quantile[i, ], col = plot.col[i])
}

bad.array <- which(male.3wk.quantile[1, ] < -5 | male.3wk.quantile[7, ] > 5)

if (length(bad.array) > 0) {

  for (i in 1:length(bad.array)) {

    k <- ifelse(abs(male.3wk.quantile[1, bad.array[i]]) > abs(male.3wk.quantile[7, bad.array[i]]), 1, 7)
    points(bad.array[i], male.3wk.quantile[k, bad.array[i]], pch = 1, cex = 2, col = plot.col[k])
    this.bad.array <- unlist(strsplit(names(bad.array[i]), split = "_"))[2:3]
    text(bad.array[i], male.3wk.quantile[k, bad.array[i]], paste("R", this.bad.array[1], "-", this.bad.array[2], sep = ""), pos = ifelse(k == 1, 1, 3), cex = 6/par("ps")/par("cex"))

  }

}

legend("topleft", bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[1:4], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[1:4], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)
legend(-14, -3, bty = "n", lty = 1, lwd = 1, col = rev(plot.col)[5:7], legend = rev(c("q1", "q5", "q25", "q50", "q75", "q95", "q99"))[5:7], cex = 6/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5)

title(main = expression(paste("Aged males")), cex.main = 9/par("ps")/par("cex"))
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5)
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = "Scaled expression", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))

# draw flow cell

male.3wk.fc <- male.flow.cell[gsub("^X", "", colnames(male.3wk.quantile)), 2]
male.3wk.fc.unique <- unique(male.flow.cell[gsub("^X", "", colnames(male.3wk.quantile)), 2])

for (i in 1:length(male.3wk.fc.unique)) {

  segments(min(which(male.3wk.fc == male.3wk.fc.unique[i])), -15, max(which(male.3wk.fc == male.3wk.fc.unique[i])), -15, lwd = 2, col = brewer.pal(8, "Dark2")[i %% 8 + 1])
  text(mean(which(male.3wk.fc == male.3wk.fc.unique[i])), ifelse(i %% 2, -14, -12.5), male.3wk.fc.unique[i], col = brewer.pal(8, "Dark2")[i %% 8 + 1])

}



box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
