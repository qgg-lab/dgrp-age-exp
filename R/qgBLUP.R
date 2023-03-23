# ============================================
# = quantitative genetics of gene expression =
# ============================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("normalize/female.gene.exp.adjust.RData", "adjustData.RData", 8)
rdata.file <- args[1]
covar.file <- args[2]
n.cpu <- as.numeric(args[3])

library("lmerTest")
library("doMC")
registerDoMC(n.cpu)

# load data
# ============================================================

load(rdata.file)
unique.id <- factor(paste(line.id, temp, sep = "-"))
line.id.18c <- line.id[temp == "18C"]; line.id.25c <- line.id[temp == "25C"]
wolba.18c <- wolba[temp == "18C"]; wolba.25c <- wolba[temp == "25C"]; 

load(covar.file)
inv.2l.t <- factor(inv[as.character(line.id), "In_2L_t"]); inv.2l.t.18c <- inv.2l.t[temp == "18C"]; inv.2l.t.25c <- inv.2l.t[temp == "25C"];
inv.2r.ns <- factor(inv[as.character(line.id), "In_2R_NS"]); inv.2r.ns.18c <- inv.2r.ns[temp == "18C"]; inv.2r.ns.25c <- inv.2r.ns[temp == "25C"];
inv.3r.k <- factor(inv[as.character(line.id), "In_3R_K"]); inv.3r.k.18c <- inv.3r.k[temp == "18C"]; inv.3r.k.25c <- inv.3r.k[temp == "25C"];
inv.3r.p <- factor(inv[as.character(line.id), "In_3R_P"]); inv.3r.p.18c <- inv.3r.p[temp == "18C"]; inv.3r.p.25c <- inv.3r.p[temp == "25C"];
inv.3r.mo <- factor(inv[as.character(line.id), "In_3R_Mo"]); inv.3r.mo.18c <- inv.3r.mo[temp == "18C"]; inv.3r.mo.25c <- inv.3r.mo[temp == "25C"];
pcs <- as.matrix(pcs[as.character(line.id), ]); pcs.18c <- pcs[temp == "18C", ]; pcs.25c <- pcs[temp == "25C", ]

# perform QG analysis
# ============================================================

gene.name <- rownames(gene.exp.adj)

blup.18c <- matrix(NA, nrow = length(gene.name), ncol = length(levels(line.id)))
blup.25c <- matrix(NA, nrow = length(gene.name), ncol = length(levels(line.id)))
blup.18c.nocovar <- blup.18c
blup.25c.nocovar <- blup.25c

all.res <- foreach(i = 1:length(gene.name)) %dopar% {
  
  y <- unlist(gene.exp.adj[i, ])
  res <- rep(NA, 4*length(levels(line.id)))
  
  y.18c <- y[temp == "18C"]
  y.25c <- y[temp == "25C"]

  # fit models
  # ============================================================
  
  lmer.18c <- lmer(y.18c ~ wolba.18c + inv.2l.t.18c + inv.2r.ns.18c + inv.3r.k.18c + inv.3r.p.18c + inv.3r.mo.18c + pcs.18c + (1|line.id.18c))
  res[1:length(levels(line.id))] <- (ranef(lmer.18c)[[1]] + fixef(lmer.18c)[1])[levels(line.id), 1]
  
  lmer.25c <- lmer(y.25c ~ wolba.25c + inv.2l.t.25c + inv.2r.ns.25c + inv.3r.k.25c + inv.3r.p.25c + inv.3r.mo.25c + pcs.25c + (1|line.id.25c))
  res[(length(levels(line.id)) + 1):(2*length(levels(line.id)))] <- (ranef(lmer.25c)[[1]] + fixef(lmer.25c)[1])[levels(line.id), 1]
  
  lmer.18c.nocovar <- lmer(y.18c ~ (1|line.id.18c))
  res[(2*length(levels(line.id)) + 1):(3*length(levels(line.id)))] <- (ranef(lmer.18c.nocovar)[[1]] + fixef(lmer.18c.nocovar)[1])[levels(line.id), 1]
  
  lmer.25c.nocovar <- lmer(y.25c ~ (1|line.id.25c))
  res[(3*length(levels(line.id)) + 1):(4*length(levels(line.id)))] <- (ranef(lmer.25c.nocovar)[[1]] + fixef(lmer.25c.nocovar)[1])[levels(line.id), 1]
  
  cat(gene.name[i], "\n")
  res
  
}

for (i in 1:length(all.res)) {
  
  blup.18c[i, ] <- all.res[[i]][1:length(levels(line.id))]
  blup.25c[i, ] <- all.res[[i]][(length(levels(line.id)) + 1):(2*length(levels(line.id)))]
  blup.18c.nocovar[i, ] <- all.res[[i]][(2*length(levels(line.id)) + 1):(3*length(levels(line.id)))]
  blup.25c.nocovar[i, ] <- all.res[[i]][(3*length(levels(line.id)) + 1):(4*length(levels(line.id)))]
  
}

line.order <- levels(line.id)

save(blup.18c, blup.25c, blup.18c.nocovar, blup.25c.nocovar, line.order, gene.name, file = args[4])

# session info
# ============================================================

sessionInfo()
