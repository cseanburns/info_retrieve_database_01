# Read search counts
source("data.R")

# combine search counts
searches <- rbind(s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12,
                  s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24,
                  s25, s26, s27, s28, s29)

searches <- as.data.frame(searches)

# Name columns/variables
colnames(searches) <- c("pubmed", "proquest", "ebscohost", "wos", "ovid")

# Cleanup
rm(s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12, s13, s14, s15,
   s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29)

# calculate modified z-score *m_i*
# function based on pubmed count and not median as center
modz <- function(searches, x) {
  s     <- as.numeric(searches[x, 1:5])
  sabs  <- abs(s - s[1])
  smad  <- median(sabs)
  smodz <- (0.6745 * (s - s[1])) / smad
  return(smodz)
}

# save output of for loop into a dataframe
# bless this package!
library("magicfor")
magic_for(print, silent = TRUE)

for (i in 1:29) {
  print(modz(searches, i))
}

modz <- magic_result_as_vector()
modz <- as.data.frame(matrix(modz, ncol = 5, byrow = T))

colnames(modz) <- c("pubmed", "proquest", "ebscohost", "wos", "ovid")

# Replace NaNs and Infs with 0 to show perfect correlation since they are
# perfectly or nearly perfectly correlated
modz[c(16, 17, 19, 20), ] <- c(0, 0, 0, 0)

# Create 
# gridExtra is for a multiplot
library("ggplot2")
library("gridExtra")

# labs(x = "") for the first two since x
mz1 <- ggplot(modz, aes(x = pubmed, y = proquest)) +
  geom_point() + labs(x = "", y = "ProQuest") + theme_bw()

mz2 <- ggplot(modz, aes(x = pubmed, y = ebscohost)) +
  geom_point() + labs(x = "", y = "EBSCOhost") + theme_bw()

mz3 <- ggplot(modz, aes(x = pubmed, y = wos)) +
  geom_point() + labs(x = "PubMed", y = "Web of Science") + theme_bw()

mz4 <- ggplot(modz, aes(x = pubmed, y = ovid)) +
  geom_point() + labs(x = "PubMed", y = "Ovid") + theme_bw()

png("plots/multi-point-plot.png", height = 9, width = 9,
    units = "in", res = 300)
grid.arrange(mz1, mz2, mz3, mz4, ncol = 2, nrow = 2)
dev.off()

# View of the scores that are greater or lesser than 1 median
# standard deviations from the PubMed scores:
# then, total how many scores are plus or minus 1 modified deviations

pqoutliers <- function(modz, proquest) {
  pqout1 <- cbind(as.matrix(modz$proquest),
                (modz$proquest >= 1 | modz$proquest <= -1))
  pqsum1 <- sum(pqout1[, 2])
  pqout2 <- cbind(as.matrix(modz$proquest),
                (modz$proquest >= 2 | modz$proquest <= -2))
  pqsum2 <- sum(pqout2[, 2])
  pqsum1 <- pqsum1 - pqsum2
  cat("result: more than one, more than two\n")
  return(c(pqsum1, pqsum2))
}

pqoutliers(modz, proquest)

eboutliers <- function(modz, ebscohost) {
  ebout1 <- cbind(as.matrix(modz$ebscohost),
                (modz$ebscohost >= 1 | modz$ebscohost <= -1))
  ebsum1 <- sum(ebout1[, 2])
  ebout2 <- cbind(as.matrix(modz$ebscohost),
                (modz$ebscohost >= 2 | modz$ebscohost <= -2))
  ebsum2 <- sum(ebout2[, 2])
  ebsum1 <- ebsum1 - ebsum2
  cat("result: more than one, more than two\n")
  return(c(ebsum1, ebsum2))
}

eboutliers(modz, ebscohost)

wosoutliers <- function(modz, wos) {
  wosout1 <- cbind(as.matrix(modz$wos),
                (modz$wos >= 1 | modz$wos <= -1))
  wossum1 <- sum(wosout1[, 2])
  wosout2 <- cbind(as.matrix(modz$wos),
                (modz$wos >= 2 | modz$wos <= -2))
  wossum2 <- sum(wosout2[, 2])
  wossum1 <- wossum1 - wossum2
  cat("result: more than one, more than two\n")
  return(c(wossum1, wossum2))
}

wosoutliers(modz, wos)

ovidoutliers <- function(modz, ovid) {
  ovidout1 <- cbind(as.matrix(modz$ovid),
                (modz$ovid >= 1 | modz$ovid <= -1))
  ovidsum1 <- sum(ovidout1[, 2])
  ovidout2 <- cbind(as.matrix(modz$ovid),
                (modz$ovid >= 2 | modz$ovid <= -2))
  ovidsum2 <- sum(ovidout2[, 2])
  ovidsum1 <- ovidsum1 - ovidsum2
  cat("result: more than one, more than two\n")
  return(c(ovidsum1, ovidsum2))
}

ovidoutliers(modz, ovid)

library("xtable")
modztable <- as.data.frame(round(modz[, -1], 3))

# convert to odt with pandoc
print(xtable(modztable), file = "docs/table1.tex")
