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

# usage: modz (row, startcol, endcol)
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

# Create a correlation matrix plot
# Remove pubmed since those are zeroed out
library("corrplot")
modcor <- modz[, -1]
modcor <- cor(modcor, method = "spearman")
png("corrplot.png", height = 9, width = 9, units = "in",
    res = 300)
corrplot.mixed(modcor, upper = "ellipse", number.cex = 1.0,
               addgrid.col = "grey", lower.col = "black",
               upper.col = grey.colors(100))
dev.off()

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

png("multi-point-plot.png", height = 9, width = 9, units = "in",
       res = 300)
grid.arrange(mz1, mz2, mz3, mz4, ncol = 2, nrow = 2)
dev.off()
