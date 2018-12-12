# Read search counts
source("data.R")

# Combine search counts
searches <- rbind(s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12,
                  s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24,
                  s25, s26, s27, s28, s29)
searches <- as.data.frame(searches)

# Name columns/variables
colnames(searches) <- c("pubmed", "proquest", "ebscohost", "wos", "ovid")

# Create column for search sets
searches$searchset <- rownames(searches)

# Cleanup
rm(s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12, s13, s14, s15,
   s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29)

# Calculate modified z-score *m_i*
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

# Name columns/variables
colnames(modz) <- c("pubmed", "proquest", "ebscohost", "wos", "ovid")

# Replace NaNs and Infs with 0 to show perfect correlation since they are
# perfectly or nearly perfectly correlated
modz[c(16, 17, 19, 20), ] <- c(0, 0, 0, 0)

# Create column for search sets
modz$searchset <- rownames(searches)

# Plots
# gridExtra is for a multiplot
library("ggplot2")
library("svglite")
library("gridExtra")

##### FIGURE 1: Raw Scores #####
# `options` will avoid scientific notation
options(scipen = 10000)
raw0 <- ggplot(searches, aes(searchset, proquest)) +
  geom_col() +
  labs(x = "Search Set", y = "PubMed") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black")) +
theme(axis.text.y = element_text(size = 10))

raw1 <- ggplot(searches, aes(searchset, proquest)) +
  geom_col() +
  labs(x = "Search Set", y = "ProQuest") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black")) +
theme(axis.text.y = element_text(size = 10))

raw2 <- ggplot(searches, aes(searchset, ebscohost)) +
  geom_col() +
  labs(x = "Search Set", y = "EBSCOhost") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black")) +
theme(axis.text.y = element_text(size = 10))

raw3 <- ggplot(searches, aes(searchset, wos)) +
  geom_col() +
  labs(x = "Search Set", y = "Web of Science") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black")) +
theme(axis.text.y = element_text(size = 10))

raw4 <- ggplot(searches, aes(searchset, ovid)) +
  geom_col() +
  labs(x = "Search Set", y = "Ovid") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black")) +
theme(axis.text.y = element_text(size = 10))

fig1 <- grid.arrange(raw0, raw1, raw2, raw3, raw4, ncol = 3, nrow = 2)
ggsave("plots/figure-1-raw.svg", plot = fig1, height = 9,
       width = 12, dpi = 300)

##### FIGURE 2: Raw Score Diffs from PubMed MEDLINE #####
raw1.diff <- ggplot(searches, aes(searchset, proquest - pubmed)) +
  geom_col() +
  labs(x = "Search Set", y = "ProQuest") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

raw2.diff <- ggplot(searches, aes(searchset, ebscohost - pubmed)) +
  geom_col() +
  labs(x = "Search Set", y = "EBSCOhost") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

raw3.diff <- ggplot(searches, aes(searchset, wos - pubmed)) +
  geom_col() +
  labs(x = "Search Set", y = "Web of Science") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

raw4.diff <- ggplot(searches, aes(searchset, ovid - pubmed)) +
  geom_col() +
  labs(x = "Search Set", y = "Ovid") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

fig2 <- grid.arrange(raw1.diff, raw2.diff, raw3.diff, raw4.diff,
                     ncol = 2, nrow = 2)
ggsave("plots/figure-2-raw-diffs.svg", plot = fig2,
       height = 9, width = 12, dpi = 300)

##### FIGURE 3: Modified z scores #####
# Remove extreme outliers
modz1 <- subset(modz, proquest <= 3.5 & proquest >= -3.5)
mz1 <- ggplot(modz1, aes(searchset, proquest)) + geom_col() +
  labs(x = "Search Set", y = "ProQuest") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

modz2 <- subset(modz, ebscohost <= 3.5 & ebscohost >= -3.5)
mz2 <- ggplot(modz2, aes(searchset, ebscohost)) + geom_col() +
  labs(x = "Search Set", y = "EBSCOhost") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

modz3 <- subset(modz, wos <= 3.5 & wos >= -3.5)
mz3 <- ggplot(modz3, aes(searchset, wos)) + geom_col() +
  labs(x = "Search Set", y = "Web of Science") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

modz4 <- subset(modz, ovid <= 3.5 & ovid >= -3.5)
mz4 <- ggplot(modz4, aes(searchset, ovid)) + geom_col() +
  labs(x = "Search Set", y = "Ovid") + theme_bw() + coord_flip() +
       theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   colour = "black"))

fig3 <- grid.arrange(mz1, mz2, mz3, mz4, ncol = 2, nrow = 2)
ggsave("plots/figure-3-modz-scores.svg", plot = fig3,
       height = 9, width = 12, dpi = 300)

##### FIGURE 4: ProQuest Extreme Outliers #####
# Plot ProQuest extreme outliers
modz1prouqestoutliers <- subset(modz, proquest >= 3.5 | proquest <= -3.5)
modzpqout <- ggplot(modz1prouqestoutliers, aes(searchset, proquest)) +
  geom_col() + labs(x = "", y = "ProQuest") + theme_bw() + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"))
modzpqout
ggsave("plots/figure-4-proquest-extreme-outliers.svg", height = 9, width = 12,
       dpi = 300)

##### FIGURE 5: Web of Science Extreme Outliers #####

modz3wosoutliers <- subset(modz, wos >= 3.5 | wos <= -3.5)
modzwosout <- ggplot(modz3wosoutliers, aes(searchset, wos)) +
  geom_col() +
  labs(x = "", y = "Web of Science") + theme_bw() + coord_flip() +
       theme(axis.text.x =
             element_text(angle = 45, hjust = 1, colour = "black"))
modzwosout
ggsave("plots/figure-5-wos-extreme-outliers.svg", height = 9, width = 12,
       dpi = 300)

# Save data as tables tables and incorporate into LibreOffice
library("xtable")
modztable <- as.data.frame(round(modz[, c(2, 3, 4, 5)], 3))

searchdiffs <- as.data.frame(cbind(searches$proquest - searches$pubmed,
                               searches$ebscohost - searches$pubmed,
                               searches$wos - searches$pubmed,
                               searches$ovid - searches$pubmed))

colnames(searchdiffs) <- c("proquest", "ebscohost", "wos", "ovid")

searchdiffs$searchset <- rownames(searches)
 
# In Bash, convert to odt with pandoc; e.g.:
# pandoc -o table1.odt table1.tex
print(xtable(searches),
      file = "tables/table1.tex")
print(xtable(modztable),
      file = "tables/table2.tex")
print(xtable(searchdiffs),
      file = "tables/table3.tex")
