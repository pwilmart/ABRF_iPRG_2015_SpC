
library(edgeR)
library(limma)
library(tidyverse)
library(psych)
library(RColorBrewer)

# try and start with the actual results file and wrangle data
# Summary table has lines before and after that are not strictly in table
# skip the leading lines when loading the table and read in just the table rows
full <- read_tsv("grouped_protein_summary_8.txt", skip = 4,
                 n_max = 2848, guess_max = 2848)

# check what are contaminants (we do not want any EXTRA accessions excluded)
check <- full %>% 
  select(Accession, Description, Filter) %>%
  filter(Filter == "contaminant")
check

# add the new columns to data frame and sort
full  <-  full %>%
  mutate(acc = sapply(strsplit(Accession, " "), `[`, 1)) %>%
  mutate(ave = ifelse(is.na(Filter), rowMeans(select(., starts_with("Corrected_"))), 0.0)) %>%
  mutate(missing = apply(select(., starts_with("Corrected_")) == 0, 1, sum)) %>%
  arrange(desc(ave))
head(full)
nrow(full)

ggplot(full, aes(x = ave, y = missing)) +
  geom_point() + geom_smooth(method = loess) +
  ggtitle("Missing versus Average SpC") + labs( x = "Ave SpC", y = "Missing")

ggplot(full, aes(x = ave, y = missing)) +
  geom_point() + geom_smooth(method = loess) +
  coord_cartesian(xlim = c(0, 10)) +
  ggtitle("Missing versus Average SpC") + labs( x = "Ave SpC", y = "Missing")

# get proteins (rows) where average SpC is greater than 2.5
# edgeR will want integer counts so round the corrected counts
counts <- full %>%
  filter(ave >= 2.5) %>%
  select(starts_with("Corrected")) %>%
  round(., 0)
accessions <- full %>%
  filter(ave >= 2.5) %>%
  select(acc)

# create a frame to hold results from edgeR
results <- accessions

S1 <- 1:3
S2 <- 4:6
S3 <- 7:9
S4 <- 10:12
pairs.panels(counts[S1], main = "Sample 1")
pairs.panels(log10(counts[S1]+1))
pairs.panels(counts[S2], main = "Sample 2")
pairs.panels(log10(counts[S2]+1))
pairs.panels(counts[S3], main = "Sample 3")
pairs.panels(log10(counts[S3]+1))
pairs.panels(counts[S4], main = "Sample 4")
pairs.panels(log10(counts[S4]+1))

# load the data into edgeR data structures
# group labels need to be factors
group = factor(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3)))
yglm <- DGEList(counts = counts, group = group, genes = accessions$acc)

# run the TMM normalization (and library size corrections)
yglm <- calcNormFactors(yglm)

# check clustering (6 different out of 1748 may not do much)
plotMDS(yglm)

# create the experimental design matrix
design <- model.matrix(~group)
rownames(design) <- colnames(yglm)
design

# extimate the dispersion parameters and check everything
yglm <- estimateDisp(yglm, design)
# yglm

# fit statistical models (design matrix already in y$design)
fit <- glmQLFit(yglm)

# this tests if any conditions caused differences
any <- glmQLFTest(fit, coef = 2:4)
topTags(any)
summary(decideTests(any))

group = c(rep("S1", 3), rep("S2", 3), rep("S3", 3), rep("S4", 3))
yc <- DGEList(counts = counts, group = group, genes = accessions$acc)
yc <- calcNormFactors(yc)
yc <- estimateDisp(yc)
plotBCV(yc)

# Compute the normalized counts (start with counts)
# sample loading adjusts each channel to the same average total
lib_facs <- mean(colSums(counts)) / colSums(counts)

# the TMM factors are library adjustment factors (so divide by them)
norm_facs <- lib_facs / yc$samples$norm.factors

# compute the normalized data as a new data frame
counts_norm <- sweep(counts, 2, norm_facs, FUN = "*")
colnames(counts_norm) <- paste(colnames(counts), "TMMNorm", sep = "_")

# add normalized counts to results
results <- cbind(results, counts_norm)

# look at count distributions across samples
boxplot(log10(counts_norm + 0.25), 
        col = c(rep("purple", 3), rep("red", 3), rep("green", 3), rep("blue", 3)),
        xlab = 'Samples', ylab = 'Normalized Counts', 
        main = 'TMM Normalized data', notch = TRUE)

# input: data frame, output: vector of CVs (%)
CV <- function(df) {
  ave <- rowMeans(df)
  sd <- apply(df, 1, sd)
  cv <- 100 * sd / ave
}

# put CVs in a data frame to simplify plots and summaries
cv_frame <- data.frame(Raw = CV(counts), Norm = CV(counts_norm))
medians <- apply(cv_frame, 2, FUN = median)
round(medians, 2)
boxplot(cv_frame, notch = TRUE, main = "SpC CV distributions", ylim = c(0, 100), ylab = "CV (%)")

# make average vectors (for results)
S1.ave <- rowMeans(counts_norm[S1])
S2.ave <- rowMeans(counts_norm[S2])

# exact test
et1_2 <- exactTest(yc, pair = c("S1", "S2"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et1_2)
tt$table

summary(decideTests(et1_2))
plotMD(et1_2, main = "S1 versus S2")
abline(h = c(-1, 1), col = "black")

# make average vectors (for results)
S3.ave <- rowMeans(counts_norm[S3])

# exact test
et1_3 <- exactTest(yc, pair = c("S1", "S3"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et1_3)
tt$table # this works at Github
summary(decideTests(et1_3))
plotMD(et1_3, main = "S1 versus S3")
abline(h = c(-1, 1), col = "black")

# make average vectors (for results)
S4.ave <- rowMeans(counts_norm[S4])

# exact test
et1_4 <- exactTest(yc, pair = c("S1", "S4"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et1_4)
tt$table
summary(decideTests(et1_4))
plotMD(et1_4, main = "S1 versus S4")
abline(h = c(-1, 1), col = "black")

# exact test
et2_3 <- exactTest(yc, pair = c("S2", "S3"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et2_3)
tt$table
summary(decideTests(et2_3))
plotMD(et2_3, main = "S2 versus S3")
abline(h = c(-1, 1), col = "black")

# exact test
et2_4 <- exactTest(yc, pair = c("S2", "S4"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et2_4)
tt$table
summary(decideTests(et2_4))
plotMD(et2_4, main = "S2 versus S4")
abline(h = c(-1, 1), col = "black")

# exact test
et3_4 <- exactTest(yc, pair = c("S3", "S4"))

# look at top DE table, see up and down, and a MA plot
tt <- topTags(et3_4)
tt$table
summary(decideTests(et3_4))
plotMD(et3_4, main = "S3 versus S4")
abline(h = c(-1, 1), col = "black")

# make a little frame for each pairwise test
df1_2 <- data.frame(S1.ave, S2.ave, topTags(et1_2, n = Inf, sort.by = "none"))
colnames(df1_2) <- paste(colnames(df1_2), "1vs2", sep = ".")

df1_3 <- data.frame(S1.ave, S3.ave, topTags(et1_3, n = Inf, sort.by = "none"))
colnames(df1_3) <- paste(colnames(df1_3), "1vs3", sep = ".")

df1_4 <- data.frame(S1.ave, S4.ave, topTags(et1_4, n = Inf, sort.by = "none"))
colnames(df1_4) <- paste(colnames(df1_4), "1vs4", sep = ".")

df2_3 <- data.frame(S2.ave, S3.ave, topTags(et2_3, n = Inf, sort.by = "none"))
colnames(df2_3) <- paste(colnames(df2_3), "2vs3", sep = ".")

df2_4 <- data.frame(S2.ave, S4.ave, topTags(et2_4, n = Inf, sort.by = "none"))
colnames(df2_4) <- paste(colnames(df2_4), "2vs4", sep = ".")

df3_4 <- data.frame(S3.ave, S4.ave, topTags(et3_4, n = Inf, sort.by = "none"))
colnames(df3_4) <- paste(colnames(df3_4), "3vs4", sep = ".")

# add to results
results <- cbind(results, df1_2, df1_3, df1_4, df2_3, df2_4, df3_4)

# drop the original count data from results (duplicated columns are issue with merging)
colnames(results)

# table "full" has more rows than "results", so we want a left join merge
full_with_stats <- left_join(full, results, by = "acc")
colnames(full_with_stats)

write.table(full_with_stats, "PAW_grouped_proteins_with_stats.txt", sep = "\t", row.names = FALSE, na = " ")
sessionInfo()


