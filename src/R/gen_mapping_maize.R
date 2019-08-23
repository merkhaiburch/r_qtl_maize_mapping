# -------------------------------
# title: R/qtl example in maize
# author: Merritt Burch
# date: 2019-08-23
# email: merrittbburch@gmail.com
# 
# Script contains an example of 
# how to map traits in NAM using
# r/qtl with the original 1106
# NAM markers
# -------------------------------
  
# Clear out the workspace and load the qtl package
# rm(list=ls())
library(qtl)

# Set directory
setwd("C:/Users/merri/git_projects/ga1s/data")

# Read in and specify the format of the data files
# Alleles not coded as numbers --> this led to issues within R/qtl
# instead convereted numeric genotypes into "B","Ky","Het", "BKy", "KyB"
# Where:
# B = homozygous B73
# Ky = homozygous Ky21
# Het = heterozygous
# BKy = more B73 than Ky21
# KyB = more Ky21 than B73
ga_ky21 <- read.cross("csv", ".", "NAM_genos_imputed_20090807_Ky21_asChars.csv",
                      na.strings =  c("-", "NA"),
                      genotypes = c("B","Ky","Het", "BKy", "KyB"))


# Summary and statistics of the data

# Set the cross type to RIL
class(ga_ky21)[1] <-"riself"

# Check
summary(ga_ky21)

# Run jittermap to move/jitter markers at the same posotion so
#     each marker is at a unique position
ga_ky21 <- jittermap(ga_ky21)

# Plot histogram of data, genetic map, and pattern of missing data
#     (shows three separate plots)
plot(ga_ky21)

# Plot only genetic map
plotMap(ga_ky21)


# --------------------------------
# Calculate genotype probabilites
# --------------------------------

# calculate genotype probabilities
#   Uses HMM for calculating genotype probs
ga_ky21 <- calc.genoprob(ga_ky21,
                         step=1, 
                         error.prob=0.05, 
                         map.function="kosambi",
                         stepwidth="fixed", 
                         off.end=0)

# Estimate recombination fractions
ga_ky21 <- est.rf(ga_ky21)


# --------------------
# Begin QTL anaylsis
# --------------------

# --------------------------
#   Single QTL Analysis
#   Single Marker Analysis
# --------------------------

# Single QTL scan/ Single Marker Analysis
# Using extended Haley-Knott regression due to its use of genotype information
#   robustness, and ability to handle slective genotyping (pg. 102)
# Set seed so results come out the same every time we run the analysis
set.seed(748215)
out.ehk_ky21 <- scanone(ga_ky21, method = "ehk")

# Obtain summaries
summary(out.ehk_ky21)
summary(out.ehk_ky21, threshold = 1.4) # set threshold to wherever you want, maybe run permutaiton test below


# ------------------------
#   Single QTL Analysis
# ------------------------

# Single QTL scan
# Same as above but testing the "HK" or Hayley-Knott method
# that gives slightly different results
out.hk <- scanone(ga_ky21, method = "hk")

# Obtain summaries
summary(out.hk)
summary(out.hk, threshold = 1.5)

# Plot this (there is a cleaner version of this at the bottom of the document)
plot(out.hk)


# ---------------------------------------------
# Begin permutation tests and select top QTLs
# ---------------------------------------------

# Permutation tests
operm_ky21 <- scanone(ga_ky21, method = "ehk", n.perm = 10000, verbose = FALSE)

# Plot all LOD scores generated from permutation tests
# LOD score from R/qtl based off histogram of this test
plot(operm_ky21)

# Obtain summaries
# We should set the LOD threshold at ~2.9
summary(operm_ky21)

# Shows QTLs that passed LOD threshold
summary(out.ehk_ky21, perms=operm_ky21, alpha=0.05, pvalues=TRUE)


# -----------------------------------
# Interval estimates of QTL location
# -----------------------------------

# Pick chromosome with prominent peak and do analysis
# Include markers in prominent peak
# Shows confidence intervals around peaks
lodint(out.ehk_ky21, chr=4, expandtomarkers=TRUE)

# ---------------------------
# More data summaries
# ---------------------------

# Set LOD threshold in summary output
# Set below reccommended threshold to study other QTL peaks in more detail
summary(out.hk, threshold=1.5)

# Show summaries of QTLs organized by chromosome and col
summary(out.ehk_ky21, threshold=1.5, format="tabByCol")


# --------------------------
#   Two way scan for QTLs
# --------------------------

# Pre-leg work
ga_ky21 <- calc.genoprob(ga_ky21, step = 1)

# 2-way scan
out2 <- scantwo(ga_ky21, method = "hk", verbose = FALSE)

# Plot these results (heatmap)
#plot(out2)

# Permutation test
operm2 <- scantwo(ga_ky21, method="hk", n.perm=5)

# Summary output
# Shows all potential significant interactions between QTL
summary(out2, perms=operm2, alpha=0.5, pvalues=TRUE)



# -----------------------
#    Find QTL Effects
# -----------------------

# Find marker closest to peak, (use summary(out.ehk_ky21) to get position estimates)
mar <- find.marker(ga_ky21, chr=1, pos=84.9)
mar2 <- find.marker(ga_ky21, chr=10, pos=22.9)

# Plot phenotypes versues marker genotypes
# Plots show marker effects for B73 and Ky21 at each of the 
# QTL peaks listed above
plotPXG(ga_ky21, marker=mar)


# Plot phenotype means against genotypes at one or two markers
# Kind of the same plot as above but shows the effect that 
# each marker has at each QTL peak between B73 and Ky21
effectplot(ga_ky21,mname1 = mar)

# Plot joint effects of the two loci with plotPXG
# Shows the effect of two markers at the same time
# Contrasts joint QTL effects, i.e. compares 
# QTL at 1S, 2S, 4S, and 8S when markers alternate
# for B73 and Ky21, only show 2 because they were the significant
# interactions in the analysis above
plotPXG(ga_ky21, marker=c(mar, mar2))
plotPXG(ga_ky21, marker=c(mar2, mar))

# Effect plot with two markers
# Show effect of two QTL peaks at one time in homozygous backgrounds
# parallel=additive effect of two markers at QTL peak
# Intersecting/crossed=potential interaction between QTL
effectplot(ga_ky21, mname1=mar, mname2=mar2)




# --------------------
#   Multi-QTL scan
# --------------------

# The following analysis narrows down the QTL search on the top QTLs found in the single marker scan. 
# Doing the multi-QTL scan is supposed to yield clearer results and the model is better fit for 
# multiple QTLs and their interactions (if any).

# Calculate genotype probabilities
ga_ky21 <- calc.genoprob(ga_ky21, step = 1)

# Pick most prominent QTLs
qtl <- makeqtl(ga_ky21, chr=c(1,10), pos=c(84.9,22.9), what="prob")

# Fit most prominent QTLs to a model
out.fq <- fitqtl(ga_ky21, qtl=qtl, method="hk")

# Obtain summary output
summary(out.fq)

# Get summary of fitted QTL model
summary(fitqtl(ga_ky21, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))

# Make models of most prominent QTLs
out.fqi <- fitqtl(ga_ky21, qtl=qtl, method="hk", formula=y~Q1*Q2)

# Get summary of models
summary(out.fqi)

# Get additive value
addint(ga_ky21, qtl=qtl, method="hk")

# Refine QTL model
rqtl <- refineqtl(ga_ky21, qtl=qtl, method="hk")
rqtl

# Summarize fitted QTL model
summary(out.fqr <- fitqtl(ga_ky21, qtl=rqtl, method="hk"))

# Plot the LOD profile
plotLodProfile(rqtl)

# Plot only most prominent QTL peak chromosomes
plot(out.hk, chr=c(1,2,4,8), col="red", add=TRUE)

# Show additivity effects
out.aq <- addqtl(ga_ky21, qtl=rqtl, method="hk")

# Calculate and show QTL penalities
print(pen <- calc.penalties(operm2))

# Do a steppwise analysis of QTLs
out.sq <- stepwiseqtl(ga_ky21, max.qtl=5, penalties=pen, method="hk", verbose=FALSE)
out.sq



# -------------------------------------
# Make nice plots of QTLs
# Maybe look at other plotting script as well
# -------------------------------------

# Packages
library(ggplot2)


# ------------------------------
#    Plot LOD scores more nicely
# -------------------------------

# Subset data if needed for only certain chromosomes
#LOD.input <- out.ehk_ky21[as.character(out.ehk_ky21$chr) %in% '4', ] # Subsets by chromosome

# Plot the thing (put hash in front of facet to remove chr number at top of plot)
# ggplot(data=out.ehk_ky21, aes(x=pos, y = lod)) + 
#   geom_line(size = 0.5, alpha = .8) + 
#   facet_wrap(~chr, nrow = 1, scales = "free_x") +
#   theme_classic() +
#   theme(strip.text = element_text(face = "bold", size = 15),
#         axis.text.y = element_text(size = 15, colour = "black"),
#         axis.text.x = element_text(angle=45,hjust=0.5, size = 12, colour = "black"),
#         strip.background = element_rect(colour = "white"), 
#         axis.line = element_line(color = "black", size = 1),
#         panel.spacing = unit(.55, "lines")) +
#   scale_y_continuous(limits = c(0, 7.5), expand = c(0, 0)) +
#   labs(x = "Chromosome 10 (Position)",
#        y = "LOD Score",
#        color = "",
#        linetype = "")



# Make function for plotting individual chromosomes
# Plot the thing (put hash in front of facet to remove chr number at top of plot)
plotbyChr <- function(pickChr, xAxis, yLimit){
  LOD.input <- out.ehk_ky21[as.character(out.ehk_ky21$chr) %in% pickChr, ] # Subsets by chromosome
  ggplot(data=LOD.input, aes(x=pos, y = lod)) + # Begins to plot 
    geom_line(size = 0.5, alpha = .8) + 
    #facet_wrap(~chr, nrow = 1, scales = "free_x") +
    theme_classic() +
    theme(strip.text = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(angle=0,hjust=0.5, size = 12, colour = "black"),
          strip.background = element_rect(colour = "white"), 
          axis.line = element_line(color = "black", size = 1),
          panel.spacing = unit(.55, "lines")) +
    scale_y_continuous(limits = c(0, yLimit), expand = c(0, 0)) +
    labs(x = xAxis,
         y = "LOD Score",
         color = "",
         linetype = "")
}

# Plot chr1
plotbyChr(pickChr = '1', xAxis = "Chromosome 1 (Position)", yLimit = 2.3)


# ---------------------------
#   Plot additivity effects
# ---------------------------

tmp <- sim.geno(ga_ky21, map.function = "haldane")
effectscans <- effectscan(tmp, chr=1:10, draw = F)

# Whole genome
ggplot(data=effectscans, aes(x=pos, y = a)) + 
  geom_line(size = 0.8, alpha = .8) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  facet_wrap(~chr, nrow = 1, scales = "free_x") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(angle=45,hjust=1, size = 12, colour = "black"),
        strip.background = element_rect(colour = "white"), 
        axis.line = element_line(color = "black", size = 1),
        panel.spacing = unit(.55, "lines")) +
  scale_y_continuous(limits = c(-.35, .35), expand = c(0, 0)) +
  labs(x = "Chromosome Position",
       y = "Additive Effect",
       color = "",
       linetype = "")

# Function for specific chromosomes
plotAddbyChr <- function(pickChr, xAxis){
  add.input <- effectscans[as.character(effectscans$chr) %in% pickChr, ] # Subsets by chromosome
  ggplot(data=add.input, aes(x=pos, y = a)) + 
    geom_line(size = 0.5, alpha = .8) + 
    geom_hline(yintercept = 0, linetype = "solid") +
    #facet_wrap(~chr, nrow = 1, scales = "free_x") +
    theme_classic() +
    theme(strip.text = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(angle=45,hjust=1, size = 12, colour = "black"),
          strip.background = element_rect(colour = "white"), 
          axis.line = element_line(color = "black", size = 1),
          panel.spacing = unit(.55, "lines")) +
    scale_y_continuous(limits = c(-.35, .35), expand = c(0, 0)) +
    labs(x = "Chromosome Position",
         y = "Additive Effect",
         color = "",
         linetype = "")
  
}

# Plot 1S
plotAddbyChr(pickChr = '1', xAxis = "Chromosome 1 Position")

#Save r workspace so I don't have to save everything again
save.image(file = "b73xKy21_Burch_analyzed.RData")


