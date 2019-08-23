# ------------------------------------------
# Merritt Burch
# 2019-07-24
# mbb262@cornell.edu
#
# Script to merge Ky21 and M162w phenotypes
# with the "correct" scoring system so that
# There are no omitted taxa/RILs
# ------------------------------------------

# Set directory
setwd("C:/Users/merri/git_projects/ga1s/data")

# All RIL names for Z014 and Z015
taxa <- read.csv("Z014_and_Z015_names.csv", header = TRUE)

# RILs present in Ga1-s study with phenotypes
phenos <- read.csv("Z014_and_Z015_phenos.csv", header = TRUE)

# Merge the two
hold <- merge(x = taxa, y = phenos,
              by.x = "Taxa", by.y = "Genotype",
              all.x = TRUE)

# Export file
write.csv(hold, file = "all_Z014_and_Z015_with_phenos.csv",
          row.names = FALSE, quote = FALSE)

# format genotype tables
ky21 <- read.csv("NAM_genos_imputed_20090807_Ky21_genosAsChars.csv", header = TRUE)
m162w <- read.csv("NAM_genos_imputed_20090807_M162w_genosAsChars.csv", header = TRUE)

# Only take needed columns for now
hold_ky21 <- ky21[3:dim(ky21)[1], 2:length(ky21)]
hold_m162w <- m162w[3:dim(m162w)[1], 2:length(m162w)]

# Do some grep replacements of genotypes
hold_ky21 <- as.data.frame(lapply(hold_ky21, function(y) gsub("0.5", "BKy", y)))
hold_m162w <- as.data.frame(lapply(hold_m162w, function(y) gsub("0.5", "BM", y)))

hold_ky21 <- as.data.frame(lapply(hold_ky21, function(y) gsub("0", "B", y)))
hold_m162w <- as.data.frame(lapply(hold_m162w, function(y) gsub("0", "B", y)))

hold_ky21 <- as.data.frame(lapply(hold_ky21, function(y) gsub("1.5", "KyB", y)))
hold_m162w <- as.data.frame(lapply(hold_m162w, function(y) gsub("1.5", "MB", y)))

hold_ky21 <- as.data.frame(lapply(hold_ky21, function(y) gsub("2", "Ky", y)))
hold_m162w <- as.data.frame(lapply(hold_m162w, function(y) gsub("2", "M", y)))

hold_ky21 <- as.data.frame(lapply(hold_ky21, function(y) gsub("1", "Het", y)))
hold_m162w <- as.data.frame(lapply(hold_m162w, function(y) gsub("1", "Het", y)))


# Merge with main (ky21)
genoPheno <- cbind(ky21[3:dim(ky21)[1],1], hold_ky21)
colnames(genoPheno) <- colnames(ky21)
final_ky21 <- rbind(ky21[1:2,], genoPheno)
final_ky21[1,1] <- ""
final_ky21[2,1] <- ""
final_ky21[,1] <- as.numeric(final_ky21[,1])
write.csv(final_ky21, file ="NAM_genos_imputed_20090807_Ky21_asChars.csv", row.names = FALSE, quote = FALSE)


# Merge with main (m162w)
genoPheno <- cbind(m162w[3:dim(m162w)[1],1], hold_m162w)
colnames(genoPheno) <- colnames(m162w)
final_m162w <- rbind(m162w[1:2,], genoPheno)
final_m162w[1,1] <- ""
final_m162w[2,1] <- ""
final_m162w[,1] <- as.numeric(final_m162w[,1])
write.csv(final_m162w, file ="NAM_genos_imputed_20090807_M162w_asChars.csv", row.names = FALSE, quote = FALSE)
