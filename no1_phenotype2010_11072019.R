## Phenotype = frm Dr. Meechai at seedling stage in 2010 both DH103 and DH212
# test
library(qtl)
library(data.table)
library(GenomicRanges)
library(ggbio)
library(ggplot2)
library(IRanges)
library(ggrepel)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)   # For using str_detect

Pheno <- read.delim("phenotype2010_10072019.csv", sep = ",",stringsAsFactors=F)
Pheno

# SNP of DH212
SNPs <- read.delim("QTLinputData_ABformat_DH212",colClasses="character")
SNPs[,1:2] <- sapply(SNPs[,1:2],as.numeric)
SNPs$cM <- SNPs$POS/1e06
colnames(SNPs)
SNPs <- SNPs[,c(1,2,120,15:119)]  # CHROM, POS, cM, CSSLs
dim(SNPs)
colnames(SNPs)
SNPs[1:5,1:8]
Keep <- apply(SNPs[,4:length(SNPs)],1,function(x) length(unique(na.omit(x))) == 2) # count unique values across rows(omit NA)
SNPsAB <- SNPs[Keep,]
SNPsAB[,4:108] <- lapply(SNPsAB[,4:108], function(x) gsub("AA","A",x))  # AA to A (kdml105)
SNPsAB[,4:108] <- lapply(SNPsAB[,4:108], function(x) gsub("BB","B",x))  # BB to B (DH)
SNPsAB[1:5,1:10]
SNPsAB$Name <- paste(SNPsAB$CHROM,SNPsAB$POS,sep="_") # add more column as CHROM_POS
SNPsAB <- SNPsAB[order(SNPsAB$CHROM,SNPsAB$POS),]   # order POS
SNPsAB[1:5,1:10]
dim(SNPsAB)
colnames(SNPsAB)
tSNP <- t(SNPsAB[,c(109,1,3,4:108)]) # col = POS, row = CSSL line
tSNP[,1:10]
QTLSample <- gsub("_GBS","",rownames(tSNP)[4:108])   # remove '_GBS' after CSSL name

UsePheno <- Pheno[match(QTLSample,Pheno$Sample),]  # Add cssl name to Phenotype file
rownames(UsePheno) <- NULL
colnames(UsePheno)
QTLPheno <- rbind(colnames(UsePheno)[2:26],rep("",26),rep("",26),UsePheno[,2:26]) # Add a row of SIS3-12

QTL_input <- cbind(QTLPheno,tSNP) # combine SIS score (phenotype) + SNPs for each line
QTL_input[, 1:30] # colnames = row no. of SNP previously (is nothing)
dim(QTL_input)

# delete rows that contain NA in Phenotypes

NA.delete <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

ans <- NA.delete(QTL_input, "TDW.C") # delete row that have NA at the first phenotype column
dim(ans)

rownames(ans)
QTL_input <- ans[c(1:3,12:4, 25:18,13),]

rownames(QTL_input)


rownames(QTL_input)
write.table(QTL_input,file="QTL_input_phenotype2010_10072019",row.names=F,col.names=F,sep=",",quote=F,na="NA")

# Read file "cross' with backcross type (BC5F3)
Salt <- read.cross(format="csv",file="QTL_input_phenotype2010_10072019",crosstype="bcsft",BC.gen=5,F.gen=3)

summary(Salt)
plotMap(Salt)  # Location of SNP on chromosome
phenames(Salt) # Phenotypic datasets (SISs)

geno.image(Salt)  # Plot grid of genotype data (A = KDML, B = DH103)
# default: missing = white, 
#AA = red, AB = blue, BB = green, 
#not BB = purple, not AA = orange
geno.image(Salt, col = c("white", "grey", "red", "#228B22"), reorder = FALSE)

# Or CMplot
library(CMplot)
dim(SNPsAB)
colnames(SNPsAB)
test <- SNPsAB[,c(109,1:2)]
table(test$CHROM)
head(test)
rownames(test) <- NULL
colnames(test) <- c("SNP","Chromosome", "Position")

SNP <- rep("X", 12)
Chromosome <- c(1:12)
p <- c(43270923, 35937250, 36413819, 35502694, 29958434,
       31248787, 29697621, 28443022, 23012720, 23207287,
       29021106, 27531856)

chrlength <- cbind(SNP, Chromosome, p)
colnames(chrlength) <- c("SNP","Chromosome", "Position")
all <- rbind.data.frame(test, chrlength)
head(all)
dim(all)
all[3180:3183,]

CMplot(all,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE, verbose=TRUE)
# Perform some intermefiate calculation before doing QTL analyses
# Simulate sequences of genotype from their joint distribution, 'given the observed marker data'
Salt <- sim.geno(Salt) 
# calculate conditional genotype probabilities 'given the multipoint marker data'
Salt <- calc.genoprob(Salt)

# Manhattan plot https://www.researchgate.net/figure/Manhattan-plot-of-QTL-LOD-scores-for-starch-viscosity-traits-The-Manhattan-plots-display_fig3_270289496 
# Manhattan plot https://www.researchgate.net/figure/Manhattan-plot-of-QTL-LOD-scores-for-starch-viscosity-traits-The-Manhattan-plots-display_fig3_270289496 
https://github.com/YinLiLin/R-CMplot#a-all-traits-in-a-axes 
library("CMplot")

out1.em.H <- scanone(Salt, method = "hk", pheno.col=1) # phenotype column1 (SIS3)

operm.em.H1 <- scanone(Salt, method = "hk",pheno.col = 1, n.perm=1000, n.clust=3)
lod_threshold <- summary(operm.em.H1, alpha = 0.10)
summary(out1.em.H, perms=operm.em.H1,alpha = 0.10, pvalues=TRUE)

labels_df <- as.data.frame(summary(out1.em.H, perms = operm.em.H1, alpha = 0.10, 
                                   pvalues = TRUE))
summary(operm.em.H1, alpha=c(0.05, 0.10))
labels_df

phenames(Salt) # Phenotypic datasets (SISs)

a <- out1.em.H
a$SNP <- rownames(a)
colnames(a) <- c("Chromosome", "pos", "lod", "SNP") 
a <- a[,c(4, 1:3)]
rownames(a) <- NULL
dim(a)
head(a)
max(a$lod)

summary(operm.em.H1, alpha=c(0.05, 0.10))

a %>% filter(lod > 4.09)

# this one
CMplot(a, plot.type = "m", band = 0.5, LOG10 = FALSE,
       col=c("#58574b","#985396"), ylab = "LOD score",
       ylim = c(0,10),
       cex = 0.7,
       threshold = 3.93, threshold.lty = 2, 
       threshold.lwd = 1, threshold.col = "black",
       amplify = TRUE, chr.den.col = c("darkgreen", "yellow", "red"),
       signal.col = "green", signal.cex = c(1.5,1.5))

