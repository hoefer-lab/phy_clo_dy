rm(list=ls())

######################## ######################## ######################## ######################## ########################
## path to data; Must contain subfolders vcfs (containg the files primary.vcf, relapse.vcf, metastasis.vcf) and cnvs (containg the files tumor_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt,
## relapse_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt, metastasis_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt)
data.location <- ""
## path to store output
output.path <- ""

######################## ######################## ######################## ######################## ########################
## load libraries and source function files

source("~/multi_sample/Phylogenetic_inference.R")
library(ggplot2)
library(Homo.sapiens)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(igraph)

rowMaxs <- function(matrix){
  apply(matrix, 1, max)
}

######################## ######################## ######################## ######################## ########################
## settings
normal.copy.number.auto <- 2
normal.copy.number.y <- 1
normal.b.allele.norm <- 1
sex <- "male"
## Extra loci we want to keep track of
extra.muts <- c("2p", "3p", "3q", "5p", "TERT", "7p", "7q", "10p", "10q", "PTEN", "11p", "16q", "17p", "7.1", "7.2", "7.3")
######################## ######################## ######################## ######################## ########################
## First we analyze the mutations present in the samples to infer a base tree

## filtered SNV files 
primary <- read.delim(paste0(data.location, "/vcfs/primary.vcf"), sep="\t", header = T, stringsAsFactors = F)
relapse <- read.delim(paste0(data.location,"/vcfs/relapse.vcf"), sep="\t", header = T, stringsAsFactors = F)
metastasis <- read.delim(paste0(data.location,"/vcfs/metastasis.vcf"), sep="\t", header = T, stringsAsFactors = F)

## prepare matrix to store read counts
readcounts.ref <- matrix(0, nrow=length(unique(c(paste(primary$GENE, primary$POS, sep="."), paste(relapse$GENE, relapse$POS, sep="."),
                                paste(metastasis$GENE, metastasis$POS, sep=".")))), ncol=3,
                         dimnames = list(unique(c(paste(primary$GENE, primary$POS, sep="."), paste(relapse$GENE, relapse$POS, sep="."),
                                                paste(metastasis$GENE, metastasis$POS, sep="."))),
                                         c("Primary", "Relapse", "Metastasis")))

rownames(readcounts.ref) <- unique(c(apply(primary,1,function(x){
  y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
  paste(y, as.character(x["POS"]), sep=".")
}),
apply(relapse,1,function(x){
  y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
  paste(y, as.character(x["POS"]), sep=".")
}),
apply(metastasis,1,function(x){
  y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
  paste(y,  as.character(x["POS"]), sep=".")
})))
readcounts.mut <- readcounts.ref
  
  ## extract read counts
primary.ref <-  sapply(primary$INFO, function(x){
 y <- x
  x <- strsplit(x, split=";")[[1]][6]
  x <- strsplit(x, split="=")[[1]][2]
  x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
  if(is.na(x)){
    x <- strsplit(y, split=";")[[1]][5]
    x <- strsplit(x, split="=")[[1]][2]
    x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
  }
  x
})
primary.mut <-  sapply(primary$INFO, function(x){
 y <- x
 x <- strsplit(x, split=";")[[1]][6]
 x <- strsplit(x, split="=")[[1]][2]
 x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 if(is.na(x)){
   x <- strsplit(y, split=";")[[1]][5]
   x <- strsplit(x, split="=")[[1]][2]
   x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 }
 x
})
relapse.ref <-  sapply(relapse$INFO, function(x){
 y <- x
 x <- strsplit(x, split=";")[[1]][6]
 x <- strsplit(x, split="=")[[1]][2]
 x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
 if(is.na(x)){
   x <- strsplit(y, split=";")[[1]][5]
   x <- strsplit(x, split="=")[[1]][2]
   x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
 }
 x
})
relapse.mut <-  sapply(relapse$INFO, function(x){
 y <- x
 x <- strsplit(x, split=";")[[1]][6]
 x <- strsplit(x, split="=")[[1]][2]
 x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 if(is.na(x)){
   x <- strsplit(y, split=";")[[1]][5]
   x <- strsplit(x, split="=")[[1]][2]
   x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 }
 x
})
metastasis.ref <-  sapply(metastasis$INFO, function(x){
 y <- x
 x <- strsplit(x, split=";")[[1]][6]
 x <- strsplit(x, split="=")[[1]][2]
 x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
 if(is.na(x)){
   x <- strsplit(y, split=";")[[1]][5]
   x <- strsplit(x, split="=")[[1]][2]
   x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(1,2)]))
 }
 x
})
metastasis.mut <-  sapply(metastasis$INFO, function(x){
 y <- x
 x <- strsplit(x, split=";")[[1]][6]
 x <- strsplit(x, split="=")[[1]][2]
 x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 if(is.na(x)){
   x <- strsplit(y, split=";")[[1]][5]
   x <- strsplit(x, split="=")[[1]][2]
   x <- sum(as.numeric(strsplit(x, split=",")[[1]][c(3,4)]))
 }
 x
})
 
primary.names <- apply(primary,1,function(x){
 y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
 paste(y, as.character(x["POS"]), sep=".")
})
relapse.names <- apply(relapse,1,function(x){
 y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
 paste(y, as.character(x["POS"]), sep=".")
})
metastasis.names <- apply(metastasis,1,function(x){
 y <- strsplit(as.character(x["GENE"]), split="\\(")[[1]][1]
 paste(y, as.character(x["POS"]), sep=".")
})

readcounts.ref[primary.names,1] <- primary.ref
readcounts.ref[relapse.names,2] <- relapse.ref
readcounts.ref[metastasis.names,3] <- metastasis.ref

readcounts.mut[primary.names,1] <- primary.mut
readcounts.mut[relapse.names,2] <- relapse.mut
readcounts.mut[metastasis.names,3] <- metastasis.mut

## add extra loci
mean.cov <- round(apply(readcounts.ref + readcounts.mut, 2, function(x){mean(x[x!=0])}))
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)

readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)
readcounts.ref <- rbind(readcounts.ref, mean.cov)

rownames(readcounts.ref)[(nrow(readcounts.ref)-15):nrow(readcounts.ref)] <- extra.muts

readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)

readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)
readcounts.mut <- rbind(readcounts.mut, 0)

rownames(readcounts.mut)[(nrow(readcounts.mut)-15):nrow(readcounts.mut)] <- extra.muts

extra.muts <- c("2p", "3p", "3q", "5p", "TERT", "7p", "7q", "10p", "10q", "PTEN", "11p", "16q", "17p", "7.1", "7.2", "7.3")

cnv.primary <- read.delim(paste0(data.location, "/cnvs/tumor_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt"), stringsAsFactors = F)
cnv.relapse <- read.delim(paste0(data.location, "/cnvs/relapse_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt"), stringsAsFactors = F)
cnv.metastasis <- read.delim(paste0(data.location, "/cnvs/metastasis_B060_Stuttgart_case_merged.mdup.bam_minipileup.pileup_ratio.txt"), stringsAsFactors = F)

cnv.files <- list(cnv.primary, cnv.relapse, cnv.metastasis)

## add as extra loci all loci where the coverage ratio changes, such that we can later plot a copy number profile
tmp <- cnv.primary[!duplicated(cnv.primary[,c(1,4)]),]
for(i in 1:nrow(tmp)){
 readcounts.mut <- rbind(readcounts.mut, 0)
 readcounts.ref <- rbind(readcounts.ref, mean.cov)
}
rownames(readcounts.mut)[(nrow(readcounts.mut)-nrow(tmp) +1):nrow(readcounts.mut)] <- paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep=".")
rownames(readcounts.ref) <- rownames(readcounts.mut)
extra.muts <- c(extra.muts, paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."))
tmp <- tmp[,c(1,2)]
dimnames(tmp) <- list(paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."), c("X.CHROM", "POS"))
locus.add <- tmp

tmp <- cnv.relapse[!duplicated(cnv.relapse[,c(1,4)]),]
for(i in 1:nrow(tmp)){
 readcounts.mut <- rbind(readcounts.mut, 0)
 readcounts.ref <- rbind(readcounts.ref, mean.cov)
}
rownames(readcounts.mut)[(nrow(readcounts.mut)-nrow(tmp) +1):nrow(readcounts.mut)] <- paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep=".")
rownames(readcounts.ref) <- rownames(readcounts.mut)
extra.muts <- c(extra.muts, paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."))
tmp <- tmp[,c(1,2)]
dimnames(tmp) <- list(paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."), c("X.CHROM", "POS"))
tmp <- tmp[!rownames(tmp) %in% rownames(locus.add),,drop=F]
locus.add <- rbind(locus.add, tmp)


tmp <- cnv.metastasis[!duplicated(cnv.metastasis[,c(1,4)]),]
for(i in 1:nrow(tmp)){
 readcounts.mut <- rbind(readcounts.mut, 0)
 readcounts.ref <- rbind(readcounts.ref, mean.cov)
}
rownames(readcounts.mut)[(nrow(readcounts.mut)-nrow(tmp) +1):nrow(readcounts.mut)] <- paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep=".")
rownames(readcounts.ref) <- rownames(readcounts.mut)
extra.muts <- c(extra.muts, paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."))
tmp <- tmp[,c(1,2)]
dimnames(tmp) <- list(paste(paste0("Chr",tmp$Chromosome), tmp$Start, sep="."), c("X.CHROM", "POS"))
tmp <- tmp[!rownames(tmp) %in% rownames(locus.add),,drop=F]
locus.add <- rbind(locus.add, tmp)


readcounts.mut <- readcounts.mut[unique(rownames(readcounts.mut)),]
readcounts.ref <- readcounts.ref[unique(rownames(readcounts.ref)),]
locus.add <- locus.add[unique(rownames(locus.add)),]
extra.muts <- unique(extra.muts)


coverage.ratios <- matrix(0, nrow(readcounts.mut), ncol(readcounts.mut), dimnames = list(rownames(readcounts.mut), colnames(readcounts.mut)))
baf <- matrix(NA, nrow(readcounts.mut), ncol(readcounts.mut), dimnames = list(rownames(readcounts.mut), colnames(readcounts.mut)))
locus <- unique(rbind(primary[,1:2], relapse[,1:2], metastasis[,1:2]))
locus <- rbind(locus, 
               matrix(c(2, mean(c(58373433, 95537299)),
                        3, mean(c(53321495, 90283666)),
                        3, mean(c(183873461, 184542551)),
                        5, 143122,
                        5, mean(c(1251287, 1297162)),
                        7, mean(c(27788121, 27886216)),
                        7, mean(c(150911230)),
                        10, mean(c(35896546, 35930313)),
                        10, mean(c(96101423, 104231144)),
                        10, mean(c(89621195, 89733687)),
                        11, mean(c(2966839,	8457769)),
                        16, mean(c(88039714, 88599311)),
                        17, mean(c(6003,	15938242)),
                        7, mean(c(21985250, 23830548)),
                        7, mean(c(64437249, 64864824)),
                        7, mean(c(102392135, 142478816))), byrow=T, ncol=2, 
                      dimnames = list( c("2p", "3p", "3q", "5p", "TERT", "7p", "7q", "10p", "10q", "PTEN", "11p", "16q", "17p", "7.1", "7.2", "7.3"),
                                       c("X.CHROM", "POS"))), locus.add)

rownames(locus) <- rownames(readcounts.mut)
  

## only keep loci with at least 30 read counts
to.keep <- apply(readcounts.mut + readcounts.ref, 1, function(x){
  if(all(x==0)){return(F)}
  if(max(x[x!=0])>=30){
    return(T)
  }else{
    return(F)
  }
})
to.keep[which(rownames(readcounts.mut)%in%extra.muts)] <- T

readcounts.mut <- readcounts.mut[to.keep,]
readcounts.ref <- readcounts.ref[to.keep,]
## extrapolate zero refs to the average in the other samples
readcounts.ref <- t(apply(readcounts.ref, 1, function(x){
  if(any(x==0)){
    x[x==0] <- round(mean(x[x!=0]))
  }
  x
}))
coverage.ratios <- coverage.ratios[rownames(readcounts.mut),]
baf <- baf[rownames(readcounts.mut),]
locus <- locus[rownames(readcounts.mut),]



for(i in 1:ncol(coverage.ratios)){
  coverage.ratios[,i] <- apply(locus, 1, function(x){
    tmp <- cnv.files[[i]][cnv.files[[i]]$Chromosome==as.character(x[1]),,drop=F]
    if(nrow(tmp)==0){return(1)}
    tmp <- tmp[which(tmp$Start <= as.numeric(x[2])),,drop=F]
    tmp <- tmp[nrow(tmp),,drop=F]#[1,,drop=F]
    if(nrow(tmp)==0 | is.na(tmp$MedianRatio)){return(1)}
    tmp$MedianRatio
  })
}

for(i in 1:ncol(coverage.ratios)){
  baf[,i] <- apply(locus, 1, function(x){
    tmp <- cnv.files[[i]][cnv.files[[i]]$Chromosome==as.character(x[1]),,drop=F]
    if(nrow(tmp)==0){return(NA)}
    if(nrow(tmp)==0){return(1)}
    tmp <- tmp[which(tmp$Start <= as.numeric(x[2])),,drop=F]
    tmp <- tmp[nrow(tmp),,drop=F]#[1,,drop=F]
    if(nrow(tmp)==0 | is.na(tmp$BAF) | abs(tmp$BAF) > 1 | tmp$BAF < 0){return(NA)}
    tmp$BAF
  })
}


## if baf is measured in only one sample, remove it from the info. If it's measured in two of three but not in the third, take the average of them for the third
baf <- t(apply(baf, 1, function(x){
  if(sum(is.na(x))==1){
    x[is.na(x)] <- mean(x[!is.na(x)])
  }else if(sum(is.na(x))==2){
    x[!is.na(x)] <- NA
  }
  x
}))

coverage.ratios <- replace(coverage.ratios, coverage.ratios < 0, 1)


readcounts.ref <- readcounts.ref[!is.na(rowSums(readcounts.mut)),,drop=F]
coverage.ratios <- coverage.ratios[!is.na(rowSums(readcounts.mut)) ,,drop=F]
readcounts.mut <- readcounts.mut[!is.na(rowSums(readcounts.mut)) ,,drop=F]
colnames(readcounts.ref) <-  LETTERS[1:(ncol(readcounts.ref))]
colnames(coverage.ratios) <- LETTERS[1:(ncol(readcounts.ref))]
colnames(readcounts.mut) <- LETTERS[1:(ncol(readcounts.ref))]


readcounts.mut <- readcounts.mut[which(rowMaxs(readcounts.mut) > 5 | rownames(readcounts.mut) %in% extra.muts) ,]
readcounts.ref <- readcounts.ref[rownames(readcounts.mut),]
coverage.ratios <- coverage.ratios[rownames(readcounts.mut),]
locus <- locus[rownames(readcounts.mut),]

## Manually re-check mutational input and correct
readcounts.mut["CDH22. 44803797",1] <- 56
readcounts.ref["CDH22. 44803797",1] <- 59
readcounts.mut["MAGI3.114165294",1] <- 41
readcounts.ref["MAGI3.114165294",1] <- 52
## throw this one out, due to low quality in primary/metastasis
readcounts.mut <- readcounts.mut[-which(rownames(readcounts.mut) %in% c("JHDM1D.139796225", "CSNK1G3.122926929")),]
readcounts.ref <- readcounts.ref[-which(rownames(readcounts.ref) %in% c("JHDM1D.139796225", "CSNK1G3.122926929")),]

readcounts.mut["RETSAT. 85572014",2] <- 5
readcounts.ref["RETSAT. 85572014",2] <- 70
readcounts.mut["CDH4. 59906707",2] <- 3
readcounts.ref["CDH4. 59906707",2] <- 10
readcounts.mut["LRP2.170082126",3] <- 14
readcounts.ref["LRP2.170082126",3] <- 24
readcounts.mut["PAK3.110435581",3] <- 5
readcounts.ref["PAK3.110435581",3] <- 33
readcounts.mut["MAGI3.114165294",3] <- 54
readcounts.ref["MAGI3.114165294",3] <- 51


all.mutations <- (readcounts.mut!=0)*1

unique.combis <- unique(all.mutations)

print(unique.combis)
## this reduces the tree to a pleasant level of complexity


readcounts.mut <- readcounts.mut[rownames(all.mutations),,drop=F]
readcounts.ref <- readcounts.ref[rownames(all.mutations),,drop=F]
coverage.ratios <- coverage.ratios[rownames(all.mutations),,drop=F]
baf <- baf[rownames(all.mutations),,drop=F]
locus <- locus[rownames(all.mutations),,drop=F]

readcounts.mut.extra <- readcounts.mut[rownames(readcounts.mut) %in% extra.muts,]
readcounts.ref.extra  <- readcounts.ref[rownames(readcounts.mut) %in% extra.muts,]
coverage.ratios.extra <- coverage.ratios[rownames(readcounts.mut) %in% extra.muts,]
locus.extra <- locus[rownames(readcounts.mut) %in% extra.muts,]


## try to find the most parsimonious tree. Start from the top, if clones are not explained by subclones, open up a new branch
unique.combis <- unique.combis[rowSums(unique.combis)>0,,drop=F]
trees.to.test <- list(0)
indicator <- list(0)
indispensable <- list(0)
for(random in 1:100){
  unique.combis <- unique.combis[sample(1:nrow(unique.combis), nrow(unique.combis), replace=F),]
  unique.combis <- unique.combis[order(rowSums(unique.combis), decreasing = T), ]
  tree <- unique.combis[1,,drop=F]
  tree.list <- list(tree)
  for(i in 2:nrow(unique.combis)){
    new.tree <- F
    for(j in 1:length(tree.list)){
      tree <- tree.list[[j]]
      to.compare <- which(rowSums(tree[,which(unique.combis[i,] != 0),drop=F]!=0)>0)
      if(any((t(tree[to.compare,,drop=F]) - unique.combis[i,]) < 0)){
        new.tree <- T
      }else{
        tree.list[[j]] <- rbind(tree.list[[j]], unique.combis[i,])
        new.tree <- F
        break
      }
    }
    if(new.tree){
      tree.list[[length(tree.list)+1]] <- unique.combis[i,,drop=F]
    }
  }
  
  tree <- tree.list[[1]]
  if(length(tree.list)>1){
    for(i in 2:length(tree.list)){
      tree <- cbind(tree, matrix(0, nrow= nrow(tree), ncol=ncol(tree.list[[i]])))
      tree <- rbind(tree, cbind(matrix(0, nrow= nrow(tree.list[[i]]), ncol= ncol(tree) - ncol(tree.list[[i]])), tree.list[[i]]))
    }
  }
  
  
  colnames(tree) <- make.unique(rep(colnames(tree.list[[1]]), length(tree.list)))
  sample.identifier <- rep(LETTERS[1:ncol(tree.list[[1]])], length(tree.list))
  sample.identifier <- sample.identifier[colSums(tree)>0]
  tree <- tree[,colSums(tree)>0]
  indispensable[[random]] <- 1:ncol(tree)
  indispensable[[random]] <- 1:ncol(readcounts.ref)
  ## add a clone for single mutation per subclone
  to.add <- matrix(0, ncol(tree), ncol(tree))
  diag(to.add) <- 1
  tree <- unique(rbind(tree, to.add))
  
  tree <- rbind(1, tree)
  
  tree <- tree[rowSums(tree)>0,,drop=F]
  
  trees.to.test[[random]] <- tree
  indicator[[random]] <- sample.identifier
}
## only test unique matrices. 
keep <- rep(T, length(trees.to.test))
for(i in 1:99){
  pseudo.tree.1 <- (sapply(unique(indicator[[i]]), function(x){
    rowSums(trees.to.test[[i]][,indicator[[i]]==x,drop=F])})!=0)*1
  for(j in (i+1):100){
    pseudo.tree.2 <- (sapply(unique(indicator[[j]]), function(x){
      rowSums(trees.to.test[[j]][,indicator[[j]]==x,drop=F])})!=0)*1
    if(sum(dim(pseudo.tree.1)==dim(pseudo.tree.2))==2){
      all.id <- T
      for(k in 1:nrow(pseudo.tree.1)){
        if(!any(colSums(t(pseudo.tree.2) == pseudo.tree.1[k,])==ncol(pseudo.tree.1))){
          all.id <- F
        }
      }
      if(all.id){
        keep[j] <- F
      }
    }
  }
}
for(i in 100:1){
  if(!keep[i]){
    trees.to.test[[i]] <- NULL
    indicator[[i]] <- NULL
    indispensable[[i]] <- NULL
  }
}

## expand these test trees to solutions of 2 subclones per sample
number.test.trees <- length(trees.to.test)
for(i in 1:number.test.trees){

  orders <- matrix(1:ncol(trees.to.test[[i]]), nrow=1)
  for(j in 2:length(orders)){
    orders <- rbind(orders, c(orders[nrow(orders),-1], orders[nrow(orders),1]))
  }
  for(k in 1:nrow(orders)){
    tmp <- trees.to.test[[i]]
    tmp.indicator <- indicator[[i]]
    for(j in orders[k,]){
      tmp <- cbind(tmp, tmp[,j])
      tmp <- rbind(tmp, 0)
      tmp[nrow(tmp), j] <- 1
      tmp <- rbind(tmp, 0)
      tmp[nrow(tmp), ncol(tmp)] <- 1
      tmp.indicator <- c(tmp.indicator, tmp.indicator[j])
      colnames(tmp) <- make.unique(tmp.indicator)
      trees.to.test[[length(trees.to.test)+1]] <- tmp
      indicator[[length(indicator)+1]] <- tmp.indicator
    }
  }


  ## and reductions thereof
}

## only test unique matrices. 
keep <- rep(T, length(trees.to.test))
number.test.trees <- length(trees.to.test)
for(i in 1:(number.test.trees-1)){
  tmp.indicator <- colnames(trees.to.test[[i]])
  for(j in (i+1):number.test.trees){
    tmp.indicator2 <- colnames(trees.to.test[[j]])
    if(all(tmp.indicator%in%tmp.indicator2) & all(tmp.indicator2%in%tmp.indicator)){
      keep[j] <- F
    }
  }
}
for(i in number.test.trees:1){
  if(!keep[i]){
    trees.to.test[[i]] <- NULL
    indicator[[i]] <- NULL
    indispensable[[i]] <- NULL
  }
}


## add three more trees that account for intermingling between A and C, the split between B and A,C, but extend to A,B - A,C, B,C - A,C, allowing with all extensions to get more complex solutions

trees.to.test <- c(trees.to.test, list(
  matrix(c(1,1,1,1,1,
           1,0,0,0,0,
           0,1,1,1,1,
           0,1,1,0,0,
           0,0,0,1,1,
           0,1,0,0,0,
           0,0,1,0,0,
           0,0,0,1,0,
           0,0,0,0,1), byrow = T, ncol=5, dimnames = list(c(), c("B", "A", "C", "A.1", "C.1"))),
  matrix(c(1,1,1,1,
           1,1,0,0,
           0,0,1,1,
           1,0,0,0,
           0,1,0,0,
           0,0,1,0,
           0,0,0,1), byrow = T, ncol=4, dimnames = list(c(), c("A", "B", "A.1", "C"))),
  matrix(c(1,1,1,1,
           1,1,0,0,
           0,0,1,1,
           1,0,0,0,
           0,1,0,0,
           0,0,1,0,
           0,0,0,1), byrow = T, ncol=4, dimnames = list(c(), c("B", "C", "A", "C.1")))
))

indicator <- c(indicator, list(c("B", "A", "C", "A", "C"), c("A", "B", "A", "C"), c("B", "C", "A", "C")))


## assume that samples are distributed over the same tree
trees.to.test.2 <- list(matrix(c(1), nrow=1, ncol=1), matrix(c(1,1,
                                                             1,0,
                                                             0,1), nrow=3, byrow=T),
                      matrix(c(1,1,1,
                               1,1,0,
                               1,0,0,
                               0,1,0,
                               0,0,1), byrow=T, nrow=5),
                      matrix(c(1,1,1,1,
                               1,1,0,0,
                               0,0,1,1,
                               1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1), ncol=4, byrow = T),
                      matrix(c(1,1,1,1,
                               1,1,1,0,
                               1,1,0,0,
                               1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,1), ncol=4, byrow=T),
                      matrix(c(1,1,1,1,1,
                               1,1,1,1,0,
                               1,1,1,0,0,
                               1,1,0,0,0,
                               1,0,0,0,0,
                               0,1,0,0,0,
                               0,0,1,0,0,
                               0,0,0,1,0,
                               0,0,0,0,1), ncol=5, byrow=T),
                      matrix(c(1,1,1,1,1,
                               1,1,1,0,0,
                               1,1,0,0,0,
                               0,0,0,1,1,
                               1,0,0,0,0,
                               0,1,0,0,0,
                               0,0,1,0,0,
                               0,0,0,1,0,
                               0,0,0,0,1), ncol=5, byrow=T),
                      matrix(c(1,1,1,1,1,
                               1,1,1,1,0,
                               1,1,0,0,0,
                               0,0,1,1,0,
                               1,0,0,0,0,
                               0,1,0,0,0,
                               0,0,1,0,0,
                               0,0,0,1,0,
                               0,0,0,0,1), ncol=5, byrow=T))

## expand to samples
indicator.2 <- list(0)
indispensable.2 <- list(0)
for(i in 1:length(trees.to.test.2)){
  trees.to.test.2[[i]] <- trees.to.test.2[[i]][,rep(1:ncol(trees.to.test.2[[i]]), each=ncol(unique.combis)),drop=F]
  colnames(trees.to.test.2[[i]]) <- rep(colnames(unique.combis), ncol(trees.to.test.2[[i]])/ncol(unique.combis))
  colnames(trees.to.test.2[[i]]) <- make.unique(colnames(trees.to.test.2[[i]]))
  indicator.2[[i]] <- rep(colnames(unique.combis), ncol(trees.to.test.2[[i]])/ncol(unique.combis))
  for(j in 1:ncol(trees.to.test.2[[i]])){
    trees.to.test.2[[i]] <- rbind(trees.to.test.2[[i]], 0)
    trees.to.test.2[[i]][nrow(trees.to.test.2[[i]]),j] <- 1
  }
  
}
tree.from.muts <- 1:length(trees.to.test)
## I want to get extensions of the very first and the 2 last trees to test from muts
tree.from.muts <- tree.from.muts[-c(1, length(tree.from.muts)-1, length(tree.from.muts))]

trees.to.test <- c(trees.to.test, trees.to.test.2)
inidispensable <- c(indispensable, indispensable.2)
indicator <- c(indicator, indicator.2)

## require the following changes to be joined: Chr2: start-end
## Chr3: start - end
## Chr4
joined.info <- as.list(rep(NA, length=nrow(readcounts.mut)))
for(i in 1:nrow(locus)){
  if(locus[i,1]==2 & sum(coverage.ratios[i,]==coverage.ratios["2p",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==2 & sum(coverage.ratios[i,]==coverage.ratios["BOC.113003315",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==3 & sum(coverage.ratios[i,]==coverage.ratios["3p",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==3 & sum(coverage.ratios[i,]==coverage.ratios["3q",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==5 & sum(coverage.ratios[i,]==coverage.ratios["5p",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==17 & sum(coverage.ratios[i,]==coverage.ratios["TERT",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==11 & sum(coverage.ratios[i,]==coverage.ratios["11p",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==16 & sum(coverage.ratios[i,]==coverage.ratios["16q",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }else if(locus[i,1]==17 & sum(coverage.ratios[i,]==coverage.ratios["17p",])==ncol(coverage.ratios)){
    joined.info[[i]] <- list(c("A", "B", "C"))
  }
}





  
  ## fit only on loci carrying a mutation
  subclones <- model.fits(tree.list = trees.to.test, readcounts.mut = readcounts.mut, 
                          readcounts.ref = readcounts.ref, 
                          coverage.ratio = coverage.ratios, 
                          baf = baf,
                          locus = locus, 
                          sex = sex, directory = output.path, sample.indicator.list = indicator,
                          normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, indispensable = indispensable,
                          normal.b.allele.norm = normal.b.allele.norm, trees.from.muts = tree.from.muts, extra.muts = extra.muts, joined.info = joined.info,
                          false.negatives = T)
  save(subclones, file=paste0(output.path, "/subclones.RData"))


  if(length(subclones[[1]])%in%c(0,1)){
    subclones[[1]] <- NULL
  }
  ## first thing we have to do is to find the most likely solution
  
    BIC <- rep(Inf, length(subclones))
    for(i in 1:length(subclones)){
      print(i)
      k <- which.min(subclones[[i]]$Likelihood)
      tree. <- subclones[[i]]$Tree
      indicator. <- subclones[[i]]$indicator
      parms. <- subclones[[i]]$Parameters[k,]
      
      tree. <- tree.[rowSums(tree.)>0,,drop=F]
      tree. <- unique(tree.)
      Tau <- Expectation(tree = tree., read.count.matrix.mut = readcounts.mut[!(rownames(readcounts.mut) %in% extra.muts),], 
                         read.count.matrix.ref = readcounts.ref[!(rownames(readcounts.mut) %in% extra.muts),],
                         parms = parms., 
                         coverage.ratio = coverage.ratios[!(rownames(readcounts.mut) %in% extra.muts),], baf = baf, locus = locus[!(rownames(readcounts.mut) %in% extra.muts),], sex = sex, 
                         sample.indicator = indicator.,
                         normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, 
                         joined.info = joined.info[!(rownames(readcounts.mut) %in% extra.muts)], false.negatives = T)
      
      BIC[i] <- BIC(-min(subclones[[i]]$Likelihood), n.parms = ncol(tree.), 
                              n.data = sum(readcounts.mut[!(rownames(readcounts.mut) %in% extra.muts),] + readcounts.ref[!(rownames(readcounts.mut) %in% extra.muts),]), 
                              Tau = Tau$number.of.possibilities, gamma = 0.9)
      
    }
 
  
  

  cand.models <- which(BIC < min(BIC) + 10)
  for(k in cand.models){
    kk <- which.min(subclones[[k]]$Likelihood)
    ## plot the tree along with the mutations per branch. Annotate the samples
    parms <- subclones[[k]]$Parameters[kk,]
    tree <- unique(subclones[[k]]$Tree)
    indicator <- subclones[[k]]$indicator
    
   
    readcounts.mut. <- readcounts.mut
    readcounts.ref. <- readcounts.ref
    locus. <- locus
    coverage.ratios. <- coverage.ratios
    
    to.keep <- readcounts.mut. + readcounts.ref.
    to.keep <- which(apply(to.keep, 1, function(x){ifelse(mean(x[x > 0])>30, T, F)}))
    readcounts.mut. <- readcounts.mut.[to.keep,,drop=F]
    readcounts.ref. <- readcounts.ref.[to.keep,]
    coverage.ratios. <- coverage.ratios.[to.keep,]
    locus. <- locus.[to.keep,,drop=F]
    joined.info. <- joined.info
    for(i in length(joined.info):1){
      if(!i %in% to.keep){
        joined.info.[[i]] <- NULL
      }
    }

    res <- Expectation(tree = tree, read.count.matrix.mut = readcounts.mut., read.count.matrix.ref = readcounts.ref.,
                       parms = parms, 
                       coverage.ratio = coverage.ratios., baf = baf, locus = locus., sex = sex, 
                       sample.indicator = indicator,
                       normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, joined.info = joined.info., false.negatives = T)
    
    if(any(parms==0)){
      res$readcounts <- res$readcounts[,-c(which(parms==0)*2 - 1, which(parms==0)*2),drop=F]
      res$norm.ploidy <- res$norm.ploidy[,-c(which(parms==0)*2 - 1, which(parms==0)*2),drop=F]
      res$ploidy <- res$ploidy[,-c(which(parms==0)*2 - 1, which(parms==0)*2),drop=F]
      tree <- tree[,-which(parms==0),drop=F]
      indicator <- indicator[-which(parms==0)]
      parms <- parms[-which(parms==0)]
    }
    


    #### Tree plotting
    library(igraph)
    tree. <- cbind(tree, 0)
    colnames(tree.) <- c(colnames(tree), "N")
    
    tree.for.plotting <- Manually.Infer.Tree(list(centers=tree.))
    
    pdf(paste0(output.path,"/cladogram_", which(cand.models==k), ".pdf"), useDingbats = F)
    
    plot(collapse.singles(tree.for.plotting), "cladogram", edge.width=2,  cex=1)
    
    parms.wo.normal <- parms
    for(i in unique(indicator)){
      parms.wo.normal[which(indicator==i)] <- parms[which(indicator==i)]/sum(parms[which(indicator==i)])
    }
    
    
    g <- graph( as.vector(t(tree.for.plotting$edge)), directed = FALSE )
    plot.igraph(g, 
                vertex.label = c(tree.for.plotting$tip.label, rep("", tree.for.plotting$Nnode+3)),
                vertex.shape = "circle",  vertex.label.cex = 1, edge.width = 4,
                layout=layout.reingold.tilford(g,root=length(tree.for.plotting$tip.label)+1))
    
    
 
    
    edge.labels <- list()
    edge.lengths <- c()
    for(jj in 1:nrow(tree)){
      print(jj)
      edge.labels[[jj]] <- rownames(res$norm.ploidy)[which(colSums(t(res$norm.ploidy[,seq(2,ncol(res$norm.ploidy), 2)]!=0)==tree[jj,])==ncol(tree) &
                                                             rowSums(res$readcounts[,seq(2,ncol(res$readcounts), 2)])!=0)]
      
      edge.lengths <- c(edge.lengths, length(edge.labels[[jj]]))
      

      
      involved.nodes <- paste(colnames(tree)[which(tree[jj,,drop=F]!=0)], collapse = ".")
      
      
      names(edge.lengths)[length(edge.lengths)] <- involved.nodes
      
      
    }
    names(edge.labels) <- names(edge.lengths)
    
    edge.lengths <- edge.lengths[tree.for.plotting$node.label[tree.for.plotting$node.label!="O"]]
    edge.lengths <- replace(edge.lengths, which(is.na(edge.lengths)),0)
    
    edge.labels <- edge.labels[tree.for.plotting$node.label[tree.for.plotting$node.label!="O"]]
    edge.labels <- replace(edge.labels, which(is.na(edge.labels)),"")
    
    for(ii in unique(indicator)){
      parms.tmp <- round(parms[which(indicator==ii)]/sum( parms[which(indicator==ii)]), digits = 2)*100
      colnames(tree.)[which(indicator==ii)] <- paste(ii, parms.tmp, "%")
    }
    
    tree.for.plotting <- Manually.Infer.Tree(list(centers=tree.))
    tree.for.plotting$edge.length <- edge.lengths
    tree.for.plotting$edge.labels <- edge.labels
    
    plot(collapse.singles(tree.for.plotting))
    
    
    
    
    plot.igraph(g, edge.label = edge.labels,
                vertex.label = c(tree.for.plotting$tip.label, rep("", tree.for.plotting$Nnode+3)),
                vertex.shape = "circle",  vertex.label.cex = 1, edge.width = 4,
                layout=layout.reingold.tilford(g,root=length(tree.for.plotting$tip.label)+1))
    
    
    dev.off()
    
    
  }
