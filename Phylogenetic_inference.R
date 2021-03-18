library(phangorn)
library(foreach)
library(doParallel)




#' Negative log-likelihood
#' 
#' This function computes the negative log-likelihood for a solution based on a multinomial distribution
#' @param readcounts a readcount table as returned by Expectation
#' @param copy.number the copy numbers as returened by Expectation
#' @param n.clones the number of subclones (tumor plus normal tissue)
#' @return The maximum likelihood estimate of the relative subclone sizes
#' @export
## This function computes the negative log-likelihood for a solution based on a multinomial distribution

NegLogLikelihood <- function(readcounts, copy.number, parms){
  
  L <- 0
  
  for(k in 1:nrow(readcounts)){
    if(sum(copy.number[k,]*parms)==0){
      return(Inf)
    }else{
      mut.cols <- seq(2, ncol(readcounts), 2)
      prob.mut <- sum((copy.number[k, ]*parms)[mut.cols])
      prob.mut <- replace(prob.mut, prob.mut>1, 1)
      prob.mut <- replace(prob.mut, prob.mut<0, 0)
      mut.counts <- sum(readcounts[k, mut.cols])
      L <- L -dbinom(mut.counts, sum(readcounts[k,]), prob.mut, log=T)
    }
  }
  return(L)
}




#' Randomize initial parameters
#' 
#' Sample n parameters between 0 and 1 that add up to 1
#' @param n The number of parameters
#' @return The randomized relative sizes
#' @export

Rand.parms <- function(n){
  prob <- runif(n,0,1)
  prob <- prob/sum(prob)
  return(prob)
}




#' Expected read counts per subclone at a given phylogenetic structure
#' 
#' This function calculates the expected read counts per subclone given a subclonal structure (an initial guess of the subclonal proportions must be provided based on which the expected read counts can be calculated). This is based on the expected values of multinomial distributions for reference and mutated read counts. 
#' @param tree a binary tree matrix as returned by Tree.matrix.R
#' @param read.count.matrix.prim a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the primary tumor
#' @param read.count.matrix.rec a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the relapse tumor
#' @param coverage.ratio.baf.prim a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the primary tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param coverage.ratio.baf.rec a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the relapse tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param which.prim vector indicating which subclones are primary subclones
#' @param which.rec vector indicating which subclones are relapse subclones
#' @param locus vector containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param parms.prim vector containing the relative subclone sizes of the primary tumor
#' @param parms.rec vector containing the relative subclone sizes of the relapse tumor
#' @param likelihood.space TRUE if additional solutions should be reported. I.e. if a SNV/CNV cannot be uniquely mapped to a specific locus in the tree, additional solutions are reported. In addition, residual sums of squares between measured and estimated variant allele frequencies, B allele frequencies and coverage ratios are reported.
#' @return Two matrices storing the expected readcounts (reference, mutated) at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones; columns represent subclones; last two columns represent normal tissue); two matrices storing the relative number of reference and mutated alleles at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones); two matrices storing the number of alleles at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones); two matrices indicating which subclones carry a CNV at each mutated locus per tumor subclone. One integer reporting the logarithmic number of possibilities to sort the mutations on the tree. If likelihood.space==T, additional information is provided: additional solutions for copy number changes are stored in separate matrices; additional solutions for SNVs are stored in the same matrices as before (they will be longer than and affected SNVs can be found by checkin which SNV name appears twice among the rownames of the matrix). In addition, residual sums of squares for the variant allele frequencies, the B allele frequences and the coverage ratios are returned
#' @export


Expectation <- function(tree, read.count.matrix.prim, read.count.matrix.rec, parms.prim, parms.rec, coverage.ratio.baf.prim, coverage.ratio.baf.rec, locus, sex = sex,
                        which.prim = c(1,3,5), which.rec=c(2,4,6), likelihood.space = F){
  
  ## distinguish shared and private mutations
  common.genes <- intersect(rownames(read.count.matrix.prim), rownames(read.count.matrix.rec))
  genes.prim <- rownames(read.count.matrix.prim)[!rownames(read.count.matrix.prim) %in% rownames(read.count.matrix.rec)]
  genes.rec <- rownames(read.count.matrix.rec)[!rownames(read.count.matrix.rec) %in% rownames(read.count.matrix.prim)]
  
  read.count.matrix.prim <- read.count.matrix.prim[c(common.genes, genes.prim), ,drop=F]
  read.count.matrix.rec <- read.count.matrix.rec[c(common.genes, genes.rec), ,drop=F]
  coverage.ratio.prim <- coverage.ratio.baf.prim[c(rownames(read.count.matrix.prim), genes.rec),,drop=F]
  if(ncol(coverage.ratio.baf.prim)==1){
    baf.prim <- rep(NA, nrow(coverage.ratio.prim))
  }else{
    baf.prim <- coverage.ratio.baf.prim[c(rownames(read.count.matrix.prim), genes.rec),2]
  }
  coverage.ratio.prim <- coverage.ratio.prim[,1,drop=F]
  
  
  
  coverage.ratio.rec <- coverage.ratio.baf.rec[c(rownames(read.count.matrix.prim), genes.rec),,drop=F]
  if(ncol(coverage.ratio.baf.rec)==1){
    baf.rec <- rep(NA, nrow(coverage.ratio.rec))
  }else{
    baf.rec <- coverage.ratio.baf.rec[c(rownames(read.count.matrix.prim), genes.rec),2]
  }
  coverage.ratio.rec <- coverage.ratio.rec[,1,drop=F]
  locus <- locus[c(rownames(read.count.matrix.prim), genes.rec), ,drop=F]
  res.prim <- matrix(0, nrow=nrow(read.count.matrix.prim), ncol=length(parms.prim)*2)
  rownames(res.prim) <- rownames(read.count.matrix.prim)
  colnames(res.prim) <- rep(c("reference", "mutated"), ncol(res.prim)/2)
  s.res.prim <- res.prim
  p.res.prim <- res.prim
  p.res.prim.index <- matrix(0, nrow=length(c(common.genes, genes.prim, genes.rec)), ncol=length(parms.prim)*2)
  rownames(p.res.prim.index) <- c(common.genes, genes.prim, genes.rec)
  colnames(p.res.prim.index) <- rep(c("reference", "mutated"), ncol(res.prim)/2)
  
  res.rec <- matrix(0, nrow=nrow(read.count.matrix.rec), ncol=length(parms.rec)*2)
  rownames(res.rec) <- rownames(read.count.matrix.rec)
  colnames(res.rec) <- rep(c("reference", "mutated"), ncol(res.rec)/2)
  s.res.rec <- res.rec
  p.res.rec <- res.rec
  p.res.rec.index <- matrix(0, nrow=length(c(common.genes, genes.prim, genes.rec)), ncol=length(parms.rec)*2)
  rownames(p.res.rec.index) <- c(common.genes, genes.prim, genes.rec)
  colnames(p.res.rec.index) <- rep(c("reference", "mutated"), ncol(res.rec)/2)
  
  cnv.indicator.all <- rep(0, nrow(p.res.prim.index))
  names(cnv.indicator.all) <- rownames(p.res.prim.index)
  
  ## if we are not only interested in a point estimate, but in the likelihood space, store the likelihood for each locus separately
  if(likelihood.space){
    like.per.loc.prim <- rep(Inf, nrow(p.res.prim.index))
    like.per.loc.rec <- rep(Inf, nrow(p.res.prim.index))
    ## store additional copy numbers in the case of ambuiguity
    additional.copy.number.prim <- c()
    additional.copy.number.rec <- c()
    ## store residuals of coverage ratios, VAFs and BAFs
    residuals.coverage.ratio.prim <- rep(NA, nrow(p.res.prim.index))
    residuals.bafs.prim <- rep(NA, nrow(p.res.prim.index))
    residuals.vafs.prim <- rep(NA, nrow(p.res.prim.index))
    residuals.coverage.ratio.rec <- rep(NA, nrow(p.res.prim.index))
    residuals.bafs.rec <- rep(NA, nrow(p.res.prim.index))
    residuals.vafs.rec <- rep(NA, nrow(p.res.prim.index))
    
    names(residuals.bafs.prim) <- c(common.genes, genes.prim, genes.rec)
    names(residuals.vafs.prim) <- c(common.genes, genes.prim, genes.rec)
    names(residuals.coverage.ratio.prim) <- c(common.genes, genes.prim, genes.rec)
    names(residuals.coverage.ratio.rec) <- c(common.genes, genes.prim, genes.rec)
    names(residuals.bafs.rec) <- c(common.genes, genes.prim, genes.rec)
    names(residuals.vafs.rec) <- c(common.genes, genes.prim, genes.rec)
  }
  
  ## If all tumor subclones are of zero size, abort the algorithm
  if(sum(parms.prim[-length(parms.prim)])==0){
    print("PT_impurity")
    return("PT_impurity")}
  if(sum(parms.rec[-length(parms.rec)])==0){
    print("RT_impurity")
    return("RT_impurity")}
  
  ## store the number of possible combinations per locus
  number.of.possibilities <- log(1)
  
  ## iterate through x genes. Each likelihood can be independently maximized (-logL minimized)
  
  for(i in seq_len((length(common.genes)+length(genes.prim)+length(genes.rec)))){
    ## we need a second iterator for the relapse tumor
    i.rec <- i
    
    
    ## the first thing we have to do is finding out what copy numbers we have to choose. We assume that evolution is simple:
    ## if the coverage ratio estimate indicates a deletion, we test combinations for copy number = 0 and copy number = 1 and choose the most likely one
    ## if the estimate indicates a duplication/ gain, we test combinations for copy number > 2
    ## if both, primary and relapse have a duplication/deletion, we test whether a common event or two separate events are more likely
    
    ## check whether you have to estimate the copy number (if it's the same as in the last segment you can reuse it)
    determine.copy.number <- T
    if(i > 1){
      if(coverage.ratio.prim[i,]==coverage.ratio.prim[i-1,] && coverage.ratio.rec[i,]==coverage.ratio.rec[i-1,] && locus[i,1]==locus[i-1,1] &&
         !(baf.prim[i]==baf.prim[i-1]) %in% c(F, NA) && 
         !(baf.rec[i]==baf.rec[i-1]) %in% c(F, NA)){
        determine.copy.number <- F
      }
    }
    
    ## the normal copy number /BAF is different on male sex chromosomes
    if(sex== "male" & locus[i,1] %in% c("X", "Y")){
      normal.copy.number <- 1
      normal.baf <- 0
    }else{
      normal.copy.number <- 2
      normal.baf <- 1
    }
    
    
    
    
    if(determine.copy.number){
      ## assume coverage ratios of 1 +- 10 % as normal
      
      if(coverage.ratio.prim[i,] < 1.1 && coverage.ratio.rec[i,] < 1.1 & coverage.ratio.prim[i,] > 0.9 && coverage.ratio.rec[i,] > 0.9 ){
        if(normal.copy.number == 2 && (is.na(baf.prim[i]) || baf.prim[i]<0.55 && baf.prim[i] > 0.45) && (is.na(baf.rec[i]) || baf.rec[i]<0.55 && baf.rec[i] > 0.45)){
          copy.number.prim <- normal.copy.number
          copy.number.rec <- normal.copy.number
          b.allele.prim <- normal.baf
          b.allele.rec <- normal.baf
          which.prim.cnv <- NA
          which.rec.cnv <- NA
          cnv.indicator <- "joined"
        }else if(normal.copy.number == 1 && (is.na(baf.prim[i]) || baf.prim[i]<0.05) && (is.na(baf.rec[i]) || baf.rec[i]<0.05)){
          copy.number.prim <- normal.copy.number
          copy.number.rec <- normal.copy.number
          b.allele.prim <- normal.baf
          b.allele.rec <- normal.baf
          which.prim.cnv <- NA
          which.rec.cnv <- NA
          cnv.indicator <- "joined"
        }else{
          
          inferred.copy.numbers <- InferCopyNumber(tree = tree, coverage.ratio.prim = coverage.ratio.prim[i,], coverage.ratio.rec = coverage.ratio.rec[i,], which.prim = which.prim, which.rec = which.rec,
                                                   locus = locus[i,], sex = sex, baf.prim = baf.prim[i], baf.rec = baf.rec[i], parms.prim, parms.rec)
          
          
          ## Choose the most likely solution:
          if(inferred.copy.numbers$joined[[7]] < inferred.copy.numbers$separated[[7]]){
            inferred.copy.numbers <- inferred.copy.numbers$joined
          }else{
            if(likelihood.space){
              ## if we are interested in the solution space, we report both, the joined and the separate solution if the logL of the joined solution is at least 50 % of the logL of the separate solution (stored in additional solutions)
              if(inferred.copy.numbers$joined[[7]] / inferred.copy.numbers$separated[[7]] <= 2){
                additional.copy.number.prim <- rbind(additional.copy.number.prim, inferred.copy.numbers$joined[[5]]*inferred.copy.numbers$joined[[1]])
                additional.copy.number.rec <- rbind(additional.copy.number.rec, inferred.copy.numbers$joined[[6]]*inferred.copy.numbers$joined[[2]])
                rownames(additional.copy.number.prim)[nrow(additional.copy.number.prim)] <- rownames(coverage.ratio.prim)[i]
                rownames(additional.copy.number.rec)[nrow(additional.copy.number.rec)] <- rownames(coverage.ratio.prim)[i]
                cnv.indicator.all <- c(cnv.indicator.all, 0)
                names(cnv.indicator.all)[length(cnv.indicator.all)] <- rownames(coverage.ratio.prim)[i]
              }
            }
            inferred.copy.numbers <- inferred.copy.numbers$separated
          }
          
          copy.number.prim <- inferred.copy.numbers[[1]]
          copy.number.rec <- inferred.copy.numbers[[2]]
          b.allele.prim <- inferred.copy.numbers[[3]]
          b.allele.rec <- inferred.copy.numbers[[4]]
          which.prim.cnv <- inferred.copy.numbers[[5]]
          which.rec.cnv <- inferred.copy.numbers[[6]]
          cnv.indicator <- inferred.copy.numbers[[8]]
        }
        
      }else{
        inferred.copy.numbers <- InferCopyNumber(tree = tree, coverage.ratio.prim = coverage.ratio.prim[i,], coverage.ratio.rec = coverage.ratio.rec[i,], which.prim = which.prim, which.rec = which.rec,
                                                 sex = sex, baf.prim = baf.prim[i], baf.rec = baf.rec[i], locus = locus[i,], parms.prim = parms.prim, parms.rec = parms.rec)
        
        
        if(inferred.copy.numbers$joined[[7]] < inferred.copy.numbers$separated[[7]]){
          inferred.copy.numbers <- inferred.copy.numbers$joined
        }else{
          if(likelihood.space){
            if(inferred.copy.numbers$joined[[7]] / inferred.copy.numbers$separated[[7]] <= 2){
              additional.copy.number.prim <- rbind(additional.copy.number.prim, inferred.copy.numbers$joined[[5]]*inferred.copy.numbers$joined[[1]])
              additional.copy.number.rec <- rbind(additional.copy.number.rec, inferred.copy.numbers$joined[[6]]*inferred.copy.numbers$joined[[2]])
              rownames(additional.copy.number.prim)[nrow(additional.copy.number.prim)] <- rownames(coverage.ratio.prim)[i]
              rownames(additional.copy.number.rec)[nrow(additional.copy.number.rec)] <- rownames(coverage.ratio.prim)[i]
            }
          }
          inferred.copy.numbers <- inferred.copy.numbers$separated
        }
        
        copy.number.prim <- inferred.copy.numbers[[1]]
        copy.number.rec <- inferred.copy.numbers[[2]]
        b.allele.prim <- inferred.copy.numbers[[3]]
        b.allele.rec <- inferred.copy.numbers[[4]]
        which.prim.cnv <- inferred.copy.numbers[[5]]
        which.rec.cnv <- inferred.copy.numbers[[6]]
        cnv.indicator <- inferred.copy.numbers[[8]]
      }
    }else{
      if(likelihood.space && length(additional.copy.number.prim)>0){
        if(rownames(coverage.ratio.prim)[i-1] %in% rownames(additional.copy.number.prim)){
          additional.copy.number.prim <- rbind(additional.copy.number.prim, additional.copy.number.prim[nrow(additional.copy.number.prim),])
          additional.copy.number.rec <- rbind(additional.copy.number.rec, additional.copy.number.rec[nrow(additional.copy.number.rec),])
          rownames(additional.copy.number.prim)[nrow(additional.copy.number.prim)] <- rownames(coverage.ratio.prim)[i]
          rownames(additional.copy.number.rec)[nrow(additional.copy.number.rec)] <- rownames(coverage.ratio.prim)[i]
        }
      }
    }
    
    
    
    if(cnv.indicator=="separated"){
      cnv.indicator.all[i] <- 1
    }
    
    
    if(likelihood.space){
      ## difference between estimated subclone sizes from entire tree and estimated subclone size at locus
      if(!is.na(which.prim.cnv[1])){
        
        fit.cell.fraction <- sum(parms.prim[which(which.prim.cnv!=0)])
        data.cell.fraction.1 <- (ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)*(coverage.ratio.prim[i,]-1))/(copy.number.prim - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
        data.cell.fraction.2 <- (baf.prim[i]*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2) - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1))/(b.allele.prim - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1) - baf.prim[i]*copy.number.prim + baf.prim[i]*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
        
        residuals.coverage.ratio.prim[i] <- fit.cell.fraction - mean(c(data.cell.fraction.1, data.cell.fraction.2))
        residuals.bafs.prim[i] <- fit.cell.fraction - data.cell.fraction.2
        
        residuals.coverage.ratio.prim[i] <- coverage.ratio.prim[i,] - (fit.cell.fraction*copy.number.prim + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))/
          ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)
        residuals.bafs.prim[i] <- baf.prim[i] - (fit.cell.fraction*b.allele.prim + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1))/
          (fit.cell.fraction*copy.number.prim + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
        
      }
      if(!is.na(which.rec.cnv[1])){
        fit.cell.fraction <- sum(parms.rec[which(which.rec.cnv!=0)])
        data.cell.fraction.1 <- (ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)*(coverage.ratio.rec[i,]-1))/(copy.number.rec - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
        data.cell.fraction.2 <- (baf.rec[i]*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2) - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1))/(b.allele.rec - ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1) - baf.rec[i]*copy.number.rec + baf.rec[i]*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
        
        residuals.coverage.ratio.rec[i] <- coverage.ratio.rec[i.rec,] - (fit.cell.fraction*copy.number.rec + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))/
          ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)
        residuals.bafs.rec[i] <- baf.rec[i.rec] - (fit.cell.fraction*b.allele.rec + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 0, 1))/
          (fit.cell.fraction*copy.number.rec + (1-fit.cell.fraction)*ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2))
      }
      
    }
    
    ## we don't need to determine the matrix of candidate solutions for a mutation distribution if the previous mutation was present in the same subset of the tumors and was at a locus of equal copy number state
    determine.s.matrix <- T
    if(i > 1){
      if(determine.copy.number==F){
        if(i <= length(common.genes) || (i > (length(common.genes) + 1) && 
                                         i <= (length(common.genes) + length(genes.prim))) || 
           (i > (length(common.genes) + length(genes.prim) + 1))){
          determine.s.matrix <- F
          if(i > (length(common.genes) + length(genes.prim) + 1)){
            i.rec <- i - length(genes.prim)
          }
        }
      }
    }
    
    
    if(determine.s.matrix){
      ## mutations that are present in both tumors require a different set of candidate solutions than mutations which are present in one tumor only
      if(i <= length(common.genes)){
        
        s.p.prim.rec <- Index.Matrix(tree = tree, copy.number.prim = copy.number.prim, copy.number.rec = copy.number.rec, which.prim = which.prim, which.rec = which.rec, mut.prim = T, mut.rec = T, locus = locus[i,],
                                     sex = sex, b.allele.prim = b.allele.prim, b.allele.rec = b.allele.rec, cnv.indicator = cnv.indicator, which.prim.cnv = which.prim.cnv, which.rec.cnv = which.rec.cnv)
        
        s.prim <- s.p.prim.rec$s.prim
        s.rec <- s.p.prim.rec$s.rec
        p.prim <- s.p.prim.rec$p.prim
        p.rec <- s.p.prim.rec$p.rec
        rm(s.p.prim.rec)
        
        
        
        ## now, only take solutions where the rows are in agreement with the previously found best solutions:
        
        if(!is.na(which.prim.cnv[1])& copy.number.prim != ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)&
           b.allele.prim !=1){
          rows.to.keep.prim <- which(colSums(t((p.prim!=0)*1) == which.prim.cnv)==(length(which.prim)))
        }else{
          rows.to.keep.prim <- seq(1:nrow(p.prim))
        }
        if(!is.na(which.rec.cnv[1]) & copy.number.rec != ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)&
           b.allele.rec !=1){
          rows.to.keep.rec <- which(colSums(t((p.rec!=0)*1) == which.rec.cnv)==(length(which.rec)))
        }else{
          rows.to.keep.rec <- seq(1:nrow(p.rec))
        }
        rows.to.keep <- intersect(rows.to.keep.prim, rows.to.keep.rec)
        
        
        p.prim <- p.prim[rows.to.keep,,drop=F]
        p.rec <- p.rec[rows.to.keep,,drop=F]
        s.prim <- s.prim[rows.to.keep,,drop=F]
        s.rec <- s.rec[rows.to.keep,,drop=F]
        
        ## now, add normal tissue to both, primary and relapse
        s.prim <- cbind(s.prim, 0)
        p.prim <- cbind(p.prim, 0)
        s.rec <- cbind(s.rec, 0)
        p.rec <- cbind(p.rec, 0)
        
        ## remove rows where no mutations are:
        rows.no.mut.prim <- which(rowSums(s.prim[,parms.prim!=0,drop=F])==0)
        rows.no.mut.rec <- which(rowSums(s.rec[,parms.rec!=0,drop=F])==0)
        if(length(rows.no.mut.prim)>0 || length(rows.no.mut.rec)>0){
          s.prim <- s.prim[-unique(c(rows.no.mut.prim, rows.no.mut.rec)),,drop=F]
          s.rec <- s.rec[-unique(c(rows.no.mut.prim, rows.no.mut.rec)),,drop=F]
          p.prim <- p.prim[-unique(c(rows.no.mut.prim, rows.no.mut.rec)),,drop=F]
          p.rec <- p.rec[-unique(c(rows.no.mut.prim, rows.no.mut.rec)),,drop=F]
        }
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.prim[,parms.prim!=0,drop=F], s.rec[,parms.rec!=0,drop=F], p.prim[,parms.prim!=0,drop=F], p.rec[,parms.rec!=0,drop=F]))
        s.prim <- s.prim[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        s.rec <- s.rec[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        
        ## remove zero rows based on the clones that are occupied
        zero.rows <- rowSums(cbind(s.prim[,parms.prim!=0,drop=F],s.rec[,parms.rec!=0,drop=F]))==0
        s.prim <- s.prim[!zero.rows,,drop=F]
        s.rec <- s.rec[!zero.rows,,drop=F]
        p.prim <- p.prim[!zero.rows,,drop=F]
        p.rec <- p.rec[!zero.rows,,drop=F]
        
        ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
        p.prim.index <- p.prim
        if(copy.number.prim == normal.copy.number && b.allele.prim == normal.baf){
          p.prim.index[,seq_len(ncol(p.prim))] <- 0
        }
        p.rec.index <- p.rec
        if(copy.number.rec == normal.copy.number && b.allele.rec == normal.baf){
          p.rec.index[,seq_len(ncol(p.rec.index))] <- 0
        }
        
        ## now, transform p.prim and p.rec into copy number table
        p.prim <- replace(p.prim, p.prim!=0, copy.number.prim)
        p.prim <- replace(p.prim, p.prim.index==0, normal.copy.number)
        if(copy.number.prim==0){
          p.prim <- replace(p.prim, p.prim.index==1, 0)
        }
        p.rec <- replace(p.rec, p.rec!=0, copy.number.rec)
        p.rec <- replace(p.rec, p.rec.index==0, normal.copy.number)
        if(copy.number.rec==0){
          p.rec <- replace(p.rec, p.rec.index==1, 0)
        }
        
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.prim[,parms.prim!=0,drop=F], s.rec[,parms.rec!=0,drop=F], p.prim[,parms.prim!=0,drop=F], p.rec[,parms.rec!=0,drop=F]))
        s.prim <- s.prim[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        s.rec <- s.rec[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        p.prim.index <- p.prim.index[!dupl.rows,,drop=F]
        p.rec.index <- p.rec.index[!dupl.rows,,drop=F]
        
        ## add normal tissue to each clone in the index matrix
        p.prim.index <- p.prim.index[,rep(seq_len(ncol(p.prim.index)),each=2),drop=F]
        p.rec.index <- p.rec.index[,rep(seq_len(ncol(p.rec.index)),each=2),drop=F]
        s.prim <- s.prim/p.prim
        s.rec <- s.rec/p.rec
        s.prim <- replace(s.prim, is.infinite(s.prim), 0)
        s.rec <- replace(s.rec, is.infinite(s.rec),0)
        s.prim <- replace(s.prim, is.na(s.prim), 0)
        s.rec <- replace(s.rec, is.na(s.rec),0)
        
        
        ## now, add the normal counterparts for each subclone and divide by the respective copy number
        s.prim <- s.prim[,rep(seq_len(ncol(s.prim)), each=2),drop=F]
        for(jj in seq(2, ncol(s.prim), 2)){
          s.prim[,jj-1] <- 1-s.prim[,jj]
        }
        s.rec <- s.rec[,rep(seq_len(ncol(s.rec)), each=2),drop=F]
        for(jj in seq(2, ncol(s.rec), 2)){
          s.rec[,jj-1] <- 1-s.rec[,jj]
        }
        p.prim <- p.prim[,rep(seq_len(ncol(p.prim)), each=2),drop=F]
        p.rec <- p.rec[,rep(seq_len(ncol(p.rec)), each=2), drop=F]
        
        ## also remove solutions where no mutation is present in either the primary or the relapse tumor
        no.mut.prim.rec <- rowSums(s.prim)==0 | rowSums(s.rec)==0
        s.prim <- s.prim[!no.mut.prim.rec,,drop=F]
        s.rec <- s.rec[!no.mut.prim.rec,,drop=F]
        p.prim <- p.prim[!no.mut.prim.rec,,drop=F]
        p.rec <- p.rec[!no.mut.prim.rec,,drop=F]
        p.prim.index <- p.prim.index[!no.mut.prim.rec,,drop=F]
        p.rec.index <- p.rec.index[!no.mut.prim.rec,,drop=F]
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.prim, s.rec, p.prim, p.rec))
        s.prim <- s.prim[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        s.rec <- s.rec[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        p.prim.index <- p.prim.index[!dupl.rows,,drop=F]
        p.rec.index <- p.rec.index[!dupl.rows,,drop=F]
        
      }else if(i <= length(common.genes)+length(genes.prim)){
        
        s.p.prim.rec <- Index.Matrix(tree = tree, copy.number.prim = copy.number.prim, copy.number.rec = copy.number.rec, which.prim = which.prim, which.rec = which.rec, mut.prim = T, mut.rec = F, locus = locus[i,],
                                     sex = sex, b.allele.prim = b.allele.prim, b.allele.rec = b.allele.rec, cnv.indicator = cnv.indicator, which.prim.cnv = which.prim.cnv, which.rec.cnv = which.rec.cnv)
        
        s.prim <- s.p.prim.rec$s.prim
        s.rec <- s.p.prim.rec$s.rec
        p.prim <- s.p.prim.rec$p.prim
        p.rec <- s.p.prim.rec$p.rec
        rm(s.p.prim.rec)
        
        
        if(!is.na(which.prim.cnv[1]) && copy.number.prim != ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)&
           b.allele.prim !=1){
          rows.to.keep.prim <- which(colSums(t((p.prim!=0)*1) == which.prim.cnv)==(length(which.prim)))# + add.prim.clone))
        }else{
          rows.to.keep.prim <- seq(1:nrow(p.prim))
        }
        p.prim <- p.prim[rows.to.keep.prim,,drop=F]
        p.rec <- p.rec[rows.to.keep.prim,,drop=F]
        s.prim <- s.prim[rows.to.keep.prim,,drop=F]
        s.rec <- s.rec[rows.to.keep.prim,,drop=F]
        
        ## add normal tissue (p.rec is important since CNVs can happen where no mutations are. s.rec is trivial since the relapse is not mutated)
        p.prim <- cbind(p.prim, 0)
        s.prim <- cbind(s.prim, 0)
        p.rec <- cbind(p.rec, 0)
        
        ## remove rows where no mutations are:
        rows.no.mut.prim <- which(rowSums(s.prim[,parms.prim!=0,drop=F])==0)
        if(length(rows.no.mut.prim)>0 ){
          s.prim <- s.prim[-rows.no.mut.prim,,drop=F]
          p.prim <- p.prim[-rows.no.mut.prim,,drop=F]
          p.rec <- p.rec[-rows.no.mut.prim,,drop=F]
        }
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.prim[,parms.prim!=0,drop=F], p.prim[,parms.prim!=0,drop=F]))
        s.prim <- s.prim[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        
        ## remove zero rows based on the clones that are occupied
        zero.rows <- rowSums(s.prim[,parms.prim!=0,drop=F])==0
        s.prim <- s.prim[!zero.rows,,drop=F]
        p.prim <- p.prim[!zero.rows,,drop=F]
        p.rec <- p.rec[!zero.rows,,drop=F]
        
        ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
        p.prim.index <- p.prim
        
        if(copy.number.prim==normal.copy.number && b.allele.prim == normal.baf){
          p.prim.index[,seq_len(ncol(p.prim))] <- 0
        }
        
        ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
        p.rec.index <- p.rec
        if(copy.number.rec == normal.copy.number && b.allele.rec == normal.baf){
          p.rec.index[,seq_len(ncol(p.rec))] <- 0
        }
        
        ## now, transform p.prim and p.rec into copy number table
        p.prim <- replace(p.prim, p.prim!=0, copy.number.prim)
        p.prim <- replace(p.prim, p.prim.index==0, normal.copy.number)
        if(copy.number.prim==0){
          p.prim <- replace(p.prim, p.prim.index[,seq_len(ncol(p.prim))]==1, 0)
        }
        ## add normal counterpart
        p.prim.index <- p.prim.index[,rep(seq_len(ncol(p.prim.index)),each=2),drop=F]
        
        ## now, transform p.prim and p.rec into copy number table
        p.rec <- replace(p.rec, p.rec!=0, copy.number.rec)
        p.rec <- replace(p.rec, p.rec.index==0, normal.copy.number)
        if(copy.number.rec==0){
          p.rec <- replace(p.rec, p.rec.index[,seq_len(ncol(p.rec))]==1, 0)
        }
        
        ## add normal counterpart
        p.rec.index <- p.rec.index[,rep(seq_len(ncol(p.rec.index)),each=2),drop=F]
        
        s.prim <- s.prim/p.prim
        s.prim <- replace(s.prim, is.infinite(s.prim),0)
        s.prim <- replace(s.prim, is.na(s.prim), 0)
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.prim[,parms.prim!=0,drop=F], p.prim[,parms.prim!=0,drop=F]))
        s.prim <- s.prim[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        p.prim.index <- p.prim.index[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        p.rec.index <- p.rec.index[!dupl.rows,,drop=F]
        
        ## now, add the normal counterparts for each subclone and divide by the respective copy number
        s.prim <- s.prim[,rep(seq_len(ncol(s.prim)), each=2),drop=F]
        for(jj in seq(2, ncol(s.prim), 2)){
          s.prim[,jj-1] <- 1-s.prim[,jj]
        }
        
        p.prim <- p.prim[,rep(seq_len(ncol(p.prim)), each=2),drop=F]
        
        ## do not evaluate anything of the relapse tumor in this case. Thus, we do not need s.rec
        s.rec <- matrix(0, nrow=nrow(s.prim), ncol=1)
        
      }else if(i <= length(common.genes)+length(genes.prim)+length(genes.rec)){
        
        ## adjust index to correct row
        i.rec <- i - length(genes.prim)
        
        
        s.p.prim.rec <- Index.Matrix(tree = tree, copy.number.prim = copy.number.prim, copy.number.rec = copy.number.rec, which.prim = which.prim, which.rec = which.rec, mut.prim = F, mut.rec = T, locus = locus[i,],
                                     sex = sex, b.allele.prim = b.allele.prim, b.allele.rec = b.allele.rec, cnv.indicator = cnv.indicator, which.prim.cnv = which.prim.cnv, which.rec.cnv = which.rec.cnv)
        s.prim <- s.p.prim.rec$s.prim
        s.rec <- s.p.prim.rec$s.rec
        p.prim <- s.p.prim.rec$p.prim
        p.rec <- s.p.prim.rec$p.rec
        
        
        
        if(!is.na(which.rec.cnv[1]) && copy.number.rec != ifelse(sex=="male" & locus[i,1] %in% c("X", "Y"), 1, 2)&
           b.allele.rec !=1){
          rows.to.keep.rec <- which(colSums(t((p.rec!=0)*1) == which.rec.cnv)==(length(which.rec)))# + add.rec.clone))
        }else{
          rows.to.keep.rec <- seq(1:nrow(p.rec))
        }
        
        p.prim <- p.prim[rows.to.keep.rec,,drop=F]
        p.rec <- p.rec[rows.to.keep.rec,,drop=F]
        s.prim <- s.prim[rows.to.keep.rec,,drop=F]
        s.rec <- s.rec[rows.to.keep.rec,,drop=F]
        
        ## add normal tissue (we don't need s.prim, since the mutation is not present in the primary)
        p.rec <- cbind(p.rec, 0)
        s.rec <- cbind(s.rec, 0)
        p.prim <- cbind(p.prim, 0)
        
        ## remove rows where no mutations are:
        
        rows.no.mut.rec <- which(rowSums(s.rec[,parms.rec!=0,drop=F])==0)
        if(length(rows.no.mut.rec)>0){
          s.rec <- s.rec[-rows.no.mut.rec,,drop=F]
          p.rec <- p.rec[-rows.no.mut.rec,,drop=F]
          p.prim <- p.prim[-rows.no.mut.rec,,drop=F]
        }
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.rec[,parms.rec!=0,drop=F], p.rec[,parms.rec!=0,drop=F]))
        s.rec <- s.rec[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        
        ## remove zero rows based on the clones that are occupied
        zero.rows <- rowSums(s.rec[,parms.rec!=0,drop=F])==0
        s.rec <- s.rec[!zero.rows,,drop=F]
        p.rec <- p.rec[!zero.rows,,drop=F]
        p.prim <- p.prim[!zero.rows,,drop=F]
        
        ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
        p.rec.index <- p.rec
        if(copy.number.rec==normal.copy.number && b.allele.rec == normal.baf){
          p.rec.index[,seq_len(ncol(p.rec.index))] <- 0
        }
        
        ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
        p.prim.index <- p.prim
        if(copy.number.prim==normal.copy.number && b.allele.prim == normal.baf){
          p.prim.index[,seq_len(ncol(p.prim))] <- 0
        }
        
        ## now, transform p.prim and p.rec into copy number table
        p.rec <- replace(p.rec, p.rec!=0, copy.number.rec)
        p.rec <- replace(p.rec, p.rec.index==0, normal.copy.number)
        if(copy.number.rec==0){
          p.rec <- replace(p.rec, p.rec.index[,seq_len(ncol(p.rec))]==1, 0)
        }
        ## add normal counterpart
        p.rec.index <- p.rec.index[,rep(seq_len(ncol(p.rec.index)),each=2),drop=F]
        
        p.prim <- replace(p.prim, p.prim!=0, copy.number.prim)
        p.prim <- replace(p.prim, p.prim.index==0, normal.copy.number)
        if(copy.number.prim==0){
          p.prim <- replace(p.prim, p.prim.index[,seq_len(ncol(p.prim))]==1, 0)
        }
        # add normal counterpart
        p.prim.index <- p.prim.index[,rep(seq_len(ncol(p.prim.index)),each=2),drop=F]
        
        s.rec <- s.rec/p.rec
        s.rec <- replace(s.rec, is.infinite(s.rec), 0)
        s.rec <- replace(s.rec, is.na(s.rec),0)
        
        ## remove duplicated rows based on the clones that are occupied
        dupl.rows <- duplicated(cbind(s.rec[,parms.rec!=0,drop=F], p.rec[,parms.rec!=0,drop=F]))
        s.rec <- s.rec[!dupl.rows,,drop=F]
        p.rec <- p.rec[!dupl.rows,,drop=F]
        p.rec.index <- p.rec.index[!dupl.rows,,drop=F]
        p.prim <- p.prim[!dupl.rows,,drop=F]
        p.prim.index <- p.prim.index[!dupl.rows,,drop=F]
        
        ## now, add the normal counterparts for each subclone and divide by the respective copy number
        s.rec <- s.rec[,rep(seq_len(ncol(s.rec)), each=2),drop=F]
        for(jj in seq(2, ncol(s.rec), 2)){
          s.rec[,jj-1] <- 1-s.rec[,jj]
        }
        
        p.rec <- p.rec[,rep(seq_len(ncol(p.rec)), each=2),drop=F]
        
        ## do not evaluate anything of the primary tumor in this case. Thus, we do not need s.prim
        s.prim <- matrix(0, nrow=nrow(s.rec), ncol=1)
        
        ## also remove solutions where no mutation is present in either the primary or the relapse tumor
        no.mut.prim.rec <- rowSums(s.rec)==0
        s.prim <- s.prim[!no.mut.prim.rec,,drop=F]
        s.rec <- s.rec[!no.mut.prim.rec,,drop=F]
        p.prim <- p.prim[!no.mut.prim.rec,,drop=F]
        p.rec <- p.rec[!no.mut.prim.rec,,drop=F]
        p.prim.index <- p.prim.index[!no.mut.prim.rec,,drop=F]
        p.rec.index <- p.rec.index[!no.mut.prim.rec,,drop=F]
        
      }
    }
    
    
    ## s is now a matrix that contains the (relative) copy numbers within each subclones for reference and mutated alleles. (already normalized by copy number)
    ## e.g. 0.5 0.5 0.0 1.0 0.0 1.0
    ##      0.0 1.0 0.5 0.5 0.0 1.0
    ##      0.5 0.5 0.5 0.5 0.0 1.0
    ## ..... (each row contains the number of reference reads in SCI, mutated reads in SCI, reference reads in SCII, mutated reads in SCII, ...)
    ## The last 2 columns correspond to normal tissue (always 0 in mutated reads, 1 in reference reads)
    
    ## for each gene iterate through each possible combination of s. Choose the most likely combination and store it in the result matrix
    ## store the likelihood
    LL <- 0
    ## store the most likely result table
    res.most.likely.prim <- 0
    res.most.likely.rec <- 0
    ## store the corresponding membership indices
    s.most.likely.prim <- 0
    s.most.likely.rec <- 0
    
    ## switch: if we are fitting, we don't do this, but later to get an idea about the overall likelihood, not only take into account the most likely solution, but the 90 % most likely ones
    if(likelihood.space){
      LL.all.sol.prim <- rep(0, nrow(s.prim))
      LL.all.sol.rec <- rep(0, nrow(s.rec))
    }
    
    if(length(read.count.matrix.rec)!=0){
      number.of.possibilities <- number.of.possibilities + log(nrow(s.prim))
    }else{
      number.of.possibilities <- number.of.possibilities + log(nrow(unique(s.prim)))
    }
    
    
    ## iterate through the possible combinations, stored in s
    for(k in seq_len(nrow(s.prim))){
      ## we only need to calculate the expected values if we have more than 1 clone. Otherwise the solution is trivial (see below). We can now optimize separately for primary and relapse
      ## start with the primary. Only do anything if there are primary genes mutated at this iteration
      if(i <= (length(common.genes)+length(genes.prim))){
        
        if(sum(s.prim[k,]>0)>1 && sum(parms.prim>0)>1){
          ## the expected value of the reference read counts per subclone is calculated as p*r_ref_ges = mu_I*(1-s_I/copy.number_I)*copy.number_I/(sum(mu_i(1-s_i/copy.number_i)*copy.number_i))*r_ref_ges
          ## calculate separately for mutated and ref
          if(sum(s.prim[k,seq(1,ncol(s.prim),2)])==0){ 
            res.mle.ref <- rep(0, ncol(s.prim)/2)
          }else if(sum(parms.prim*s.prim[k,seq(1,ncol(s.prim),2)])==0){
            res.mle.ref <- rep(0, ncol(s.prim)/2)
          }else{
            res.mle.ref <- parms.prim*(s.prim[k,seq(1,ncol(s.prim),2)] * p.prim[k, seq(1, ncol(p.prim), 2)]/
                                         sum(parms.prim*s.prim[k, seq(1, ncol(s.prim), 2)] * p.prim[k, seq(1, ncol(p.prim), 2)])) * 
              read.count.matrix.prim[i,1]
          }
          if(sum(parms.prim*s.prim[k,seq(2,ncol(s.prim),2)])==0){
            res.mle.mut <- rep(0, ncol(s.prim)/2)
          }else{
            res.mle.mut <- parms.prim * (s.prim[k,seq(2,ncol(s.prim),2)] * p.prim[k, seq(2, ncol(p.prim), 2)]/
                                           sum(parms.prim * s.prim[k, seq(2, ncol(s.prim), 2)] * p.prim[k, seq(2, ncol(p.prim), 2)])) *
              read.count.matrix.prim[i,2]
          }
          ## combine mutated and ref, sort, such that ref and mutated intercalate (ref clone I, mut clone I ref clone II, mut clone II, ... as in s-matrix etc.)
          res.mle.prim <- cbind(res.mle.ref, res.mle.mut)
          res.mle.prim <- as.vector(t(res.mle.prim))
          res.mle.prim <- round(res.mle.prim)
          
          
        }else{
          ## trivial solution if we only have 1 clone
          res.mle.mut <- 0
          res.mle.prim <- round(c(read.count.matrix.prim[i,], rep(0, ncol(s.prim) -2)))
        }
        ## corresponding negative log likelihood (based on a multinomial. Here, we are only interested in how likely the subclonal structure is, not the mutated/ref structure within the subclones)
        if(sum(res.mle.mut)==0){ ## set -logL to infinity if no mutated reads are expected (since we have them!)
          LL.new.prim <- Inf
        }else if(sum(parms.prim*s.prim[k,seq(1,ncol(s.prim),2)])==0){
          LL.new.prim <- Inf
        }else{
          LL.new.prim <- NegLogLikelihood(readcounts = matrix(res.mle.prim, nrow=1), 
                                          copy.number = s.prim[k, ] * matrix(p.prim[k, ] / sum(parms.prim * p.prim[k,seq(1, ncol(p.prim), 2)]), nrow=1),
                                          parms = rep(parms.prim, each=2))
        }
        
        if(likelihood.space){
          LL.all.sol.prim[k] <- LL.new.prim
        }
      }else{LL.new.prim <- 0
      if(likelihood.space){
        LL.all.sol.prim[k] <- LL.new.prim
      }
      res.mle.prim <- c()
      }
      ## now the same for the recurrence. Only do anything if there are primary genes mutated at this iteration
      if(i <= (length(common.genes)) || (i > (length(common.genes)+length(genes.prim)) && i<=(length(common.genes)+length(genes.prim)+length(genes.rec)))){
        if(sum(s.rec[k,]>0)>1 && sum(parms.rec>0)>1){
          ## the expected value of the reference read counts per subclone is calculated as p*r_ref_ges = mu_I*(1-s_I/copy.number_I)*copy.number_I/(sum(mu_i(1-s_i/copy.number_i)*copy.number_i))*r_ref_ges
          ## calculate separately for mutated and ref
          if(sum(s.rec[k,seq(1,ncol(s.rec),2)])==0 || sum(parms.rec*s.rec[k,seq(1,ncol(s.rec),2)])==0){
            res.mle.ref <- rep(0, ncol(s.rec)/2)
          }else{
            res.mle.ref <- parms.rec*(s.rec[k,seq(1,ncol(s.rec),2)] * p.rec[k, seq(1, ncol(p.rec), 2)]/
                                        sum(parms.rec*s.rec[k, seq(1, ncol(s.rec), 2)] * p.rec[k, seq(1, ncol(p.rec), 2)])) * 
              read.count.matrix.rec[i.rec,1]
          }
          if(sum(parms.rec*s.rec[k,seq(2,ncol(s.rec),2)])==0){
            res.mle.mut <- rep(0, ncol(s.rec)/2)
          }else{
            res.mle.mut <- parms.rec*(s.rec[k,seq(2,ncol(s.rec),2)] * p.rec[k, seq(2, ncol(p.rec), 2)]/
                                        sum(parms.rec*s.rec[k, seq(2, ncol(s.rec), 2)] * p.rec[k, seq(2, ncol(p.rec), 2)])) * 
              read.count.matrix.rec[i.rec,2]
          }
          ## combine mutated and ref, sort, such that ref and mutated intercalate (ref clone I, mut clone I ref clone II, mut clone II, ... as in s-matrix etc.)
          res.mle.rec <- cbind(res.mle.ref, res.mle.mut)
          res.mle.rec <- as.vector(t(res.mle.rec))
          res.mle.rec <- round(res.mle.rec)
        }else{
          ## trivial solution if we only have 1 clone
          res.mle.mut <- 0
          res.mle.rec <- round(c(read.count.matrix.rec[i.rec,], rep(0, ncol(s.rec) -2)))
        }
        ## corresponding negative log likelihood (based on a multinomial. Here, we are only interested in how likely the subclonal structure is, not the mutated/ref structure within the subclones)
        if(sum(res.mle.mut)==0 || sum(parms.rec*s.rec[k,seq(1,ncol(s.rec),2)])==0){ ## set -logL to infinity if no mutated reads are expected (since we have them!)
          LL.new.rec <- Inf
        }else{
          LL.new.rec <- NegLogLikelihood(readcounts = matrix(res.mle.rec, nrow=1), 
                                         copy.number = s.rec[k, ] * matrix(p.rec[k, ] / sum(parms.rec * p.rec[k,seq(1, ncol(p.rec), 2)]), nrow=1),
                                         parms = rep(parms.rec, each=2))
        }
        if(likelihood.space){
          LL.all.sol.rec[k] <- LL.new.rec
        }
      }else{LL.new.rec <- 0
      if(likelihood.space){
        LL.all.sol.rec[k] <- LL.new.rec
      }
      res.mle.rec <- c()}
      
      ## initialize LL in the first iteration, res.most.likely, s.most.likely
      if(k==1){
        LL <- LL.new.prim + LL.new.rec
        s.most.likely.prim <- s.prim[k,]
        p.most.likely.prim <- p.prim[k,]
        p.most.likely.prim.index <- p.prim.index[k,]
        res.most.likely.prim <- res.mle.prim
        s.most.likely.rec <- s.rec[k,]
        p.most.likely.rec <- p.rec[k,]
        p.most.likely.rec.index <- p.rec.index[k,]
        res.most.likely.rec <- res.mle.rec
        
      }
      
      ## if the current result is more likely than any previous ones, replace res.most.likely with this result
      if((LL.new.prim + LL.new.rec )<LL & !any(c(res.mle.prim, res.mle.rec)<0)){
        res.most.likely.prim <- res.mle.prim
        res.most.likely.rec <- res.mle.rec
        s.most.likely.prim <- s.prim[k,]
        s.most.likely.rec <- s.rec[k,]
        p.most.likely.prim <- p.prim[k,]
        p.most.likely.rec <- p.rec[k,]
        p.most.likely.prim.index <- p.prim.index[k,]
        p.most.likely.rec.index <- p.rec.index[k,]
        LL <- LL.new.prim + LL.new.rec
      }
    }
    
    ## store the most likely result for this gene in the result table that summarizes the most likely results for all genes
    if(i <= (length(common.genes)+length(genes.prim))){
      res.prim[i,] <- res.most.likely.prim
      s.res.prim[i,] <- s.most.likely.prim
      p.res.prim[i,] <- p.most.likely.prim
      p.res.prim.index[i,] <- p.most.likely.prim.index
      if(i > length(common.genes) && length(genes.rec)>0){
        p.res.rec.index[i,] <- p.most.likely.rec.index
      }
      ## if we are interested in a likelihood space (all solutions within the 95 % LL), then add additional solutions and store their likelihood
      if(likelihood.space){
        
        ## Residuals between measured and expected variant allele frequencies
        residuals.vafs.prim[i] <- (sum(parms.prim*s.most.likely.prim[seq(2,length(s.most.likely.prim),2)]*p.most.likely.prim[seq(2,length(s.most.likely.prim),2)])/sum(parms.prim*p.most.likely.prim[seq(2,length(s.most.likely.prim),2)])*sum(read.count.matrix.prim[i,])-
                                     read.count.matrix.prim[i,2])/sum(read.count.matrix.prim[i,])
        
        
        like.per.loc.prim[i] <- exp(-LL)/(sum(exp(-LL.all.sol.prim - LL.all.sol.rec)))
        ## if the total likelihood is too small we cannot estimate a likelihood per solution
        if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
          like.per.loc.prim[i] <- NA
        }
        
        if(i > length(common.genes)){
          if(length(read.count.matrix.rec)>0){
            like.per.loc.rec[i] <- exp(-LL)/(sum(exp(-LL.all.sol.prim -LL.all.sol.rec)))
          }else{
            like.per.loc.rec[i] <-0
          }
          
          if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
            like.per.loc.rec[i] <- NA
          }
        }
        ## accept additional solutions until 90 % of the total likelihood are explained
        additional.solutions <- min(which(log(cumsum(sort(exp(-LL.all.sol.prim-LL.all.sol.rec), decreasing=T))) - log(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))) > log(0.9)))
        
        if(length(additional.solutions)> 0){
          if(min(additional.solutions)==1 | is.na(additional.solutions) | additional.solutions==Inf){
            additional.solutions <- c()
          }
          
        }
        additional.solutions <- order(exp(-LL.all.sol.prim - LL.all.sol.rec), decreasing=T)[additional.solutions]
        
        if(length(additional.solutions)>0){
          like.per.loc.prim <- c(like.per.loc.prim, exp((-LL.all.sol.prim - LL.all.sol.rec)[additional.solutions] - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec)))))
          if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
            like.per.loc.prim[length(like.per.loc.prim)] <- NA
          }
          s.res.prim <- rbind(s.res.prim, s.prim[additional.solutions,])
          
          rownames(s.res.prim) <- c(rownames(s.res.prim)[-((nrow(s.res.prim)-length(additional.solutions)+1):nrow(s.res.prim))], rep(rownames(s.res.prim)[i], length(additional.solutions)))
          p.res.prim <- rbind(p.res.prim, p.prim[additional.solutions,])
          rownames(p.res.prim) <- rownames(s.res.prim)
          
          p.res.prim.index <- rbind(p.res.prim.index, p.prim.index[additional.solutions,])
          rownames(p.res.prim.index) <- c(rownames(p.res.prim.index)[-((nrow(p.res.prim.index)-length(additional.solutions)+1):nrow(p.res.prim.index))], rep(rownames(p.res.prim.index)[i], length(additional.solutions)))
          if(i > length(common.genes)){
            if(length(read.count.matrix.rec)>0){
              p.res.rec.index <- rbind(p.res.rec.index, p.rec.index[additional.solutions,])
              rownames(p.res.rec.index) <- rownames(p.res.prim.index)
              like.per.loc.rec <- c(like.per.loc.rec, exp((-LL.all.sol.prim - LL.all.sol.rec)[additional.solutions] - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec)))))
              if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
                like.per.loc.rec[length(like.per.loc.rec)] <- NA
              }   
            }
            
          }
          
        }
      }
      
    }
    
    if(i <= (length(common.genes)) || (i > (length(common.genes)+length(genes.prim)) && i<=(length(common.genes)+length(genes.prim)+length(genes.rec)))){
      res.rec[i.rec,] <- res.most.likely.rec
      s.res.rec[i.rec,] <- s.most.likely.rec
      p.res.rec[i.rec,] <- p.most.likely.rec
      p.res.rec.index[i,] <- p.most.likely.rec.index
      if(i > length(common.genes)){
        p.res.prim.index[i,] <- p.most.likely.prim.index
      }
      
      
      ## if we are interested in a likelihood space (all solutions within the 90 % LL), then add additional solutions and store their likelihood
      if(likelihood.space){
        
        residuals.vafs.rec[i] <- (sum(parms.rec*s.most.likely.rec[seq(2,length(s.most.likely.rec),2)]*p.most.likely.rec[seq(2,length(p.most.likely.rec),2)])/sum(parms.rec*p.most.likely.rec[seq(2,length(s.most.likely.rec),2)])*sum(read.count.matrix.rec[i.rec,]) -
                                    read.count.matrix.rec[i.rec,2])/sum(read.count.matrix.rec[i.rec,])
        
        like.per.loc.rec[i] <- exp(-LL - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec))))
        if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
          like.per.loc.rec[i] <- NA
        }    
        if(i > length(common.genes)){
          like.per.loc.prim[i] <- exp(-LL - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec))))
          if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
            like.per.loc.prim[i] <- NA
          }    
        }
        additional.solutions <- min(which(log(cumsum(sort(exp(-LL.all.sol.prim-LL.all.sol.rec), decreasing=T))) - log(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))) > log(0.9)))
        
        
        if(length(additional.solutions)> 0){
          if(min(additional.solutions)  == 1 | is.na(additional.solutions)| additional.solutions==Inf){
            additional.solutions <- c()
          }
        }
        
        additional.solutions <- order(-LL.all.sol.prim - LL.all.sol.rec, decreasing=T)[additional.solutions]
        
        if(length(additional.solutions)>0){
          
          like.per.loc.rec <- c(like.per.loc.rec, exp((-LL.all.sol.prim - LL.all.sol.rec)[additional.solutions] - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec)))))
          if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
            like.per.loc.rec[length(like.per.loc.rec)] <- NA
          }    
          s.res.rec <- rbind(s.res.rec, s.rec[additional.solutions,])
          rownames(s.res.rec) <- c(rownames(s.res.rec)[-((nrow(s.res.rec)-length(additional.solutions)+1):nrow(s.res.rec))], rep(rownames(s.res.rec)[i.rec], length(additional.solutions)))
          p.res.rec <- rbind(p.res.rec, p.rec[additional.solutions,])
          rownames(p.res.rec) <- rownames(s.res.rec)
          p.res.rec.index <- rbind(p.res.rec.index, p.rec.index[additional.solutions,])
          rownames(p.res.rec.index) <- c(rownames(p.res.rec.index)[-((nrow(p.res.rec.index)-length(additional.solutions)+1):nrow(p.res.rec.index))], rep(rownames(p.res.rec.index)[i], length(additional.solutions)))
          if(i > length(common.genes)){
            p.res.prim.index <- rbind(p.res.prim.index, p.prim.index[additional.solutions,])
            rownames(p.res.prim.index) <- rownames(p.res.rec.index)
            like.per.loc.prim <- c(like.per.loc.prim, exp((-LL.all.sol.prim - LL.all.sol.rec)[additional.solutions] - log(sum(exp(-LL.all.sol.prim -LL.all.sol.rec)))))
            if(sum(exp(-LL.all.sol.prim - LL.all.sol.rec))==0){
              like.per.loc.prim[length(like.per.loc.prim)] <- NA
            }    
          }
        }
      }
      
    }
    
    
    
    
    
  }
  
  
  res.list <- list(readcounts.prim=res.prim, rel.mutated.alleles.prim=s.res.prim, copy.number.prim = p.res.prim, readcounts.rec=res.rec, rel.mutated.alleles.rec=s.res.rec, copy.number.rec = p.res.rec,
                   copy.number.index.prim = p.res.prim.index, copy.number.index.rec = p.res.rec.index, number.of.possibilities = number.of.possibilities)
  
  if(likelihood.space){
    names(like.per.loc.prim) <- rownames(p.res.prim.index)
    names(like.per.loc.rec) <- rownames(p.res.rec.index)
    res.list <- list(readcounts.prim=res.prim, rel.mutated.alleles.prim=s.res.prim, copy.number.prim = p.res.prim, readcounts.rec=res.rec, rel.mutated.alleles.rec=s.res.rec, copy.number.rec = p.res.rec,
                     copy.number.index.prim = p.res.prim.index, copy.number.index.rec = p.res.rec.index, like.per.loc.prim = like.per.loc.prim, like.per.loc.rec = like.per.loc.rec,
                     number.of.possibilities = number.of.possibilities, cnv.indicator = cnv.indicator.all, additional.copy.number.prim = additional.copy.number.prim, additional.copy.number.rec = additional.copy.number.rec,
                     residuals.coverage.ratio.prim = residuals.coverage.ratio.prim, residuals.coverage.ratio.rec = residuals.coverage.ratio.rec,
                     residuals.vafs.prim = residuals.vafs.prim, residuals.vafs.rec = residuals.vafs.rec, residuals.bafs.prim = residuals.bafs.prim,
                     residuals.bafs.rec = residuals.bafs.rec)
  }
  
  return(res.list)
}





#' Calculate the maximum likelihood estimate for the relative subclone sizes at given expected readcounts/copy numbers
#' 
#' Maximizes the log-likelihood with respect to the subclone sizes. 
#' @param readcounts a readcount table as returned by Expectation
#' @param copy.number the copy numbers as returened by Expectation
#' @param n.clones the number of subclones (tumor plus normal tissue)
#' @return The maximum likelihood estimate of the relative subclone sizes
#' @export

Maximization <- function(readcounts, copy.number, n.clones){
  ## to maximize the logL we introduce a lagrange multiplier lambda) to account for the constraint that all probabilities add up to 1.
  parms <- rep(0, n.clones)
  parms <- colSums((readcounts[,seq(1, n.clones*2, 2), drop=F] + readcounts[,seq(2, n.clones*2, 2), drop=F])/replace(copy.number[,seq(1, n.clones*2,2)], copy.number[,seq(1, n.clones*2,2)] ==0, 1))/
    sum(readcounts/replace(copy.number, copy.number==0, 1))
  parms <- replace(parms, is.na(parms), 0)
  return(parms)
}





#' EM algorithm for phylogenetic inference
#' 
#' This function iteratively repeats the expectation and maximization steps until convergence is reached (based on a cutoff value). 
#' @param tree A binary tree matrix as returned by Tree.matrix.R
#' @param read.count.matrix.prim a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the primary tumor
#' @param read.count.matrix.rec a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the relapse tumor
#' @param coverage.ratio.baf.prim a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the primary tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param coverage.ratio.baf.rec a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the relapse tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param which.prim vector indicating which subclones are primary subclones
#' @param which.rec vector indicating which subclones are relapse subclones
#' @param locus vector containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param parms.prim vector containing the relative subclone sizes of the primary tumor
#' @param parms.rec vector containing the relative subclone sizes of the relapse tumor
#' @param n.clones.prim number of primary subclones (excluding normal tissue)
#' @param n.clones.rec number of secondary subclones (excluding normal tissue)
#' @param cutoff cutoff for EM algorithm (RSS between two subsequent iterations)
#' @param maxit maximal number of iterations
#' @return The relative subclone sizes of primary and relapse tumor (last value corresponds to normal tissue).
#' @export

EM.alg <- function(tree, read.count.matrix.prim, read.count.matrix.rec, parms.prim, parms.rec, coverage.ratio.baf.prim, coverage.ratio.baf.rec, locus, 
                   n.clones.prim, n.clones.rec, cutoff, maxit = 25, sex = sex, which.prim, which.rec){
  
  ## initialization
  ## start with expectation
  prev.it <- Expectation(tree = tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = parms.prim, parms.rec = parms.rec,
                         coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, 
                         sex = sex, which.prim = which.prim, which.rec = which.rec)
  
  
  if(prev.it[[1]][1]=="PT_impurity"){return(NA)}
  if(prev.it[[1]][1]=="RT_impurity"){return(NA)}
  
  
  ## corresponding negative log likelihood
  
  prev.ll <- NegLogLikelihood(readcounts=prev.it$readcounts.prim, copy.number=prev.it$rel.mutated.alleles.prim*prev.it$copy.number.prim/
                                colSums(parms.prim*t(prev.it$copy.number.prim[,seq(1, n.clones.prim*2+2, 2), drop=F])), parms = rep(parms.prim, each=2)) +
    NegLogLikelihood(readcounts=prev.it$readcounts.rec, copy.number=prev.it$rel.mutated.alleles.rec*prev.it$copy.number.rec/
                       colSums(parms.rec*t(prev.it$copy.number.rec[,seq(1, n.clones.rec*2+2, 2), drop=F])), 
                     parms = rep(parms.rec, each=2))
  
  ## maximization
  prev.it.prim <- Maximization(readcounts = prev.it$readcounts.prim, copy.number = prev.it$copy.number.prim, n.clones = n.clones.prim + 1)
  prev.it.rec <- Maximization(readcounts = prev.it$readcounts.rec, copy.number = prev.it$copy.number.rec, n.clones = n.clones.rec + 1)
  
  rm(prev.it)
  gc()
  it.prim <- prev.it.prim
  it.rec <- prev.it.rec
  
  ## iterate
  it.counter <- 0
  while(it.counter < maxit){
    
    prev.it.prim[prev.it.prim < 0.01] <- 0
    prev.it.prim <- prev.it.prim/sum(prev.it.prim)
    prev.it.rec[prev.it.rec <0.01] <- 0
    prev.it.rec <- prev.it.rec/sum(prev.it.rec)
    
    it.clus <- Expectation(tree = tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = prev.it.prim, parms.rec = prev.it.rec,
                           coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, 
                           sex = sex, which.prim = which.prim, which.rec = which.rec)
    
    if(it.clus[[1]][1]=="PT_impurity"){return(NA)}
    if(it.clus[[1]][1]=="RT_impurity"){return(NA)}
    
    it.prim <- Maximization(readcounts = it.clus$readcounts.prim, copy.number = it.clus$copy.number.prim, n.clones = n.clones.prim + 1)
    it.rec <- Maximization(readcounts = it.clus$readcounts.rec, copy.number = it.clus$copy.number.rec, n.clones = n.clones.rec + 1)
    if(sum(it.prim[-length(it.prim)]==0)==(length(it.prim)-1)){return(NA)}
    if(sum(it.rec[-length(it.rec)]==0)==(length(it.rec)-1)){return(NA)}
    
    ll <- NegLogLikelihood(readcounts=it.clus$readcounts.prim, copy.number=it.clus$rel.mutated.alleles.prim*it.clus$copy.number.prim/
                             colSums(it.prim*t(it.clus$copy.number.prim[,seq(1, n.clones.prim*2+2, 2), drop=F])),
                           parms = rep(it.prim, each=2))  +
      NegLogLikelihood(readcounts=it.clus$readcounts.rec, 
                       copy.number=it.clus$rel.mutated.alleles.rec*it.clus$copy.number.rec/
                         colSums(it.rec*t(it.clus$copy.number.rec[, seq(1, n.clones.rec*2+2, 2), drop=F])), parms = rep(it.rec, each=2))
    rm(it.clus)
    gc()
    
    ## comparison between parameter sets is only possible for equally sized data frames. If we remove a cluster due to a mean of zero, we cannot compare the mean
    ## values to the former parameter set and thus, allow one more iteration before we do a comparison again
    
    if(length(c(it.prim, it.rec))==length(c(prev.it.prim, prev.it.rec))){
      if(sum((c(it.prim-prev.it.prim, it.rec - prev.it.rec))^2)<cutoff){
        break}
      
      if(ll > prev.ll){
        
        it.prim <- prev.it.prim
        it.rec <- prev.it.rec
        break}
    }
    
    prev.it.prim <- it.prim
    prev.it.rec <- it.rec
    
    prev.ll <- ll
    
    it.counter <- it.counter + 1
  }
  return(list(prim=it.prim, rec=it.rec))
  
}




#' Repeats the EM algorithm often with random starting parameters to find the global optimum (EM itself does only local optimization)
#' 
#' All possible trees for up to three primary and relapse subclones are automatically tested. More complex trees can be tested as well, but information must be provided on how these trees should look like
#' @param read.count.matrix.prim a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the primary tumor
#' @param read.count.matrix.rec a matrix storing the number of reference reads (column 1) and the number of mutated reads (column 2) at mutated loci in the relapse tumor
#' @param coverage.ratio.baf.prim a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the primary tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param coverage.ratio.baf.rec a matrix storing the measured coverage number ratio (column 1) and B allele frequency (column2) of the relapse tumor (these values should be reported at all loci that are mutated in either of the tumor samples)
#' @param locus vector containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param cutoff the cutoff for the EM algorithm
#' @param it the number of iterations per candidate tree (random starting conditions)
#' @param directory directory for output file (not saving directory)
#' @param tree.best.fit if we already have a good tree and want to know whether an additional subclone improves the fit, provide a tree matrix here. Else, leave as NA
#' @param which.prim.best.fit if tree.best.fit is provided, provide a vector indicating which subclones are primary subclones here
#' @param which.rec.best.fit if tree.best.fit is provided, provide a vector indicating which subclones are relapse subclones here
#' @param add.prim.clone if tree.best.fit is provided, indicate here, whether a primary subclone should be added (add.rec.clone must be FALSE in this case)
#' @param add.rec.clone if tree.best.fit is provided, indicate here, whether a relapse subclone should be added (add.prim.clone must be FALSE in this case)
#' @return A list containing for each candidate tree matrices which contain the relative subclone sizes for the primary and the relapse tumor for each fit (rows correspond to fits) and a vector containing the corresponding likelihood. In the case of more complex trees, an additional list element reports the candidate tree and vectors indicating which subclones are primary clones and which are relapse clones.
#' @export


model.fits <- function(read.count.matrix.prim, read.count.matrix.rec, coverage.ratio.baf.prim, coverage.ratio.baf.rec, locus,
                       cutoff=5e-04, it=100, sex = "female", directory = "",
                       tree.best.fit = NA, which.prim.best.fit = NA, which.rec.best.fit = NA, add.prim.clone = F, add.rec.clone = F,
                       clusterComputation = F, n.cores = 8){
  resultlist <- list(0)
  
  if(clusterComputation){
    library(doParallel)
    cl<-makeCluster(n.cores, type = "FORK", outfile=paste0(directory, "/Log.txt"))
    registerDoParallel(cores=n.cores)
  }
  
  
  ## test all combinations of primary and relapse subclones for up to three subclones each (trees defined in function Tree.matrix):
  
  n.clones.prim.list <- c(1, 2, 1, 2, 2, 3, 3, 1, 1, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3)
  n.clones.rec.list <-  c(1, 1, 2, 2, 2, 1, 1, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3)
  trees <- c("classic", "classic", "classic", "classic", "asymmetric", "classic", "asymmetric",  "classic", "asymmetric", "classic", "asymmetric", "classic",
             "classic", "asymmetric", "classic", "symmetric", "classic", "asymmetric", "asymmetric2", "asymmetric3", "asymmetric4")
  which.prim <- list(list(c(1)), 
                     list(c(3,5), c(3,4)), 
                     list(c(5), c(4)), 
                     list(c(1, 2), c(1, 3)), 
                     list(c(1,3), c(1,5), c(3,6), c(5,6)), 
                     list(c(1,2,3)), 
                     list(c(1,3,5), c(1,5,6), c(3,5,6)),
                     list(c(5)), 
                     list(c(6), c(3), c(1)), 
                     list(c(1,2,3), c(1,3,5), c(1,5,6), c(3,5,6), c(1,2,5)), 
                     list(c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4), c(2,4,5), c(1,4,5), c(3,4,5)),
                     list(c(1,3,4), c(1,4,5), c(3,4,5)), 
                     list(c(5,6), c(2,6), c(2,3), c(1,2), c(3,6)), 
                     list(c(4,5), c(3,5), c(2,5), c(1,5), c(1,3), c(2,3), c(1,2)), 
                     list(c(5,6), c(3,6), c(1,6)), 
                     list(c(1,2,3), c(2,3,4), c(2,3,5), c(1,3,4), c(1,3,5)), 
                     list(c(1,3,5), c(1,2,4), c(3,5,6), c(1,3,4)), 
                     list(c(1,2,3), c(1,2,4), c(1,2,5), c(1,3,4),
                          c(1,3,5), c(1,4,5), c(2,3,4), c(2,4,5), 
                          c(2,5,6), c(3,4,5), c(3,5,6), c(4,5,6)),
                     list(c(1,2,3), c(1,2,4), c(1,2,5), c(1,3,4), c(1,3,5), c(1,4,5), c(1,5,6), 
                          c(4,5,6), c(3, 5, 6), c(3,4,6)),
                     list(c(1,2,3), c(1,3,4), c(1,3,5), c(2,3,4), c(2,3,5), c(3,4,5)),
                     list(c(1,2,3), c(1,2,4), c(1,2,5), c(1,4,5), c(1,5,6), c(2,3,4), c(2,3,5), c(3,4,5), c(3,5,6), c(4,5,6)))
  which.rec <-  list(list(c(2)), 
                     list(c(4), c(5)), 
                     list(c(3,4), c(3,5)),
                     list(c(3, 5), c(2, 5)), 
                     list(c(5,6), c(3,6), c(1,5), c(1,3)), 
                     list(c(5)), 
                     list(c(6), c(3), c(1)),
                     list(c(1,2,3)), 
                     list(c(1,3,5), c(1,5,6), c(3,5,6)), 
                     list(c(5,6), c(2,6), c(2,3), c(1,2), c(3,6)), 
                     list(c(4,5), c(3,5), c(2,5), c(1,5), c(1,3), c(2,3), c(1,2)), 
                     list(c(5,6), c(3,6), c(1,6)), list(c(1,2,3), c(1,3,5), c(1,5,6), c(3,5,6), c(1,2,5)), 
                     list(c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4), c(2,4,5), c(1,4,5), c(3,4,5)),
                     list(c(1,3,4), c(1,4,5), c(3,4,5)), 
                     list(c(4,5,6), c(1,5,6), c(1,4,6), c(2,5,6), c(2,4,6)), 
                     list(c(2,4,6), c(3,5,6), c(1,2,4), c(2,5,6)), 
                     list(c(4,5,6), c(3,5,6), c(3,4,6), c(2,5,6),
                          c(2,4,6), c(2,3,6), c(1,5,6), c(1,3,6),
                          c(1,3,4), c(1,2,6), c(1,2,4), c(1,2,3)),
                     list(c(4,5,6), c(3,5,6), c(3,4,6), c(2,5,6), c(2,4,6), c(2,3,6), c(2,3,4), 
                          c(1,2,3), c(1,2,4), c(1,2,5)),
                     list(c(4,5,6), c(2,5,6), c(2,4,6), c(1,5,6), c(1,4,6), c(1,2,6)),
                     list(c(4,5,6), c(3,5,6), c(3,4,6), c(2,3,6), c(2,3,4), c(1,5,6), c(1,4,6), c(1,2,6), c(1,2,4), c(1,2,3)))
  
  ## if we already tested up to 3 subclones and now want to add another one, we run into this part of the program:
  
  if(add.prim.clone | add.rec.clone){
    n.clones.prim <- length(which.prim.best.fit) + ifelse(add.prim.clone, 1, 0)
    n.clones.rec <- length(which.rec.best.fit) + ifelse(add.rec.clone, 1, 0)
    
    j <- 1
    while(tree.matrix(tree="classic", add.prim.clone = add.prim.clone, add.rec.clone = add.rec.clone, which.clone.add = j, 
                      which.prim = which.prim.best.fit, which.rec = which.rec.best.fit, best.tree = tree.best.fit)!="stop"){
      
      print(j)
      
      sim.res.prim <- c()
      sim.res.rec <- c()
      sim.like <- c()
      print(paste(n.clones.prim," primary,", n.clones.rec, " relapse"))
      
      tree.for.fit <-  tree.matrix(tree="classic", add.prim.clone = add.prim.clone, add.rec.clone = add.rec.clone, which.clone.add = j, 
                                   which.prim = which.prim.best.fit, which.rec = which.rec.best.fit, best.tree = tree.best.fit)
      
      fits <- foreach(i=1:it, .combine=rbind) %dopar%{
        set.seed(i)
        
        parms.prim <- Rand.parms(n.clones.prim + 1)
        parms.rec <- Rand.parms(n.clones.rec + 1)
        
        res <- EM.alg(tree = tree.for.fit$tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = parms.prim, parms.rec = parms.rec,
                      coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, n.clones.prim = n.clones.prim, 
                      n.clones.rec = n.clones.rec, cutoff, sex = sex, maxit = 50, which.prim = tree.for.fit$which.prim, which.rec = tree.for.fit$which.rec)
        
        if(!is.na(res[[1]])){
          if(res[[1]]!="PT_impurity"){
            if(res[[1]]!="RT_impurity"){
              
              res. <- Expectation(tree = tree.for.fit$tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = res$prim, 
                                  parms.rec = res$rec, coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, 
                                  sex = sex, which.prim = tree.for.fit$which.prim, which.rec = tree.for.fit$which.rec)
              
              
              if(!any(c(res$prim, res$rec)<0) & !any(res.$readcounts.prim<0) & !any(res.$readcounts.rec<0)){
                
                
                like <- NegLogLikelihood(readcounts =  res.$readcounts.prim, 
                                         copy.number = res.$rel.mutated.alleles.prim*res.$copy.number.prim/
                                           colSums(res$prim*t(res.$copy.number.prim[,seq(1, n.clones.prim*2+2,2),drop=F])),
                                         parms = rep(res$prim, each=2)) +
                  NegLogLikelihood(readcounts =  res.$readcounts.rec, 
                                   copy.number = res.$rel.mutated.alleles.rec*res.$copy.number.rec/
                                     colSums(res$rec*t(res.$copy.number.rec[,seq(1, n.clones.rec*2+2,2),drop=F])),
                                   parms = rep(res$rec, each=2))
                
                
                return.value <- c(res$prim, res$rec, like)
                
              }else{
                return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
              }
              rm(res.)
              gc()
            }else{
              return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
            }
          }else{
            return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
          }
        }else{
          return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
        }
        
        return.value
        
      }
      
      if(length(fits)>0){
        for(i.fits in 1:nrow(fits)){
          if(!is.infinite(fits[i.fits, ncol(fits)])){
            sim.like <- c(sim.like, fits[i.fits, ncol(fits)])
            sim.res.prim <- rbind(sim.res.prim, fits[i.fits, 1:(n.clones.prim+1)])
            sim.res.rec <- rbind(sim.res.rec, fits[i.fits, (n.clones.prim+2):(ncol(fits)-1)])
          }
        }
      }
      
      tree. <- tree.for.fit$tree
      
      
      if(add.prim.clone){
        new.best.tree <- tree.[,order(c(which.prim.best.fit, ncol(tree.), which.rec.best.fit))]
        new.which.prim <- c(which.prim.best.fit, n.clones.prim + n.clones.rec)
      }else{
        new.which.prim <- which.prim.best.fit
      }
      
      if(add.rec.clone){
        new.best.tree <- tree.[,order(c(which.prim.best.fit, which.rec.best.fit, ncol(tree.)))]
        new.which.rec <- c(which.rec.best.fit, n.clones.prim + n.clones.rec)
      }else{
        new.which.rec <- which.rec.best.fit
      }      
      
      resultlist[[j]] <- list(Parameters.prim=sim.res.prim, Parameters.rec=sim.res.rec, Likelihood=sim.like, tree = new.best.tree,
                              which.prim = new.which.prim, which.rec = new.which.rec)
      
      j <- j+1
    }
    
    if(clusterComputation){
      stopCluster(cl)
    }
    
    return(resultlist)
  }
  
  ii.start <- 1
  if(file.exists(paste0(directory,"/resultlist.RData")) & directory != "./test"){
    load(paste0(directory,"/resultlist.RData"))
    
    if(length(resultlist[[length(resultlist)]][[3]]==it)){
      ii.start <- length(resultlist) + 1
      
    }else{
      ii.start <- length(resultlist) 
      
    }
  }
  ## we now test 1 to 3 subclones in each tumor
  result.counter <- ii.start - 1
  for(ii in ii.start:length(n.clones.prim.list)){
    n.clones.prim <- n.clones.prim.list[ii]
    n.clones.rec <- n.clones.rec.list[ii]
    
    print(paste(n.clones.prim," primary,", n.clones.rec, " relapse"))
    
    for(j in 1:length(which.prim[[ii]])){
      result.counter <- result.counter + 1
      
      tree.for.fit <-  tree.matrix(tree=trees[ii], add.prim.clone = F, add.rec.clone = F, which.clone.add = NA, 
                                   which.prim = which.prim[[ii]][[j]], which.rec = which.rec[[ii]][[j]])
      
      fits <- foreach(i=1:it, .combine=rbind) %dopar%{
        set.seed(i)
        
        parms.prim <- Rand.parms(n.clones.prim + 1)
        parms.rec <- Rand.parms(n.clones.rec + 1)
        
        
        res <- EM.alg(tree = tree.for.fit$tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = parms.prim, parms.rec = parms.rec,
                      coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, n.clones.prim = n.clones.prim, n.clones.rec = n.clones.rec, cutoff,
                      sex = sex, maxit = 50, which.prim = which.prim[[ii]][[j]], which.rec = which.rec[[ii]][[j]])
        
        parms.prim <- res$prim
        parms.rec <- res$rec
        
        
        res <- EM.alg(tree = tree.for.fit$tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = parms.prim, parms.rec = parms.rec,
                      coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, n.clones.prim = n.clones.prim, n.clones.rec = n.clones.rec, cutoff,
                      sex = sex, maxit = 50, which.prim = which.prim[[ii]][[j]], which.rec = which.rec[[ii]][[j]])
        if(!is.na(res[[1]])){
          if(res[[1]][1]!="PT_impurity"){
            if(res[[1]][1]!="RT_impurity"){
              res. <- Expectation(tree = tree.for.fit$tree, read.count.matrix.prim = read.count.matrix.prim, read.count.matrix.rec = read.count.matrix.rec, parms.prim = res$prim, 
                                  parms.rec = res$rec, coverage.ratio.baf.prim = coverage.ratio.baf.prim, coverage.ratio.baf.rec = coverage.ratio.baf.rec, locus = locus, 
                                  sex = sex, which.prim = which.prim[[ii]][[j]], which.rec = which.rec[[ii]][[j]])
              
              
              if(!any(c(res$prim, res$rec)<0) & !any(res.$readcounts.prim<0) & !any(res.$readcounts.rec<0)){
                
                
                like <- NegLogLikelihood(readcounts =  res.$readcounts.prim, 
                                         copy.number = res.$rel.mutated.alleles.prim*res.$copy.number.prim/
                                           colSums(res$prim*t(res.$copy.number.prim[,seq(1, n.clones.prim*2+2,2),drop=F])),
                                         parms = rep(res$prim, each=2)) +
                  NegLogLikelihood(readcounts =  res.$readcounts.rec, 
                                   copy.number = res.$rel.mutated.alleles.rec*res.$copy.number.rec/
                                     colSums(res$rec*t(res.$copy.number.rec[,seq(1, n.clones.rec*2+2,2),drop=F])),
                                   parms = rep(res$rec, each=2))
                
                
                return.value <- c(res$prim, res$rec, like)
                
              }else{
                return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
              }
              rm(res.)
              gc()
            }else{
              return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
            }
          }else{
            return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
          }
        }else{
          return.value <- c(rep(0, length(parms.prim) + length(parms.rec)), Inf)
        }
        
        return.value
        
      }
      
      tree. <- tree.for.fit$tree
      tree. <- tree.[,c(tree.for.fit$which.prim, tree.for.fit$which.rec)]
      which.prim. <- 1:length(tree.for.fit$which.prim)
      which.rec. <- (length(tree.for.fit$which.prim)+1):ncol(tree.)
      
      
      
      resultlist[[result.counter]] <- list(Parameters.prim=fits[,1:(n.clones.prim+1)], Parameters.rec=fits[,(n.clones.prim+2):(ncol(fits)-1)], 
                                           Likelihood=fits[,ncol(fits)], tree = tree., which.prim = which.prim.,
                                           which.rec = which.rec.)
      save(resultlist, file=paste0(directory,"/resultlist.RData"))
    }
  }
  
  if(clusterComputation){
    stopImplicitCluster()
    stopCluster(cl)
  }
  
  return(resultlist)
}





#' Function to generate all hypothetical distributions of a mutation among subclones.
#' 
#' This function generates matrices which indicate the combinations of subclones in a tree that can share a mutation and the number of alleles on which this mutation is present.
#' @param tree A binary tree matrix as returned by Tree.matrix.R
#' @param copy.number.prim the copy number of the primary tumor
#' @param copy.number.rec the copy number of the relapse tumor
#' @param which.prim vector indicating which subclones are primary subclones
#' @param which.rec vector indicating which subclones are relapse subclones
#' @param mut.prim set to TRUE if the mutation is present in the primary tumor
#' @param mut.rec set to TRUE if the mutation is present in the relapse tumor
#' @param locus vector containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param b.allele.prim the number of b alleles in the primary tumor
#' @param b.allele.rec the number of b alleles in the relapse tumor
#' @param cnv.indicator either "joined" or "separated", indicating whether the same CNV is present in both, primary and relapsed tumors
#' @param which.prim.cnv binary vector indicating which primary subclones carry the CNV. Set to NA, if there is no CNV at this locus
#' @param which.rec.cnv binary vector indicating which relapse subclones carry the CNV. Set to NA, if there is no CNV at this locus
#' @return Two matrices containing the different possibilities for mutated alleles in primary and relapsed subclones as well as two matrices indicating which subclones carried a copy number change
#' @export
#' @examples 
#' basic.tree <- tree.matrix(tree = "asymmetric", add.prim.clone = F, add.rec.clone = F, which.clone.add = NA, which.prim=c(1,3,4), which.rec = c(2,5), best.tree = NA) # Basic tree of three primary and two relapse specific subclones
#' copy.number.prim <- 2
#' copy.number.rec <- 3
#' which.prim <- basic.tree$which.prim
#' which.rec <- basic.tree$which.rec
#' mut.prim <- T
#' mut.rec <- T
#' locus <- c(1, 111000) # Mutation in both tumors on Chromosome 1
#' sex <- "male"
#' b.allele.prim <- 1
#' b.allele.rec <- 1 ## normal b allele (1) in both tumors
#' cnv.indicator <- "separated"
#' which.prim.cnv <- NA
#' which.rec.cnv <- c(0,1) # subclonal gain in one relapse subclone
#' Index.Matrix(basic.tree, copy.number.prim, copy.number.rec, which.prim, which.rec, mut.prim, mut.rec, locus, sex, b.allele.prim, b.allele.rec, cnv.indicator, which.prim.cnv, which.rec.cnv)


Index.Matrix <- function(tree = "classic", copy.number.prim, copy.number.rec, which.prim, which.rec, mut.prim=T, mut.rec=T, locus,  sex, 
                         b.allele.prim = 1, b.allele.rec = 1, cnv.indicator, which.prim.cnv = NA, which.rec.cnv = NA){
  
  
  s <- tree
  ## remove duplicates, store possibilities of mutated alleles in s and copy numbers in p
  s <- s[!duplicated(s[,c(which.prim,which.rec)]),,drop=F]
  
  p <- s
  
  
  ## if the mutation is only present in the relapse tumor, only take solutions where the two tumors do not share the mutation and vice versa:
  if(mut.prim == F){
    s <- s[rowSums(s[,which.prim,drop=F])==0,,drop=F]
    s <- s[rowSums(s[,which.rec,drop=F])!=0,,drop=F]
  }else if(mut.rec==F){
    s <- s[rowSums(s[,which.rec,drop=F])==0,,drop=F]
    s <- s[rowSums(s[,which.prim,drop=F])!=0,,drop=F]
  }else{
    s <- s[rowSums(s[,which.prim,drop=F])!=0 & rowSums(s[,which.rec,drop=F])!=0,,drop=F]
  }
  
  ## if we have a male patient and a locus on the x chromosome, change the normal copy number to 1
  normal.copy.number <- 2
  b.allele.norm <- 1
  if(locus[1] == "X" && sex == "male"){
    normal.copy.number <- 1
    b.allele.norm <- 0
  }
  ## if only one of the two tumors carries a CNV, only take solutions in agreement with this:
  if((copy.number.prim == normal.copy.number && copy.number.rec != normal.copy.number & b.allele.prim == b.allele.norm) ||
     (b.allele.prim == b.allele.norm && b.allele.rec > 1 & copy.number.prim == normal.copy.number)){
    p <- p[rowSums(p[,which.prim,drop=F])==0,,drop=F]
  }else if((copy.number.prim != normal.copy.number && copy.number.rec == normal.copy.number && b.allele.rec == b.allele.norm)||
           (b.allele.prim >1 && b.allele.rec==b.allele.norm && copy.number.rec == normal.copy.number)){
    p <- p[rowSums(p[,which.rec,drop=F])==0,,drop=F]
  }
  ## if both have the same CNV we need solutions in which the copy number is changed in both 
  
  if(cnv.indicator == "joined"){
    p <- p[rowSums(p[,which.rec,drop=F])!=0 & rowSums(p[,which.prim,drop=F]!=0),,drop=F]
  }
  
  
  if(cnv.indicator != "joined" && ((copy.number.prim != normal.copy.number || b.allele.prim != b.allele.norm) &&
                                   (copy.number.rec != normal.copy.number || b.allele.rec != b.allele.norm))){
    ## only take solutions that are in agreement with one CNV per measured coverage ratio only
    
    p.prim <- p[rowSums(p[,which.rec,drop=F])==0,which.prim,drop=F]
    p.rec <- p[rowSums(p[,which.prim,drop=F])==0,which.rec,drop=F]
    
    if(length(p.prim)==0){
      p.prim <- matrix(0, ncol=length(which.prim), nrow=1)
    }
    if(length(p.rec)==0){
      p.rec <- matrix(0, ncol=length(which.rec), nrow=1)
    }
    
    p.new <- matrix(0, ncol= ncol(p), nrow=nrow(p.prim)*nrow(p.rec))
    for(i in 1:nrow(p.prim)){
      for(j in 1:nrow(p.rec)){
        p.new[(i-1)*nrow(p.rec) + j,][which.prim] <- p.prim[i,]
        p.new[(i-1)*nrow(p.rec) + j,][which.rec] <- p.rec[j,]
      }
    }
    p <- p.new
    rm(p.new)
  }
  ## now intermingle s with all ps
  p.new <- do.call(rbind, replicate(nrow(s), p, simplify=F))
  s.new <- matrix(rep(as.vector(s), each=nrow(p)), ncol=ncol(s))
  s <- s.new
  p <- p.new
  rm(s.new)
  rm(p.new)
  
  
  
  ## only take ps that are in agreement with the previously identified best solutions
  
  if(!is.na(which.prim.cnv[1]) & (copy.number.prim != normal.copy.number | b.allele.prim != b.allele.norm)){
    rows.to.keep.prim <- which(colSums(t((p[,which.prim,drop=F]!=0)*1) == which.prim.cnv)==(length(which.prim)))
  }else{
    rows.to.keep.prim <- seq(1:nrow(p))
  }
  if(!is.na(which.rec.cnv[1]) & (copy.number.rec != normal.copy.number | b.allele.rec != b.allele.norm)){
    rows.to.keep.rec <- which(colSums(t((p[,which.rec,drop=F]!=0)*1) == which.rec.cnv)==(length(which.rec)))
  }else{
    rows.to.keep.rec <- seq(1:nrow(p))
  }
  
  rows.to.keep <- intersect(rows.to.keep.prim, rows.to.keep.rec)
  p <- p[rows.to.keep,,drop=F]
  s <- s[rows.to.keep,,drop=F]
  
  ## remove duplicates
  dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
  
  s <- s[!dupl.rows,,drop=F]
  p <- p[!dupl.rows,,drop=F]
  
  ## in case of normal copy number states: 
  if(((copy.number.prim==normal.copy.number && copy.number.rec==normal.copy.number && b.allele.prim == b.allele.norm && b.allele.rec == b.allele.norm))||
     (copy.number.prim==1 && copy.number.rec==1)){
    ## remove duplicates
    dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
    
    s <- s[!dupl.rows,,drop=F]
    p <- p[!dupl.rows,,drop=F]
    
    
    ## extract the primary and relapse tumor clones
    s.prim <- s[,which.prim,drop=F]
    s.rec <- s[,which.rec,drop=F]
    p.prim <- p[,which.prim, drop=F]
    p.rec <- p[,which.rec, drop=F]
    
    ## returns the mutations and the copy numbers. The number stored in mutations refers to the number of mutated alleles. 1s in the p-matrix indicate that a copy number change happened, 0 indicates that the copy number is normal
    return(list(s.prim = s.prim, s.rec = s.rec, p.prim = p.prim, p.rec = p.rec))
  }
  
  ## now, if the copy number is abnormal in more clones than mutations happened, but at least 1 of the clones carries the mutation and a CNV,
  ## then the CNV happened before the mutation. Only 1 allele is thus affected and nothing needs to be changed in our matrix
  ## if the copy number is abnormal in fewer clones than the mutations, but at least 1 clone shares both, then the CNV happened after the mutation, such that 
  ## all copies of the mutated allele carry the mutation (= either 1 or copy number-1)
  ## case A: same copy number in both tumors or only CNV in one of them.
  
  ## CNV after mutation --> mutation is either on all A alleles or on all B alleles 
  
  cand.clones <- (rowSums(s) > rowSums(p)) & (rowSums(s == p)>0)
  
  if(length(cand.clones)>0){
    row.counter <- nrow(s)+1
    s <- rbind(s, s[which(cand.clones==T),])
    p <- rbind(p, p[which(cand.clones==T),])
    
    for(j in which(cand.clones==T)){
      ## A allele
      
      if(copy.number.prim - b.allele.prim > 0){
        s[row.counter, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- copy.number.prim - b.allele.prim
      }else{
        s[row.counter, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- b.allele.prim
      }
      if(copy.number.rec - b.allele.rec > 0){
        s[row.counter, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <- copy.number.rec - b.allele.rec
      }else{
        s[row.counter, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <- b.allele.rec
      }
      
      row.counter <- row.counter + 1
      
      ## B allele
      if(b.allele.prim > 0){
        s[j, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- b.allele.prim
      }else{
        s[j, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- copy.number.prim - b.allele.prim
      }
      if(b.allele.rec > 0){
        s[j, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <- b.allele.rec
      }else{
        s[j, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <- copy.number.rec - b.allele.rec
      }
      
    }
    ## remove duplicates
    dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
    
    s <- s[!dupl.rows,,drop=F]
    p <- p[!dupl.rows,,drop=F]  
  }
  
  ## case B: two different CNVs in primary and relapse tumors and mutation before copy number change
  
  if(cnv.indicator != "joined"){
    cand.clones.prim <- (rowSums(s[,which.prim,drop=F]) > rowSums(p[,which.prim,drop=F]) & rowSums(s[,which.prim,drop=F]==p[,which.prim,drop=F])>0)
    cand.clones.rec <- (rowSums(s[,which.rec,drop=F]) > rowSums(p[,which.rec,drop=F]) & rowSums(s[,which.rec,drop=F]==p[,which.rec,drop=F])>0)
    ## this only matters if copy.number.prim - b.allele.prim or b.allele.prim is > 1:
    if(copy.number.prim - b.allele.prim > 1 | b.allele.prim > 1){
      if(length(cand.clones.prim)>0){
        row.counter <- nrow(s)+1
        s <- rbind(s, s[which(cand.clones.prim==T),])
        p <- rbind(p, p[which(cand.clones.prim==T),])
        for(j in which(cand.clones.prim==T)){
          ## A allele
          if(copy.number.prim-b.allele.prim > 0){
            s[row.counter, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- copy.number.prim - b.allele.prim
          }else{
            s[row.counter, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- b.allele.prim
          }
          
          row.counter <- row.counter + 1
          ## B allele
          if(b.allele.prim > 0){
            s[j, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- b.allele.prim
          }else{
            s[j, which.prim][p[j,which.prim]!=0 & s[j,which.prim]!=0] <- copy.number.prim - b.allele.prim
          }
          
        }
      }
    }
    
    if(copy.number.rec - b.allele.rec > 1 | b.allele.rec > 1){
      if(length(cand.clones.rec)>0){
        row.counter <- nrow(s)+1
        s <- rbind(s, s[which(cand.clones.rec==T),])
        p <- rbind(p, p[which(cand.clones.rec==T),])
        for(j in which(cand.clones.rec==T)){
          ## A allele
          if(copy.number.rec - b.allele.rec > 0){
            s[row.counter, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <-  copy.number.rec - b.allele.rec
          }else{
            s[row.counter, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <-  b.allele.rec
            
          }
          row.counter <- row.counter + 1
          
          ## B allele
          if(b.allele.rec > 0){
            s[j, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <-  b.allele.rec
          }else{
            s[j, which.rec][p[j,which.rec]!=0 & s[j,which.rec]!=0] <- copy.number.rec - b.allele.rec
            
          }
          
        }
      }
    }
    
    
    ## remove duplicates
    dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
    
    s <- s[!dupl.rows,,drop=F]
    p <- p[!dupl.rows,,drop=F]  
    
  }
  
  
  
  ## case C: unclear what came first --> mutation either on 1 copy or several copies, max. all of either A or B allele
  
  ## if the copy number is abnormal in and only in the clones that carry the mutation, then either it happened before, such that the mutation is present on all alleles
  ## or it happened afterwards, such that it is only present on 1 allele. Or in between (if several copy number changes happened), such that part of them carry it
  if(copy.number.prim != normal.copy.number || copy.number.rec != normal.copy.number |
     b.allele.prim != b.allele.norm || b.allele.rec != b.allele.norm){
    
    if((copy.number.prim > normal.copy.number && copy.number.rec > normal.copy.number && copy.number.prim == copy.number.rec && b.allele.prim == b.allele.rec)||
       (copy.number.prim == normal.copy.number && copy.number.rec == normal.copy.number && b.allele.prim == b.allele.rec && b.allele.prim != b.allele.norm)){
      cand.clones <- rowSums(s == p)==ncol(s)
      
      for(j in which(cand.clones==T)){
        row.counter <- nrow(s)+1
        if(copy.number.prim>1){
          s <- rbind(s, matrix(s[j,], nrow= max(copy.number.prim - b.allele.prim, b.allele.prim),ncol = ncol(s),
                               byrow=T))
          p <- rbind(p, matrix(p[j,], nrow=max(copy.number.prim - b.allele.prim, b.allele.prim),ncol = ncol(s),
                               byrow=T))
          for(k in 1: max(copy.number.prim - b.allele.prim, b.allele.prim)){
            s[row.counter,] <- replace(s[row.counter,], s[row.counter,]!=0,k)
            row.counter <- row.counter + 1
          }
        }
        
        
        
      }
      
      ## remove duplicates
      dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
      
      s <- s[!dupl.rows,,drop=F]
      p <- p[!dupl.rows,,drop=F]  
    }else if(((copy.number.prim != normal.copy.number || copy.number.rec != normal.copy.number) && copy.number.prim != copy.number.rec) ||
             ((b.allele.prim != b.allele.norm || b.allele.rec != b.allele.norm) && b.allele.prim != b.allele.rec)){
      cand.clones.prim <- rowSums(s[,which.prim,drop=F]==p[,which.prim,drop=F])==length(which.prim)
      cand.clones.rec <- rowSums(s[,which.rec,drop=F]==p[,which.rec,drop=F])==length(which.rec)
      for(j in which(cand.clones.prim==T)){
        
        row.counter <- nrow(s)+1
        if(copy.number.prim > 1){
          s <- rbind(s, matrix(s[j,], nrow= max(copy.number.prim - b.allele.prim, b.allele.prim),ncol = ncol(s),
                               byrow=T))
          p <- rbind(p, matrix(p[j,], nrow=max(copy.number.prim - b.allele.prim, b.allele.prim),ncol = ncol(s),
                               byrow=T))
          for(k in 1: max(copy.number.prim - b.allele.prim, b.allele.prim)){
            s[row.counter,] <- replace(s[row.counter,], which.prim[s[row.counter,which.prim]!=0],k)
            row.counter <- row.counter + 1
          }
        }
        
        
        
        
        
      }
      for(j in which(cand.clones.rec==T)){
        
        row.counter <- nrow(s)+1
        if(copy.number.rec > 1){
          s <- rbind(s, matrix(s[j,], nrow= max(copy.number.rec - b.allele.rec, b.allele.rec),ncol = ncol(s),
                               byrow=T))
          p <- rbind(p, matrix(p[j,], nrow=max(copy.number.rec - b.allele.rec, b.allele.rec),ncol = ncol(s),
                               byrow=T))
          for(k in 1:max(copy.number.rec - b.allele.rec, b.allele.rec)){
            s[row.counter,] <- replace(s[row.counter,], which.rec[s[j,which.rec]!=0],k)
            row.counter <- row.counter + 1
          }
        }
      }
      
      ## remove duplicates
      dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
      
      s <- s[!dupl.rows,,drop=F]
      p <- p[!dupl.rows,,drop=F]  
    }
  }
  
  
  ## if copy number = 0 (full deletion), remove solutions where full deletion happened in all subclones (we also found a mutation here, so it needs to exist)
  ## also remove the solutions where the mutation is present on the deleted allele only
  if(copy.number.prim==0 && mut.prim ==T){
    s <- s[rowSums(p[,which.prim,drop=F])!=ncol(p[,which.prim,drop=F]),,drop=F]
    p <- p[rowSums(p[,which.prim,drop=F])!=ncol(p[,which.prim,drop=F]),,drop=F]
    mut.only.at.del <- apply(cbind(s[,which.prim,drop=F],p[,which.prim,drop=F]), 1, function(x){
      ifelse(sum((x[1:ncol(s[,which.prim,drop=F])]!=0)==(x[(ncol(s[,which.prim,drop=F])+1):(2*ncol(s[,which.prim,drop=F]))]!=0))==ncol(s[,which.prim,drop=F]),F,T)
    })
    s <- s[mut.only.at.del,,drop=F]
    p <- p[mut.only.at.del,,drop=F]
    p <- p[rowSums(s[,which.prim,drop=F])!=0,,drop=F]
    s <- s[rowSums(s[,which.prim,drop=F])!=0,,drop=F]
  }
  
  if(copy.number.rec==0 && mut.rec==T){
    s <- s[rowSums(p[,which.rec,drop=F])!=ncol(p[,which.rec,drop=F]),,drop=F]
    p <- p[rowSums(p[,which.rec,drop=F])!=ncol(p[,which.rec,drop=F]),,drop=F]
    mut.only.at.del <- apply(cbind(s[,which.rec,drop=F],p[,which.rec,drop=F]), 1, function(x){
      ifelse(sum((x[1:ncol(s[,which.rec,drop=F])]!=0)==(x[(ncol(s[,which.rec,drop=F])+1):(2*ncol(s[,which.rec,drop=F]))]!=0))==ncol(s[,which.rec,drop=F]),F,T)
    })
    s <- s[mut.only.at.del,,drop=F]
    p <- p[mut.only.at.del,,drop=F]
    p <- p[rowSums(s[,which.rec,drop=F])!=0,,drop=F]
    s <- s[rowSums(s[,which.rec,drop=F])!=0,,drop=F]
  }
  
  
  ## remove duplicates
  dupl.rows <- duplicated(cbind(s[,c(which.prim, which.rec),drop=F], p[,c(which.prim, which.rec),drop=F]))
  
  s <- s[!dupl.rows,,drop=F]
  p <- p[!dupl.rows,,drop=F]
  
  ## remove rows without CNVs although there should be one
  if(copy.number.prim != normal.copy.number){
    s <- s[rowSums(p[,which.prim,drop=F])!=0,,drop=F]
    p <- p[rowSums(p[,which.prim,drop=F])!=0,,drop=F]
  }
  if(copy.number.rec != normal.copy.number){
    s <- s[rowSums(p[,which.rec,drop=F])!=0,,drop=F]
    p <- p[rowSums(p[,which.rec,drop=F])!=0,,drop=F]
  }
  
  ## extract the primary and relapse tumor clones
  s.prim <- s[,which.prim,drop=F]
  s.rec <- s[,which.rec,drop=F]
  p.prim <- p[,which.prim, drop=F]
  p.rec <- p[,which.rec, drop=F]
  ## returns the mutations and the CNVs The number stored in mutations refers to the number of mutated alleles. 1s in the CNV-matrix indicate CNVs, 0s indicate normal copy number states
  
  return(list(s.prim = s.prim, s.rec = s.rec, p.prim = p.prim, p.rec = p.rec))
}





#' Computes a modified Bayesian Information Criterion (compare Chen, 2008, Bioinformatika)
#' 
#' @param logL the log-Likelihood 
#' @param n.parms the number of parameters 
#' @param n.data the number of data points
#' @param Tau a measure of model complexity
#' @param gamma how strictly should increasing model complexity be penalized? (value between 0 and 1)
#' @return The value of the modified Bayesian Information Criterion
#' @export

BIC <- function(logL, n.parms, n.data,  Tau = Tau, gamma = 0.1){
  res <- -2*logL + log(n.data)*n.parms + 2*gamma * Tau
  return(res)
}



#' Function to infer the most likely copy number state per subclone.
#' 
#' This function infers the most likely combination of CNVs in primary/relapsed tumors based on a given tree structure and measaured copy number ratios / B allele frequencies. Only one CNV per tumor sample can be inferred.
#' @param tree A binary tree matrix as returned by Tree.matrix.R
#' @param coverage.ratio.prim measured coverage number ratio of the primary tumor
#' @param coverage.ratio.rec measured coverage number ratio of the relapse tumor
#' @param which.prim vector indicating which subclones are primary subclones
#' @param which.rec vector indicating which subclones are relapse subclones
#' @param locus vector containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param baf.prim measured B allele fraction in the primary tumor
#' @param baf.rec measured B allele fraction in the relapse tumor
#' @param parms.prim vector containing the relative subclone sizes of the primary tumor
#' @param parms.rec vector containing the relative subclone sizes of the relapse tumor
#' @return The best solution if there was a single CNV in both tumors and two different CNVs (one in each tumor). For each solution the copy number states, the B alleles, vectors indicating the subclones carrying the CNV, the squared errors and an indicator indicating whether it is the joined or the separate solution are returned.
#' @export
#' @examples 
#' basic.tree <- tree.matrix(tree = "asymmetric", add.prim.clone = F, add.rec.clone = F, which.clone.add = NA, which.prim=c(1,3,4), which.rec = c(2,5), best.tree = NA) # Basic tree of three primary and two relapse specific subclones
#' coverage.ratio.prim <- 1.5
#' coverage.ratio.rec <- 1.3
#' which.prim <- c(1,2,3)
#' which.rec <- c(4,5) 
#' locus <- c(1, 111000) # Mutation in both tumors on Chromosome 1
#' sex <- "male"
#' baf.prim <- 0.7
#' baf.rec <- 0.6 ## elevated copy number ratios and B allele fractions in both tumors
#' InferCopyNumber(basic.tree, coverage.ratio.prim, coverage.ratio.rec, which.prim, which.rec, locus, sex, baf.prim, baf.rec, parms.prim, parms.rec)


InferCopyNumber <- function(tree, coverage.ratio.prim, coverage.ratio.rec, which.prim, which.rec, locus, sex, baf.prim, baf.rec, parms.prim, parms.rec){
  
  ## start with 0 copynumber:
  if(locus[1]=="Y" | (locus[1]=="X" & sex == "male")){
    copynumber <- 1
    copynumber. <- 0
    b.allele <- 0
    b.allele. <- 0
    b.allele.norm <- 0
    copynumber.norm <- 1
  }else{
    copynumber <- 1
    copynumber. <- 0
    b.allele <- 1
    b.allele. <- 1
    b.allele.norm <- 1
    copynumber.norm <- 2
  }
  
  
  ## We consider all possibilities: same CNV in both tumors, a CNV in only one of them, different CNVs in both of them
  ## The possibility with the same CNV has to be fit separately from the rest. We start with this possibility (measured ratios between 0.9 and 1.1 are considered normal copy number states):
  if((coverage.ratio.prim > 1.1 & coverage.ratio.rec > 1.1 )| (coverage.ratio.prim < 0.9 & coverage.ratio.rec < 0.9)| 
     (!is.na(baf.prim) & !is.na(baf.rec) && ((baf.prim < b.allele/copynumber.norm & baf.rec < b.allele/copynumber.norm) | 
                                             (baf.prim > b.allele/copynumber.norm & baf.rec > b.allele/copynumber.norm)))){
    
    if(coverage.ratio.prim < 0.9){
      copynumber <- 1
    }else{
      copynumber <- copynumber.norm
    }
    
    ## get the possibilities how CNVs can be distributed on the tree from the Index.Matrix function
    s.p.prim.rec <- Index.Matrix(tree = tree, copy.number.prim = 1, copy.number.rec = 1, which.prim = which.prim, which.rec = which.rec, mut.prim = T, mut.rec = T, locus = 1,
                                 sex = sex, cnv.indicator = "joined")
    p.prim <- s.p.prim.rec$p.prim
    p.rec <- s.p.prim.rec$p.rec
    ## only take unique combis of p.prim, p.rec
    to.remove <- duplicated(cbind(p.prim, p.rec)) | 
      rowSums(p.prim)== 0 | rowSums(p.rec)==0
    p.prim <- p.prim[!to.remove,, drop=F]
    p.rec <- p.rec[!to.remove,,drop=F]
    
    
    ## initialize the parameter estimation (we start at normal copy numbers + 1, test all possibilities where the CNV was, take the most likely one. Then we increase the copy number again and repeat this with).
    ## We repeat the procedure until convergence
    squared.error.cn <- Inf
    which.prim.cnv <- p.prim[1,]
    which.rec.cnv <- p.rec[1,]
    
    while(copynumber <= 100 ){
      squared.error <- Inf
      best.sol <- 1
      
      for(j in 1:nrow(p.prim)){
        ## squared error of the copy number change
        squared.error. <- ((sum(p.prim[j,]*copynumber*parms.prim[-length(parms.prim)]) + sum(c(1-p.prim[j,],1)*parms.prim*copynumber.norm))/copynumber.norm - coverage.ratio.prim[1])^2 + 
          ((sum(p.rec[j,]*copynumber*parms.rec[-length(parms.rec)]) + sum(c(1-p.rec[j,],1)*parms.rec*copynumber.norm))/copynumber.norm - coverage.ratio.rec[1])^2
        
        ## test all possible b.allele numbers
        if(!is.na(baf.prim) && !is.na(baf.rec)){
          b.alleles <- seq(0, copynumber)
          baf.prims <- rep(0, length(b.alleles))
          if(sum(c(p.prim[j,],1)*parms.prim)!=0){
            baf.prims <- sapply(b.alleles, function(x){
              (sum(p.prim[j,]*x*parms.prim[-length(parms.prim)]) + sum((1-p.prim[j,])*b.allele.norm*parms.prim[-length(parms.prim)]) + b.allele.norm*parms.prim[length(parms.prim)])/
                (sum(p.prim[j,]*copynumber*parms.prim[-length(parms.prim)]) + sum((1-p.prim[j,])*copynumber.norm*parms.prim[-length(parms.prim)]) + copynumber.norm*parms.prim[length(parms.prim)])
            })
          }
          baf.recs <- rep(0, length(b.alleles))
          if(sum(c(p.rec[j,],1)*parms.rec)!=0){
            
            baf.recs <- sapply(b.alleles, function(x){
              (sum(p.rec[j,]*x*parms.rec[-length(parms.rec)]) + sum((1-p.rec[j,])*b.allele.norm*parms.rec[-length(parms.rec)]) + b.allele.norm*parms.rec[length(parms.rec)])/
                (sum(p.rec[j,]*copynumber*parms.rec[-length(parms.rec)]) + sum((1-p.rec[j,])*copynumber.norm*parms.rec[-length(parms.rec)]) + copynumber.norm*parms.rec[length(parms.rec)])})
          }
          
          b.allele.sq.err <- (baf.prims - baf.prim)^2 + (baf.recs - baf.rec)^2
          squared.error. <- squared.error. + min(b.allele.sq.err)
          
        }else{
          b.allele.sq.err <- NA
        }
        
        
        
        if(squared.error. < squared.error){
          best.sol <- j
          squared.error <- squared.error.
          if(!is.na(b.allele.sq.err)){
            b.allele.. <- b.alleles[which.min(b.allele.sq.err)]
          }else{
            b.allele.. <- 1
          }
          
        }
        
      }
      
      ## level of resolution: 10 % copy number change of +-1, corresponds to a squared error difference of 0.01 (default is no change)
      
      
      if(squared.error.cn > squared.error & abs(squared.error.cn - squared.error) >= 0.01){
        squared.error.cn <- squared.error
        which.prim.cnv <- p.prim[best.sol,]
        which.rec.cnv <- p.rec[best.sol,]
        copynumber. <- copynumber
        b.allele. <- b.allele..
        
      }else{
        break
      }
      
      copynumber <- copynumber + 1
      
    }
    
    copy.number.prim <- copynumber.
    copy.number.rec <- copynumber.
    b.allele.prim <- b.allele.
    b.allele.rec <- b.allele.
    
    
    if(copy.number.prim == copynumber.norm & b.allele.prim == b.allele.norm){
      which.prim.cnv <- NA
      which.rec.cnv <- NA
    }
    
    if(is.na(baf.prim)){
      b.allele.prim <- 1
      b.allele.rec <- 1
    }
    
    solution1 <- list(copy.number.prim = copy.number.prim, copy.number.rec=copy.number.rec, b.allele.prim = b.allele.prim, b.allele.rec = b.allele.rec,
                      which.prim.cnv = which.prim.cnv, which.rec.cnv = which.rec.cnv, squared.error.cn = squared.error.cn, indicator = "joined")
    
  }else{
    solution1 <- list(copy.number.prim = NA, copy.number.rec=NA, b.allele.prim = NA, b.allele.rec = NA,
                      which.prim.cnv = NA, which.rec.cnv = NA, squared.error.cn = Inf, indicator = "joined")
    
  }
  
  
  ## In the second part we estimate primary and relapse tumors separately
  
  s.p.prim.rec <- Index.Matrix(tree = tree, copy.number.prim = 1, copy.number.rec = 1, which.prim = which.prim, which.rec = which.rec, mut.prim = T, mut.rec = T, locus = 1,
                               sex = sex, cnv.indicator = "separated")
  p.prim <- s.p.prim.rec$p.prim
  p.rec <- s.p.prim.rec$p.rec
  
  
  ## only take unique combis of p.prim, p.rec
  to.remove <- duplicated(cbind(p.prim, p.rec)) | 
    rowSums(p.prim)== 0 | rowSums(p.rec)==0
  p.prim <- p.prim[!to.remove,, drop=F]
  p.rec <- p.rec[!to.remove,,drop=F]
  
  
  ## start with 0 copynumber:
  
  
  if(coverage.ratio.prim > 0.9){
    copynumber.prim <- copynumber.norm
  }else{
    copynumber.prim <- 1
  }
  
  if(coverage.ratio.rec > 0.9){
    copynumber.rec <- copynumber.norm
  }else{
    copynumber.rec <- 1
  }
  
  
  
  break.prim <- F
  break.rec <- F
  
  squared.error.cn.prim <- Inf
  squared.error.cn.rec <- Inf
  which.prim.cnv <- p.prim[1,]
  which.rec.cnv <- p.rec[1,]
  while(copynumber.prim <=100  | copynumber.rec <=100){
    squared.error.prim <- Inf
    squared.error.rec <- Inf
    best.sol.prim <- 1
    best.sol.rec <- 1
    
    
    for(j in 1:nrow(p.prim)){
      squared.error.prim. <- ((sum(p.prim[j,]*copynumber.prim*parms.prim[-length(parms.prim)]) + sum(c(1-p.prim[j,],1)*parms.prim*copynumber.norm))/copynumber.norm - coverage.ratio.prim[1])^2 
      squared.error.rec. <-  ((sum(p.rec[j,]*copynumber.rec*parms.rec[-length(parms.rec)]) + sum(c(1-p.rec[j,],1)*parms.rec*copynumber.norm))/copynumber.norm - coverage.ratio.rec[1])^2
      
      ## test all possible b.allele numbers
      if(!is.na(baf.prim)){
        b.alleles.prim <- seq(0, copynumber.prim)
        baf.prims <- rep(0, length(b.alleles.prim))
        if(sum(c(p.prim[j,],1)*parms.prim)!=0){
          baf.prims <- sapply(b.alleles.prim, function(x){
            (sum(p.prim[j,]*x*parms.prim[-length(parms.prim)]) + sum((1-p.prim[j,])*b.allele.norm*parms.prim[-length(parms.prim)]) +  b.allele.norm*parms.prim[length(parms.prim)])/
              (sum(p.prim[j,]*copynumber.prim*parms.prim[-length(parms.prim)]) + sum((1-p.prim[j,])*copynumber.norm*parms.prim[-length(parms.prim)])  +copynumber.norm*parms.prim[length(parms.prim)])
          })
        }
        b.allele.sq.err.prim <- (baf.prims - baf.prim)^2
        squared.error.prim. <- squared.error.prim. + min(b.allele.sq.err.prim)
        
      }else{
        b.allele.sq.err.prim <- NA
      }
      
      if(!is.na(baf.rec)){
        b.alleles.rec <- seq(0, copynumber.rec)
        baf.recs <- rep(0, length(b.alleles.rec))
        if(sum(c(p.rec[j,],1)*parms.rec)!=0){
          
          baf.recs <- sapply(b.alleles.rec, function(x){
            (sum(p.rec[j,]*x*parms.rec[-length(parms.rec)]) + sum((1-p.rec[j,])*b.allele.norm*parms.rec[-length(parms.rec)]) + b.allele.norm*parms.rec[length(parms.rec)])/
              (sum(p.rec[j,]*copynumber.rec*parms.rec[-length(parms.rec)]) + sum(p.rec[j,]*copynumber.norm*parms.rec[-length(parms.rec)]) + copynumber.norm*parms.rec[length(parms.rec)])})
        }
        b.allele.sq.err.rec <-  (baf.recs - baf.rec)^2
        squared.error.rec. <- squared.error.rec. + min(b.allele.sq.err.rec)
        
      }else{
        b.allele.sq.err.rec <- NA
      }
      
      
      
      if(squared.error.prim. < squared.error.prim){
        best.sol.prim <- j
        squared.error.prim <- squared.error.prim.
        if(!is.na(b.allele.sq.err.prim)){
          b.allele.prim.. <- b.alleles.prim[which.min(b.allele.sq.err.prim)]
          
        }else{
          b.allele.prim.. <- 1
        }
      }
      
      if(squared.error.rec. < squared.error.rec){
        best.sol.rec <- j
        squared.error.rec <- squared.error.rec.
        if(!is.na(b.allele.sq.err.rec)){
          b.allele.rec.. <- b.alleles.rec[which.min(b.allele.sq.err.rec)]
          
        }else{
          b.allele.rec.. <- 1
        }
      }
      
    }
    
    if(squared.error.rec >= squared.error.cn.rec &
       squared.error.prim >= squared.error.cn.prim ){
      break
    }
    if(squared.error.rec >= squared.error.cn.rec | abs(squared.error.cn.rec - squared.error.rec) <= 0.01){
      break.rec <- T
    }
    if(squared.error.prim >= squared.error.cn.prim  | abs(squared.error.cn.prim - squared.error.prim) <= 0.01){
      break.prim <- T
    }
    
    if(break.prim & break.rec){break}
    ## level of resolution: 10 % copy number change of +-1, corresponds to a squared error difference of 0.01 (default is no change)
    
    if(squared.error.prim < squared.error.cn.prim ){
      squared.error.cn.prim <- squared.error.prim
      which.prim.cnv <- p.prim[best.sol.prim,]
      copynumber.prim. <- copynumber.prim
      copynumber.prim <- copynumber.prim + 1
      b.allele.prim. <- b.allele.prim..
    }
    
    
    if(squared.error.rec < squared.error.cn.rec ){
      squared.error.cn.rec <- squared.error.rec
      which.rec.cnv <- p.rec[best.sol.rec,]
      copynumber.rec. <- copynumber.rec
      copynumber.rec <- copynumber.rec + 1
      b.allele.rec. <- b.allele.rec..
    }
    
    
    
    
  }
  
  if(coverage.ratio.prim<=1.1 & coverage.ratio.prim >=0.9 & (is.na(baf.prim) || (baf.prim <= 0.55 & baf.prim >= 0.45))){
    copynumber.prim. <- copynumber.norm
    b.allele.prim. <- b.allele.norm
  }
  
  if(coverage.ratio.rec<=1.1 & coverage.ratio.rec >=0.9 & (is.na(baf.rec) || (baf.rec <= 0.55 & baf.rec >= 0.45))){
    copynumber.rec. <- copynumber.norm
    b.allele.rec. <- b.allele.norm
  }
  
  copy.number.prim <- copynumber.prim.
  copy.number.rec <- copynumber.rec.
  b.allele.prim <- b.allele.prim.
  b.allele.rec <- b.allele.rec.
  
  if(copy.number.prim == copynumber.norm & b.allele.prim == b.allele.norm){
    which.prim.cnv <- NA
  }
  
  if(copy.number.rec == copynumber.norm & b.allele.rec == b.allele.norm){
    which.rec.cnv <- NA
  }
  
  indicator = "separated"
  
  
  if(is.na(baf.prim)){
    b.allele.prim <- 1
    b.allele.rec <- 1
  }
  
  
  solution2 <- list(copy.number.prim = copy.number.prim, copy.number.rec=copy.number.rec, b.allele.prim = b.allele.prim, b.allele.rec = b.allele.rec,
                    which.prim.cnv = which.prim.cnv, which.rec.cnv = which.rec.cnv, squared.error.cn = squared.error.cn.prim + squared.error.cn.rec, indicator = indicator)
  
  
  return(list(joined = solution1, separated = solution2))
  
}






#' Binary tree matrices
#' 
#' This function returns binary tree matrices for different tree structures and orders of primary and relapse subclones. For trees of up to three subclones per tumor, a tree matrix can directly chosen from a predefined set. If higher orders are required, these can iteratively be built by providing a basic tree matrix and indicating at which branch an additional subclones should be added. 
#' @param tree should be either of "classic", "symmetric", "asymmetric", "asymmetric2", "asymmetric3", "asymmetric4" and indicates the basic structure of the tree. Should be used if up to three subclones per tumor are desired. If more subclones are requested, set this parameter to NA.
#' @param add.prim.clone if an additional primary subclone should be added to a tree structure, set this parameter to TRUE, else FALSE. Should be FALSE if add.rec.clone is TRUE
#' @param add.rec.clone if an additional relapse subclone should be added to a tree structures, set this parameter to TRUE, else FALSE. Should be FALSE if add.prim.clone is TRUE
#' @param which.clone.add if either an additional primary or relapse subclone should be added, indicate here to which branch the clone should be added (indicating the row of the tree matrix)
#' @param which.prim vector indicating which subclones are primary subclones
#' @param which.rec vector indicating which subclones are relapse subclones
#' @param best.tree tree structure of the basic tree to which an additional subclone whould be added. Set to NA if you want to use a basic tree of up to 3 primary/relapse subclones only
#' @return A binary matrix representing the tree structure (columns correspond to subclones, rows indicate the phylogenetic relationships between subclones), two vectors indicating which subclones are primary and which are relapsed subclones.
#' @export
#' @examples 
#' basic.tree <- tree.matrix(tree = "asymmetric", add.prim.clone = F, add.rec.clone = F, which.clone.add = NA, which.prim=c(1,3,4), which.rec = c(2,5), best.tree = NA) # Basic tree of three primary and two relapse specific subclones
#' extended.tree <- tree.matrix(tree = NA, add.prim.clone = T, add.rec.clone = F, which.clone.add = 5, which.prim = c(1,3,4), which.rec = c(2,5), best.tree = basic.tree) # Add a primary subclone to this basic tree 


tree.matrix <- function(tree, add.prim.clone = F, add.rec.clone = F, which.clone.add = NA, which.prim, which.rec,
                        best.tree){
  if(tree == "classic"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  1, 1, 0, 0, 0, 0,
                  0, 0, 1, 1, 1, 1,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 0, 0, 1, 1), ncol = 6, byrow = T)
  }else if(tree == "symmetric"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  1, 1, 1, 0, 0, 0,
                  0, 1, 1, 0, 0, 0,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 1, 1), ncol = 6, byrow = T)
  }else if(tree == "asymmetric"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  0, 1, 1, 1, 1, 1,
                  0, 0, 1, 1, 1, 1,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 1, 1), ncol = 6, byrow = T)
  }else if(tree == "asymmetric2"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 1, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  0, 0, 1, 1, 1, 1,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 1, 1), ncol = 6, byrow=T)
  }else if(tree == "asymmetric3"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  0, 0, 1, 1, 1, 1,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 0, 0, 1, 1), ncol = 6, byrow = T)
  }else if(tree == "asymmetric4"){
    s <- matrix(c(1, 1, 1, 1, 1, 1,
                  1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1,
                  0, 1, 1, 0, 0, 0,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 1, 1,
                  0, 1, 1, 1, 1, 1), ncol =6, byrow=T)
  }
  
  ## This part is executed if an existing tree (best.tree) shall be extended
  ## The clone can be added at specific branches
  if(add.prim.clone | add.rec.clone){
    
    s <- best.tree
    s <- cbind(s, 0)
    s <- rbind(s, c(rep(0, ncol(s)-1), 1))
    if(which.clone.add > (nrow(s)-1)){
      return("stop")
    }
    ## several possibilities, where to add the new clone
    j <- which.clone.add
    s <- rbind(s,c(s[j,-ncol(s)], 0))
    ## extract the subset, to which the clone will be added
    ## above which clone do we add?
    above.clone <- which(colSums(t(s[,s[j,]!=0,drop=F])==s[j,][s[j,]!=0])==sum(s[j,]))
    above.clone <- above.clone[-length(above.clone)]
    s[above.clone,ncol(s)] <- 1
    s <- unique(s)
    if(add.prim.clone){
      which.prim <- c(which.prim, ncol(s))
    }else{
      which.rec <- c(which.rec, ncol(s))
    }
  }
  
  
  return(list(tree = s, which.prim = which.prim, which.rec = which.rec))
}




#' Best solution per category
#' 
#' This function determines the best fit per candidate tree and subsequently infers the most likely tree
#' @param tree.fitting.result A fitting result as returned by fit.model
#' @param min.non.ambiguous.clonals the minimal fraction of clonal mutations which should fit unambiguously in order to accept a fit
#' @return A matrix containing the likelihood and the modified Bayseian Information Criterion of the most likely fit per candidate tree
#' @export
#' 

Best.solution <- function(tree.fitting.result, min.non.ambiguous.clonals = 0.5, sex="male"){
  best.sol.per.cat <- as.data.frame(matrix(0, nrow= length(tree.fitting.result), ncol = 9, dimnames = list(c(), c("n.parms.prim", "n.parms.rec", "Index", "LL", "mBIC", "mean sq error clonal", "mean sq error all", "ambiguous clonals", "clonals"))))
  best.sol.per.cat[,1] <- unlist(lapply(tree.fitting.result, function(x){length(x$which.prim)}) )
  best.sol.per.cat[,2] <- unlist(lapply(tree.fitting.result, function(x){length(x$which.rec)}) )
  best.sol.per.cat[,c(4:9)] <- Inf
  
  
  
  for(j in 1:length(tree.fitting.result)){
    if(length(tree.fitting.result[[j]]$Parameters.prim)==0){next}
    
    ## we get the logarithmic number of possibilities from running the expectation function
    Tau <- Expectation(tree = tree.fitting.result[[j]]$tree, read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec,
                       parms.prim = tree.fitting.result[[j]]$Parameters.prim[1,], 
                       parms.rec = tree.fitting.result[[j]]$Parameters.rec[1,], 
                       coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim,
                       coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec, locus = example.data$Locus.prim, sex = sex,
                       which.prim = tree.fitting.result[[j]]$which.prim, which.rec = tree.fitting.result[[j]]$which.rec)
    
    
    Tau <- Tau$number.of.possibilities
    
    for(k in 1:nrow(tree.fitting.result[[j]]$Parameters.prim)){
      ## number of non-zero clones
      n.prim <- sum(tree.fitting.result[[j]]$Parameters.prim[k,-ncol(tree.fitting.result[[j]]$Parameters.prim)]!=0)
      n.rec <- sum(tree.fitting.result[[j]]$Parameters.rec[k,-ncol(tree.fitting.result[[j]]$Parameters.rec)]!=0)
      
      mBIC <- BIC(logL = -tree.fitting.result[[j]]$Likelihood[k], n.parms = n.prim + n.rec, n.data = sum(example.data$read.count.matrix.prim) + sum(example.data$read.count.matrix.rec),  Tau = Tau, gamma = 0.9)
      
      
      if(n.prim == best.sol.per.cat[j,1] & n.rec == best.sol.per.cat[j,2]){
        if(mBIC < best.sol.per.cat[j, "mBIC"]){
          
          ## only use this solution if clonal mutations fit well (unambiguous). The information is retained if running the expectation function with option likelihood.space = T. Now, all fits are reported (in case of ambiguous mutations, all solutions that fit equally well are reported if likelihood.space=T)
          sqerr <- Expectation(tree = tree.fitting.result[[j]]$tree, read.count.matrix.prim = example.data$read.count.matrix.prim, read.count.matrix.rec = example.data$read.count.matrix.rec,
                               parms.prim = tree.fitting.result[[j]]$Parameters.prim[k,], 
                               parms.rec = tree.fitting.result[[j]]$Parameters.rec[k,], 
                               coverage.ratio.baf.prim = example.data$Coverage.B.alleles.prim,
                               coverage.ratio.baf.rec = example.data$Coverage.B.alleles.rec, locus = example.data$Locus.prim, sex = sex,
                               which.prim = tree.fitting.result[[j]]$which.prim, which.rec = tree.fitting.result[[j]]$which.rec, likelihood.space = T)
          ## mutations clonal in primary sample
          clonal <- rownames(sqerr$rel.mutated.alleles.prim[rownames(sqerr$readcounts.prim),])[rowSums(sqerr$rel.mutated.alleles.prim[rownames(sqerr$readcounts.prim),seq(2, ncol(sqerr$rel.mutated.alleles.prim),2)]!=0) ==
                                                                                                 (ncol(sqerr$rel.mutated.alleles.prim)/2-1) &
                                                                                                 rowSums(sqerr$readcounts.prim[,seq(2, ncol(sqerr$rel.mutated.alleles.prim),2)])!=0]
          ## mutations clonal in relapse sample
          clonal. <- rownames(sqerr$rel.mutated.alleles.rec[rownames(sqerr$readcounts.rec),])[rowSums(sqerr$rel.mutated.alleles.rec[rownames(sqerr$readcounts.rec),seq(2, ncol(sqerr$rel.mutated.alleles.rec),2)]!=0) ==
                                                                                                (ncol(sqerr$rel.mutated.alleles.rec)/2-1) &
                                                                                                rowSums(sqerr$readcounts.rec[,seq(2, ncol(sqerr$rel.mutated.alleles.rec),2)])!=0]
          ## mutations clonal in both
          clonal <- intersect(clonal, clonal.)
          
          ## look at how many of the clonal genes have an ambiguous solution. Count these as not truly clonal
          non.clonal.prim <- rownames(sqerr$rel.mutated.alleles.prim)[(rowSums(sqerr$rel.mutated.alleles.prim[,seq(2,ncol(sqerr$rel.mutated.alleles.prim),2)]!=0)!=(ncol(sqerr$rel.mutated.alleles.prim)/2-1))]
          non.clonal.rec <-  rownames(sqerr$rel.mutated.alleles.rec)[(rowSums(sqerr$rel.mutated.alleles.rec[,seq(2,ncol(sqerr$rel.mutated.alleles.rec),2)]!=0)!=(ncol(sqerr$rel.mutated.alleles.rec)/2-1))]
          non.clonal.prim <- intersect(clonal, non.clonal.prim)
          non.clonal.rec <- intersect(clonal, non.clonal.rec)
          ## percentage of ambiguous clonal mutations
          percent.non.unique.clonals <- length(unique(c(non.clonal.prim, non.clonal.rec)))/length(clonal)
          
          ## in addition compute the mean squared error of measured vs predicted variant allele frequencies of mutations present in both samples
          mean.sqerr <- sum(sqerr$residuals.vafs.prim[names(sqerr$residuals.vafs.prim) %in% names(sqerr$residuals.vafs.rec)]^2 +
                              sqerr$residuals.vafs.rec[names(sqerr$residuals.vafs.rec) %in% names(sqerr$residuals.vafs.prim)]^2, na.rm = T)/
            sum(names(sqerr$residuals.vafs.rec) %in% names(sqerr$residuals.vafs.prim))
          
          ## plus account for the mean squared errors of the clonal mutations 
          if(length(clonal)==0){
            sqerr <- Inf
          }else{
            sqerr <- sum(sqerr$residuals.vafs.prim[clonal]^2 + sqerr$residuals.vafs.rec[clonal]^2)/length(clonal)
            
          }
          
          ## if more than 50 % of the clonal mutations are ambiguous this is considered a bad fit (threshold can be chosen differently)
          if(!is.na(percent.non.unique.clonals) & !is.infinite(percent.non.unique.clonals) && percent.non.unique.clonals < min.non.ambiguous.clonals){
            best.sol.per.cat[j ,"LL"] <- tree.fitting.result[[j]]$Likelihood[k]
            best.sol.per.cat[j , "mBIC"] <- mBIC
            best.sol.per.cat[j,"Index"] <- k 
            best.sol.per.cat[j ,"mean sq error clonal"] <- sqerr
            best.sol.per.cat[j ,"clonals"] <- length(clonal)
            best.sol.per.cat[j ,"mean sq error all"] <- mean.sqerr
            best.sol.per.cat[j ,"ambiguous clonals"] <- percent.non.unique.clonals
            
          }
        }
      }
    }
  }
  
  best.sol.per.cat <- cbind(best.sol.per.cat, best.sol.per.cat[,5] - min(best.sol.per.cat[,5][best.sol.per.cat[,5] != 0 & best.sol.per.cat[,5] > 0]))
  colnames(best.sol.per.cat)[10] <- "delta.mBIC"
  
  return(best.sol.per.cat)
  
}



