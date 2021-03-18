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
#' @param tree a binary tree matrix; each column corresponds to one subclone; each row to possible relationships between mutations (1: mutation is present in clone, 0: mutation is absent in clone)
#' @param read.count.matrix.mut a matrix storing the number of mutant reads in each sample (row: location, column: sample)
#' @param read.count.matrix.ref a matrix storing the number of reference reads in each sample (row: location, column: sample)
#' @param coverage.ratio a matrix storing the measured coverage ratioper sample
#' @param baf a matrix storing the measured BAF per sample
#' @param locus matrix containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param parms vector containing the relative subclone sizes of each clone, corresponding to the clades in the tree
#' @param sample.indicator vector assigning each clade in the tree to a sample (e.g. c("A", "B", "C", "A", "B) for a tree with 5 clades)
#' @param normal.copy.number.auto normal copy number on autosomes 
#' @param normal.copy.number.y normal copy number on the y chromosome
#' @param normal.b.allele.norm normal B allele count
#' @param joined.info List of lists containing prior information on copy number changes that occurred together. Length of list must correspond to number of loci. For example, if it is known that the copy number change at locus x is present in sample 'A' and 'B' set joined.info[[x]] <- c("A", "B). If there's no information, set joined.info[[x]] <- NA
#' @param false.negatives logical; should the possibility of false negative mutations be accounted for?
#' @param likelihood.space TRUE if additional solutions should be reported. I.e. if a SNV/CNV cannot be uniquely mapped to a specific locus in the tree, additional solutions are reported. In addition, residual sums of squares between measured and estimated variant allele frequencies, B allele frequencies and coverage ratios are reported.
#' @return Two matrices storing the expected readcounts (reference, mutated) at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones; columns represent subclones; last two columns represent normal tissue); two matrices storing the relative number of reference and mutated alleles at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones); two matrices storing the number of alleles at each mutated locus per tumor subclone (one matrix for primary, one for relapse subclones); two matrices indicating which subclones carry a CNV at each mutated locus per tumor subclone. One integer reporting the logarithmic number of possibilities to sort the mutations on the tree. If likelihood.space==T, additional information is provided: additional solutions for copy number changes are stored in separate matrices; additional solutions for SNVs are stored in the same matrices as before (they will be longer than and affected SNVs can be found by checkin which SNV name appears twice among the rownames of the matrix). In addition, residual sums of squares for the variant allele frequencies, the B allele frequences and the coverage ratios are returned
#' @export


Expectation <- function(tree = "classic", read.count.matrix.mut, read.count.matrix.ref, parms, coverage.ratio, baf,
                        locus, sex = sex, likelihood.space = F, sample.indicator, normal.copy.number.auto = 2, normal.copy.number.y = 1, normal.b.allele.norm = 1, joined.info = NA, false.negatives = F){
  
  ## should the entire procedure be aborted in case of a problem?
  break.function <- F
  
  ## If there's no information whether the copy number state is equal across samples, set to NA
  if(all(is.na(joined.info))){
    joined.info <- as.list(rep(NA, nrow(read.count.matrix.mut)))
  }
  
  ## if any parameter is set to zero, remove the clade from the tree
  zero <- F
  if(any(parms ==0)){
    tree <- tree[,-which(parms==0),drop=F]
    tree <- tree[rowSums(tree)>0,,drop=F]
    tree <- unique(tree)
    sample.indicator <- sample.indicator[-which(parms==0)]
    zeros <- which(parms==0)
    parms <- parms[-which(parms==0)]
    zero <- T
  }
  
  
  ## the output will be a readcount matrix (res), containing reference and mutated reads assigned to each subclone, an inidcator matrix (s.res)
  ## containing binary information whether a mutation is present in a subclone or not; a matrix containing information on the copy number state (cn.res) 
  ## and a corresponding binary indicator matrix (cn.res.index) as well as a matrix containing information of the number of reads inferred to be due to contamination by normal tissue (normal.tissue)
  res <- matrix(0, nrow=nrow(read.count.matrix.mut), ncol=length(parms)*2)
  rownames(res) <- rownames(read.count.matrix.mut)
  colnames(res) <- rep(c("reference", "mutated"), ncol(res)/2)
  s.res <- res
  cn.res <- res
  cn.res.index <- res
  normal.tissue <- matrix(0, nrow = nrow(read.count.matrix.mut), ncol=ncol(read.count.matrix.mut))
  rownames(normal.tissue) <- rownames(read.count.matrix.mut)
  
  ## do we need this?
  cn.indicator.all <- rep(0, nrow(cn.res.index))
  names(cn.indicator.all) <- rownames(cn.res.index)
  
  ## if we are not only interested in a point estimate, but in the likelihood space, store the likelihood for each locus separately
  if(likelihood.space){
    like.per.loc <- rep(Inf, nrow(cn.res.index))
  }
  
  if(sum(parms)==0){ ## If all subclones are zero, the other samples are not fitted. Abort
    print("Impurity")
    return("Impurity")}
  
  
  ## store the number of possible combinations per locus
  number.of.possibilities <- log(1)
  
  ## iterate through x genes. Each likelihood can be independently maximized (-logL minimized)
  
  for(i in seq_len(nrow(read.count.matrix.mut))){
    
    
    ## the first thing we have to do is finding out what copy numbers we have to choose. We assume that evolution is simple:
    ## if the coverage ratio estimate indicates a deletion, we test combinations for copy number = 0 and copy number = 1 and choose the most likely one
    ## if the estimate indicates a duplication/ gain, we test combinations for copy number > 2
    ## if both, primary and relapse have a duplication/deletion, we test whether a common event or two separate events are more likely
    
    
    ## check whether you have to estimate the copy number (if it's the same as in the last segment you can reuse it)
    determine.copy.number <- T
    if(i > 1){
      if(all(coverage.ratio[i,]==coverage.ratio[i-1,]) && locus[i,1]==locus[i-1,1] && 
         !any((baf[i,]==baf[i-1,]) %in% c(F, NA))){
        determine.copy.number <- F
      }
    }
    
    ## the normal copy number /BAF is different on male sex chromosomes
    if(sex== "male" & locus[i,1] %in% c("X", "Y")){
      normal.copy.number <- normal.copy.number.y
      normal.baf <- 0
    }else{
      normal.copy.number <- normal.copy.number.auto
      normal.baf <- normal.copy.number.auto/2
    }
    
    
    if(determine.copy.number){
      
      ## assume coverage ratios of 1 +- 5 % as normal
      if(all(coverage.ratio[i,] < 1.05) && all(coverage.ratio[i,] > 0.95)){
        if( (all(is.na(baf[i,])) || all(baf[i,]<(normal.baf/normal.copy.number +0.05))) && 
            (all(is.na(baf[i,])) || all(baf[i,] > (normal.baf/normal.copy.number - 0.05)))){
          
          copy.number.list <- lapply(1:ncol(coverage.ratio), function(x){normal.copy.number})
          b.allele.list <- lapply(1:ncol(coverage.ratio), function(x){normal.baf})
          best.cn.sol.all.list <- NA
          
        }else if(normal.copy.number == 1 && (all(is.na(baf[i,])) || all(baf[i,]<0.05))){
          
          copy.number.list <- lapply(1:ncol(coverage.ratio), function(x){normal.copy.number})
          b.allele.list <- lapply(1:ncol(coverage.ratio), function(x){normal.baf})
          best.cn.sol.all.list <- NA
          
        }else{
          
          inferred.copy.number <- InferCopyNumber(tree = tree, coverage.ratios = coverage.ratio[i,], bafs = baf[i,], sample.indicator = sample.indicator, sex = sex,
                                                  locus = locus[i,], parms = parms, normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm,
                                                  joined.info = joined.info[[i]])
          
          if(inferred.copy.number[[1]][1]=="no fit"){
            return(NA)
          }
          
          
          copy.number.list <- inferred.copy.number[[1]]
          b.allele.list <- inferred.copy.number[[2]]
          best.cn.sol.all.list <- inferred.copy.number[[3]]
          
        }
        
      }else{
        inferred.copy.number <- InferCopyNumber(tree = tree, coverage.ratios = coverage.ratio[i,], bafs = baf[i,], sample.indicator = sample.indicator, sex = sex,
                                                locus = locus[i,], parms = parms, normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm,
                                                joined.info = joined.info[[i]])
        
        if(inferred.copy.number[[1]][1]=="no fit"){
          
          return(NA)
        }
        
        copy.number.list <- inferred.copy.number[[1]]
        b.allele.list <- inferred.copy.number[[2]]
        best.cn.sol.all.list <- inferred.copy.number[[3]]
        
      }
      
      
      
      
    }
    
    ## we don't need to determine the matrix of candidate solutions for a mutation distribution if the previous mutation was present in the same subset of the tumors and was at a locus of equal copy number state
    determine.s.matrix <- T
    if(i > 1){
      if(determine.copy.number==F){
        if(all((read.count.matrix.mut[i,]!=0) == (read.count.matrix.mut[i-1,]!=0))){
          determine.s.matrix <- F
        }
      }
    }
    
    
    
    
    if(determine.s.matrix){
      
      s.p <- Index.Matrix(tree = tree, coverage.ratio = copy.number.list, sample.identifier = sample.indicator, readcounts.mut = read.count.matrix.mut[i,], 
                          locus = locus[i,], sex = sex, best.cn.sol.all = best.cn.sol.all.list, normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, false.negatives = false.negatives)
      
      s <- s.p$s
      p <- s.p$p
      
      
      
      if(length(s)==0){
        print(rownames(read.count.matrix.mut)[i])
        return("Impurity")}
      
      
      rm(s.p)
      
      ## remove rows where no mutations are:
      
      rows.no.mut <- which(rowSums(s[,parms!=0,drop=F])==0)
      
      if(length(rows.no.mut)>0 & rowSums(read.count.matrix.mut)[i]!=0){
        s <- s[-unique(rows.no.mut),,drop=F]
        p <- p[-unique(rows.no.mut),,drop=F]
      }
      
      ## remove duplicated rows based on the clones that are occupied
      dupl.rows <- duplicated(cbind(s[,parms!=0,drop=F], p[,parms!=0,drop=F]))
      
      s <- s[!dupl.rows,,drop=F]
      p <- p[!dupl.rows,,drop=F]
      
      ## remove zero rows based on the clones that are occupied
      zero.rows <- rowSums(s[,parms!=0,drop=F])==0
      
      if(rowSums(read.count.matrix.mut)[i]!=0){
        s <- s[!zero.rows,,drop=F]
        p <- p[!zero.rows,,drop=F]
      }
      
      
      ## store an extra matrix as a copy number indicating matrix (1: copy number change, 0: no copy number change)
      cn.index <- p
      
      if(all(copy.number.list == normal.copy.number) && (all(is.na(b.allele.list)) ||all(b.allele.list == normal.baf))){
        cn.index[,seq_len(ncol(cn.index))] <- 0
      }
      
      ## remove duplicated rows based on the clones that are occupied
      dupl.rows <- duplicated(cbind(s[,parms!=0,drop=F], p[,parms!=0,drop=F]))
      
      s <- s[!dupl.rows,,drop=F]
      p <- p[!dupl.rows,,drop=F]
      cn.index <- cn.index[!dupl.rows,,drop=F]
      
      ## add normal tissue to each clone in the index matrix 
      cn.index <- cn.index[,rep(seq_len(ncol(cn.index)),each=2),drop=F]
      
      ## divide the indicator matrix by the respective copy number
      s <- s/p
      s <- replace(s, is.infinite(s), 0)
      s <- replace(s, is.na(s), 0)
      
      ## add the normal counterparts for each subclone
      s <- s[,rep(seq_len(ncol(s)), each=2),drop=F]
      for(jj in seq(2, ncol(s), 2)){
        s[,jj-1] <- 1-s[,jj]
      }
      p <- p[,rep(seq_len(ncol(p)), each=2),drop=F]
      
      ## also remove solutions where no mutation is present in any sample
      mutated.samples <- names(read.count.matrix.mut[i,])[read.count.matrix.mut[i,]!=0]
      non.mutated.samples <- names(read.count.matrix.mut[i,])[read.count.matrix.mut[i,]==0]
      to.keep <- rep(T, nrow(s))
      for(k in mutated.samples){
        to.keep <- to.keep * (rowSums(s[,which(sample.indicator %in% k)*2, drop=F])!=0)
      }
      
      s <- s[as.logical(to.keep),,drop=F]
      p <- p[as.logical(to.keep),,drop=F]
      cn.index <- cn.index[as.logical(to.keep),,drop=F]
      
      
      ## remove duplicated rows based on the clones that are occupied
      dupl.rows <- duplicated(cbind(s, p))
      
      s <- s[!dupl.rows,,drop=F]
      p <- p[!dupl.rows,,drop=F]
      cn.index <- cn.index[!dupl.rows,,drop=F]
      
    }
    
    ## if no solution reamins, abort
    if(length(s)==0){
      print("here")
      print("Impurity")
      return("Impurity")
    }
    
    if(break.function==T){
      break}
    
    
    ## s is now a matrix that contains the (relative) coverage.ratio within each subclone for reference and mutated alleles. (already normalized by total copy number)
    ## e.g. 0.5 0.5 0.0 1.0 0.0 1.0
    ##      0.0 1.0 0.5 0.5 0.0 1.0
    ##      0.5 0.5 0.5 0.5 0.0 1.0
    ## ..... (each row contains the number of reference reads in SCI, mutated reads in SCI, reference reads in SCII, mutated reads in SCII, ...)
    ## The last 2 columns correspond to normal tissue (always 0 in mutated reads, 1 in reference reads)
    
    ## for each gene iterate through each possible combination of s. Choose the most likely combination and store it in the result matrix
    ## store the likelihood
    LL <- 0
    ## store the most likely result table
    res.most.likely <- 0
    ## store the corresponding membership indices
    s.most.likely <- 0
    
    ## switch: if we are fitting, we don't do this, but later to get an idea about the overall likelihood, not only take into account the most likely solution, but the 95 % most likely ones
    if(likelihood.space){
      LL.all.sol <- rep(0, nrow(s))
    }
    
    number.of.possibilities <- number.of.possibilities + log(nrow(s))
    
    normal.copy.number <- ifelse(sex=="male" & locus[1] %in% c("X", "Y"), normal.copy.number.y, normal.copy.number.auto)
    b.allele.norm <- ifelse(sex=="male" & locus[1] %in% c("X", "Y"), 0, floor(normal.copy.number.auto/2))
    
    
    ## in case of a clonal mutation add the solution that it might be present on all alleles A and all alleles B of the clonal copy number
    if(sum(read.count.matrix.mut[i,]!=0)==ncol(read.count.matrix.mut)){
      if(any(rowSums(s[,seq(2,ncol(s),2),drop=F]!=0)==(ncol(s)/2))){
        for(k in which(rowSums(s[,seq(2,ncol(s),2),drop=F]!=0)==(ncol(s)/2))){
          to.add <- s[k,]*max(b.allele.norm, normal.copy.number - b.allele.norm)
          to.add[to.add > 1] <- 1
          to.add[seq(1,length(to.add),2)] <- 1 - to.add[seq(2, length(to.add),2)]
          s <- rbind(s, to.add)
          p <- rbind(p, p[k,])
          cn.index <- rbind(cn.index, cn.index[k,])
        }
      }
    }
    
    ## iterate through the possible combinations, stored in s
    for(k in seq_len(nrow(s))){
      ## we only need to calculate the expected values if we have more than 1 clone. Otherwise the solution is trivial (see below). We can now optimize separately per sample
      
      res.mle.ref <- rep(0, ncol(s)/2)
      res.mle.mut <- res.mle.ref
      normal.tissue.tmp <- rep(0, ncol(coverage.ratio))
      
      ## go through all combinations that belong together copy-number-wise
      s. <- s[k,seq(2, ncol(s), 2),drop=F]
      p. <- p[k,,drop=F]
      
      
      LL.at.k <- 0
      
      ## extract the inferred copy number state per sample for solutions with joined or independent copy number events
      best.cn.sol.all.list. <- lapply(best.cn.sol.all.list, function(x){if(is.na(x[1])){
        return(rep(0, length(sample.indicator)))}else{return(x)}})
      
      best.cn.sol.all. <- unlist(best.cn.sol.all.list)
      best.cn.sol.all. <- best.cn.sol.all.[!is.na(best.cn.sol.all.)]
      
      
      if(length(best.cn.sol.all.)==0){
        best.cn.sol.all. <- rep(1, length(sample.indicator))
      }
      
      best.cn.sol.all. <- matrix(best.cn.sol.all., byrow = T, ncol=length(sample.indicator))
      best.cn.sol.all. <- unique(best.cn.sol.all.)
      best.cn.sol.all. <- rbind(replace(rep(0, ncol(best.cn.sol.all.)), colSums(best.cn.sol.all.)==0,1), best.cn.sol.all.)
      
      ## iterate through the different possible combinations of distributing the CNVs across the tree
      for(l in 1:nrow(best.cn.sol.all.)){
        s. <- s[k,seq(2, ncol(s), 2),drop=F]
        
        ## Test if there is a CNV co-occuring with a mutation. In this case, the order of the mutation is relevant (SNV < CNV or CNV < SNV)
        ## This leads to an expansion of s.
        ## This is only relevant, if there is a CNV at all and if the CNV is in a clone that also harbors the SNV. The first case in best., l=1 encodes the subclones without a CNV --> we don't need to adjust anything
        ## If the CNV is present in all clones, 
        if(length(best.cn.sol.all.[l,])>0 & any(which(best.cn.sol.all.[l,]!=0) %in% which(s.!=0)) &
           sum(best.cn.sol.all.[l,])!=ncol(best.cn.sol.all.) & l > 1){
          
          ## which clones harbor the copy number change?
          which.clones <- which(best.cn.sol.all.[l,]>0)
          
          ## in which copy number list do we have to look?
          
          if(length(which.clones)==ncol(best.cn.sol.all.) & 
             length(unique(unlist(copy.number.list)[ which(LETTERS %in% sample.indicator[which(best.cn.sol.all.[l,]>0)])]))>1){
            which.copy.number <- which(unlist(copy.number.list)==normal.copy.number)[1]
          }else{
            which.copy.number <- which(LETTERS %in% sample.indicator[which(best.cn.sol.all.[l,]>0)])[1]
            
          }
          
          ## copy number change after mutation --> mutation on 1 (2 if copy number=4 etc) or all alleles
          cnv.after.mutation <- sum(s.[1,]!=0) > sum(best.cn.sol.all.[l,]) & sum((s.!=0)==best.cn.sol.all.[l,])>0
          
          
          ## only applies if we have mutations at all
          if(cnv.after.mutation & rowSums(read.count.matrix.mut)[i]!=0){
            ##duplicate s. to account for the 2nd allele
            s. <-  rbind(s.,s.)
            
            if(-min(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm), 
                    b.allele.list[[which.copy.number]] - b.allele.norm)<max(s.[1, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[1,which.clones]!=0]*copy.number.list[[which.copy.number]])){
              
              ## if the copy number change is a multiple of the original copy number, and the original copy number is not 1 or 2, then the copy number may have changed the VAF of the mutant, as it occurred afterwards
              ## e.g. if the CNV is 3 it may have elevated the mutation to 2/3 instead of 1/3.
              if((copy.number.list[[which.copy.number]]%%normal.copy.number)!=0 | normal.copy.number %in% c(1,2)){
                
                ## A allele
                s.[1, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[1,which.clones]!=0] <-(s.[1, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[1,which.clones]!=0]*copy.number.list[[which.copy.number]]+
                                                                                                      (copy.number.list[[which.copy.number]]-b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm)))/copy.number.list[[which.copy.number]]
                
                
                ## B allele
                
                s.[2, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[2,which.clones]!=0] <- (s.[2, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[2,which.clones]!=0]*copy.number.list[[which.copy.number]] +
                                                                                                       (b.allele.list[[which.copy.number]]-b.allele.norm))/copy.number.list[[which.copy.number]]
                
                
                ## if anything happened to be now < 0, take this solution out
                
                s. <- s.[rowSums(s.>=0)==ncol(s.),,drop=F]
              }
              ## if the original copy.number is 2 and the copy.number change a multiple, additionally allow for 50% 50%
              if((copy.number.list[[which.copy.number]]%%normal.copy.number)==0 & normal.copy.number==2){
                s. <- rbind(s., s.[1,])
                
                ## A allele
                s.[3, which.clones][best.cn.sol.all.[l,which.clones]!=0 & s.[1,which.clones]!=0] <- rep(0.5, length(best.cn.sol.all.[l,which.clones]!=0 & s.[1,which.clones]!=0))
                
                
                ## if anything happened to be now < 0, take this solution out
                s. <- s.[rowSums(s.>=0)==ncol(s.),,drop=F]
              }
              
              
              
            }
            
          }
          
          ## if it is unclear whether the SNV preceded the CNV or vice versa, the mutation can lie on all, 1 and all intermediate steps of the copy number change
          cand.combinations <- which(colSums(t(s.!=0)==best.cn.sol.all.[l,])==ncol(s.))
          
          if(length(cand.combinations & rowSums(read.count.matrix.mut)[i]!=0)>0){
            
            row.counter <- nrow(s.)+1
            s. <- rbind(s., s.[which(cand.combinations==T),])
            
            for(j in which(cand.combinations==T)){
              
              ## in case we don't have mutant copy numbers, don't need to do anything
              if(max(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm),
                     copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - b.allele.norm,
                     b.allele.list[[which.copy.number]] - b.allele.norm)==0){next}
              
              ## make space for each possible state
              s.<- rbind(s., matrix(s.[j,], nrow=max(abs(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm)),
                                                     abs(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - b.allele.norm),
                                                     abs(b.allele.list[[which.copy.number]] - b.allele.norm)), ncol = ncol(s.), byrow=T))
              
              
              ## for gains count upwards, for losses count downwards
              sign <- 1
              
              if((copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm)) <0 &
                 (copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - b.allele.norm) <0 &
                 (b.allele.list[[which.copy.number]] - b.allele.norm)<0){
                sign <- -1
              }
              
              
              for(kk in 1: max(abs(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - (normal.copy.number - b.allele.norm)), 
                               abs(copy.number.list[[which.copy.number]] - b.allele.list[[which.copy.number]] - b.allele.norm),
                               abs(b.allele.list[[which.copy.number]] - b.allele.norm))){
                k. <- sign*kk
                
                
                
                s.[row.counter,which.clones] <- replace(s.[row.counter,which.clones], s.[row.counter,which.clones]!=0,
                                                        (s.[row.counter,which.clones][s.[row.counter,which.clones]!=0]*copy.number.list[[which.copy.number]]+
                                                           k.)/copy.number.list[[which.copy.number]])
                s.[row.counter,which.clones] <- replace(s.[row.counter,which.clones], s.[row.counter,which.clones]>1, 1)
                row.counter <- row.counter + 1
              }
              
              
              s. <- s.[rowSums(s.>=0)==ncol(s.),,drop=F]
              
            }
          }
          
          
          
        }else{
          which.clones <- 1:ncol(best.cn.sol.all.)
        }
        
        
        ## remove duplicates
        dupl.rows <- duplicated(s.)
        
        s. <- s.[!dupl.rows,,drop=F]
        
        
        
        ## now that we have s. ready, we test all samples that show this copy number state (stored in which.clones)
        
        LL.tmp <- 0
        res.mle.ref. <- res.mle.ref
        res.mle.mut. <- res.mle.mut
        normal.tissue. <- rep(0, ncol(coverage.ratio))
        
        s. <- s.[,rep(1:ncol(s.), each=2),drop=F]
        
        s.[,seq(1, ncol(s.)-1,2)] <- 1 - s.[,seq(2, ncol(s.),2)]
        
        p. <- do.call(rbind, replicate(nrow(s.), p., simplify=F))
        
        
        
        for(it in 1:nrow(s.)){
          
          LL.tmp.2 <- 0
          
          ## iterate through all samples corresponding to the copy number state
          for(sample in which(LETTERS %in% sample.indicator[which.clones])){
            
            
            which.columns <- which(sample.indicator == LETTERS[sample])*2
            parms. <- parms[which.columns/2]
            
            
            ## expected reference read count
            if(sum(s.[it,seq(1,ncol(s.),2)])==0){ 
              res.mle.ref.[which.columns/2] <- rep(0, length(which.columns/2))
            }else if(sum(parms.*s.[it,which.columns-1])==0){
              res.mle.ref.[which.columns/2] <- rep(0, length(which.columns/2))
            }else{
              
              tmp <- parms.*(s.[it,which.columns-1] * p.[it, which.columns]/
                               (sum(parms.*s.[it, which.columns-1] * p.[it, which.columns]) + (1-sum(parms.))*ifelse(locus[i,1] %in% c("X", "Y") & sex=="male",1,2) )) * 
                read.count.matrix.ref[i,sample]
              
              res.mle.ref.[which.columns/2] <- tmp
              
              
            }
            
            
            
            ## expected mutant read count
            if(sum(parms.*s.[it,which.columns])==0){
              res.mle.mut.[which.columns/2] <- rep(0, length(which.columns/2))
            }else{
              
              tmp <- parms. * (s.[it,which.columns] * p.[it, which.columns]/
                                 (sum(parms.*s.[it, which.columns] * p.[it, which.columns]) )) *
                read.count.matrix.mut[i,sample]
              
              res.mle.mut.[which.columns/2] <- tmp
              
              
            }
            
            
            
            if(sum(res.mle.mut.[which.columns/2])==0 & read.count.matrix.mut[i,sample]>0){ ## set -logL to infinity if no mutated reads are expected (since we have them!)
              LL.tmp.2 <- Inf
              
            }else if(sum(parms.*s.[it,which.columns])==0 & read.count.matrix.mut[i,sample]>0){
              LL.tmp.2 <- Inf
              
            }else{
              
              res.mle. <- cbind(res.mle.ref., res.mle.mut.)
              res.mle. <- as.vector(t(res.mle.))
              res.mle. <- round(res.mle.)
              
              
              normal.read.counts <- sum(read.count.matrix.ref[i,sample]) - round(sum(res.mle.ref.[which.columns/2]))
              
              normal.tissue.[sample] <- normal.read.counts
              
              ## compute the negative log-likelihood of this solution
              LL.tmp.2 <- LL.tmp.2 + 
                NegLogLikelihood(readcounts = matrix(c(res.mle.[as.vector(rbind(which.columns-1, which.columns))], normal.read.counts, 0), nrow=1), 
                                 coverage.ratio = c(s.[it,as.vector(rbind(which.columns-1, which.columns))],
                                                    rep(ifelse(locus[1]%in%c("X","Y") & sex=="male",1,0),2))* 
                                   matrix(c(p.[it,as.vector(rbind(which.columns-1, which.columns))], rep(ifelse(locus[1]%in%c("X","Y") & sex=="male",1,2)),2)/ 
                                            sum(sum(parms. * (p.[it,seq(1, ncol(p.), 2)])[which.columns/2]) +
                                                  (1-sum(parms.))*ifelse(locus[1]%in%c("X","Y") & sex=="male",1,2)), nrow=1),
                                 parms = rep(c(parms., 1-sum(parms.)), each=2), multiple.samples = T)
              
              
              
              if(read.count.matrix.mut[i,sample]==0){
                
                
                if(NegLogLikelihood(readcounts = matrix(c(res.mle.[as.vector(rbind(which.columns-1, which.columns))], normal.read.counts, 0), nrow=1), 
                                    coverage.ratio = c(s.[it,as.vector(rbind(which.columns-1, which.columns))],
                                                       rep(ifelse(locus[1]%in%c("X","Y") & sex=="male",1,0),2))* 
                                    matrix(c(p.[it,as.vector(rbind(which.columns-1, which.columns))], rep(ifelse(locus[1]%in%c("X","Y") & sex=="male",1,2)),2)/ 
                                           sum(sum(parms. * (p.[it,seq(1, ncol(p.), 2)])[which.columns/2]) +
                                               (1-sum(parms.))*ifelse(locus[1]%in%c("X","Y") & sex=="male",1,2)), nrow=1),
                                    parms = rep(c(parms., 1-sum(parms.)), each=2), multiple.samples = T) > 5){
                  
                  
                }
              }
            }
          }
          
          
          
          
          if(it==1){
            LL.tmp <- LL.tmp.2
            
            s[k,as.vector(rbind(2*which.clones-1, 2*which.clones))] <- s.[it,as.vector(rbind(2*which.clones-1,2* which.clones))]
            res.mle.mut[which(sample.indicator %in% sample.indicator[which.clones])] <- res.mle.mut.[which(sample.indicator %in% sample.indicator[which.clones])]
            res.mle.ref[which(sample.indicator %in% sample.indicator[which.clones])] <- res.mle.ref.[which(sample.indicator %in% sample.indicator[which.clones])]
            
            
            normal.tissue.tmp[which(LETTERS %in% sample.indicator[which.clones])] <- normal.tissue.[which(LETTERS %in% sample.indicator[which.clones])]
            
          }
          
          ## update in case we got a better result
          if(LL.tmp.2 < LL.tmp){
            LL.tmp <- LL.tmp.2
            s[k,as.vector(rbind(2*which.clones-1, 2*which.clones))] <- s.[it,as.vector(rbind(2*which.clones-1, 2*which.clones))]
            res.mle.mut[which(sample.indicator %in% sample.indicator[which.clones])] <- res.mle.mut.[which(sample.indicator %in% sample.indicator[which.clones])]
            res.mle.ref[which(sample.indicator %in% sample.indicator[which.clones])] <- res.mle.ref.[which(sample.indicator %in% sample.indicator[which.clones])]
            normal.tissue.tmp[which(LETTERS %in% sample.indicator[which.clones])] <- normal.tissue.[which(LETTERS %in% sample.indicator[which.clones])]
          }
          
          
        }
        LL.at.k <- LL.at.k + LL.tmp
        
      }
      
      
      
      
      if(likelihood.space){
        LL.all.sol[k] <- LL.at.k
      }
      
      
      res.mle <- cbind(res.mle.ref, res.mle.mut)
      res.mle <- as.vector(t(res.mle))
      res.mle <- round(res.mle)
      
      
      
      
      ## initialize LL in the first iteration, res.most.likely, s.most.likely
      if(k==1){
        LL <- LL.at.k 
        s.most.likely <- s[k,]
        cn.most.likely <- p[k,]
        cn.most.likely.index <- cn.index[k,]
        res.most.likely <- res.mle
        normal.tissue[i,] <- normal.tissue.tmp
      }
      
      
      
      ## if the current result is more likely than any previous ones, replace res.most.likely with this result
      
      if(LL.at.k<LL & !any(res.mle<0)){
        res.most.likely <- res.mle
        s.most.likely <- s[k,]
        cn.most.likely <- p[k,]
        cn.most.likely.index <- cn.index[k,]
        LL <- LL.at.k
        normal.tissue[i,] <- normal.tissue.tmp
      }
      
    }
    
    
    ## store the most likely result for this gene in the result table that summarizes the most likely results for all genes
    
    res[i,] <- res.most.likely
    s.res[i,] <- s.most.likely
    cn.res[i,] <- cn.most.likely
    cn.res.index[i,] <- cn.most.likely.index
    
    
    ## if we are interested in a likelihood space (all solutions within the 95 % LL), then add additional solutions and store their likelihood
    if(likelihood.space){
      ## the total likelihood is so small that we cannot estimate a likelihood per solution
      
      like.per.loc[i] <- exp(-LL)/(sum(exp(-LL.all.sol)))
      
      if(sum(exp(-LL.all.sol))==0){
        like.per.loc[i] <- NA
      }
      
      
      additional.solutions <- min(which(log(cumsum(sort(exp(-LL.all.sol), decreasing=T))) - log(sum(exp(-LL.all.sol))) > log(0.9)))
      
      
      if(length(additional.solutions)> 0){
        if(min(additional.solutions)==1 | is.na(additional.solutions) | additional.solutions==Inf){
          additional.solutions <- c()
        }
        
      }
      additional.solutions <- order(exp(-LL.all.sol.prim), decreasing=T)[additional.solutions]
      
      
      
      
      if(length(additional.solutions)>0){
        like.per.loc <- c(like.per.loc, exp((-LL.all.sol)[additional.solutions] - log(sum(exp(-LL.all.sol)))))
        if(sum(exp(-LL.all.sol))==0){
          like.per.loc[length(like.per.loc)] <- NA
        }
        s.res <- rbind(s.res, s[additional.solutions,])
        rownames(s.res) <- c(rownames(s.res)[-((nrow(s.res)-length(additional.solutions)+1):nrow(s.res))], rep(rownames(s.res)[i], length(additional.solutions)))
        cn.res <- rbind(cn.res, p[additional.solutions,])
        rownames(cn.res) <- rownames(s.res)
        
        cn.res.index <- rbind(cn.res.index, cn.index[additional.solutions,])
        rownames(cn.res.index) <- c(rownames(cn.res.index)[-((nrow(cn.res.index)-length(additional.solutions)+1):nrow(cn.res.index))], rep(rownames(cn.res.index)[i], length(additional.solutions)))
        
      }
    }
    
    
    
  }
  
  
  
  ## if there are any clones wo alterations, add them to normal tissue
  for(l in LETTERS[1:ncol(coverage.ratio)]){
    which.columns <- which(sample.indicator == l)*2
    s.res. <- s.res[,as.vector(rbind(which.columns-1, which.columns)),drop=F]
    cn.res.index. <- cn.res.index[,as.vector(rbind(which.columns-1, which.columns)),drop=F]
    res. <- res[,as.vector(rbind(which.columns-1, which.columns)),drop=F]
    if(any(colSums(s.res.[,seq(2, ncol(s.res.),2),drop=F])==0 & colSums(cn.res.index.[,seq(2, ncol(s.res.),2),drop=F])==0)){
      normal.clone <- which(colSums(s.res.[,seq(2, ncol(s.res.),2),drop=F])==0 & colSums(cn.res.index.[,seq(2, ncol(s.res.),2),drop=F])==0)
      for(i in normal.clone){
        res.[,c(i*2-1,i*2)] <- matrix(0, nrow(res.), 2)
        res[,which.columns] <- res.
      }
    }
  }
  
  
  if(zero){
    
    for(i in zeros){
      if(i==1){
        res <- cbind(matrix(0, nrow=nrow(res), ncol=2), res)
        s.res <- cbind(matrix(c(1, 0), byrow=T, nrow=nrow(res), ncol=2), s.res)
        cn.res <- cbind(matrix(normal.copy.number.auto, nrow=nrow(res), ncol=2), cn.res)
      }else if(i>=(ncol(res)/2)){
        res <- cbind(res, matrix(0, nrow=nrow(res), ncol=2))
        s.res <- cbind(s.res, matrix(c(1,0), byrow=T, nrow=nrow(res), ncol=2))
        cn.res <- cbind(cn.res, matrix(normal.copy.number.auto, nrow=nrow(res), ncol=2))
      }else{
        res <- cbind(res[,1:(2*(i-1))], matrix(0, nrow=nrow(res), ncol=2), res[,(2*i -1 ):ncol(res)])
        s.res <- cbind(s.res[,1:(2*(i-1))], matrix(c(1,0), byrow=T, nrow=nrow(res), ncol=2), s.res[,(2*i -1 ):ncol(s.res)])
        cn.res <- cbind(cn.res[,1:(2*(i-1))], matrix(normal.copy.number.auto, nrow=nrow(res), ncol=2), cn.res[,(2*i -1):ncol(cn.res)])
      }
    }
  }
  
  res.list <- list(readcounts=res, norm.copy.number=s.res, copy.number = cn.res, copy.number.index = cn.res.index, 
                   break.function = break.function, number.of.possibilities = number.of.possibilities, normal.tissue = normal.tissue)
  
  if(likelihood.space){
    names(like.per.loc) <- rownames(cn.res.index)
    res.list <- list(readcounts=res, norm.copy.number=s.res, copy.number = cn.res, 
                     copy.number.index = cn.res.index,  break.function = break.function, like.per.loc = like.per.loc, 
                     number.of.possibilities = number.of.possibilities, cn.indicator = cn.indicator.all, normal.tissue = normal.tissue)
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
#' @param tree a binary tree matrix; each column corresponds to one subclone; each row to possible relationships between mutations (1: mutation is present in clone, 0: mutation is absent in clone)
#' @param read.count.matrix.mut a matrix storing the number of mutant reads in each sample (row: location, column: sample)
#' @param read.count.matrix.ref a matrix storing the number of reference reads in each sample (row: location, column: sample)
#' @param coverage.ratio a matrix storing the measured coverage ratioper sample
#' @param baf a matrix storing the measured BAF per sample
#' @param locus matrix containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param parms vector containing the relative subclone sizes of each clone, corresponding to the clades in the tree
#' @param sample.indicator vector assigning each clade in the tree to a sample (e.g. c("A", "B", "C", "A", "B) for a tree with 5 clades)
#' @param normal.copy.number.auto normal copy number on autosomes 
#' @param normal.copy.number.y normal copy number on the y chromosome
#' @param normal.b.allele.norm normal B allele count
#' @param joined.info List of lists containing prior information on copy number changes that occurred together. Length of list must correspond to number of loci. For example, if it is known that the copy number change at locus x is present in sample 'A' and 'B' set joined.info[[x]] <- c("A", "B). If there's no information, set joined.info[[x]] <- NA
#' @param false.negatives logical; should the possibility of false negative mutations be accounted for?
#' @param cutoff cutoff for EM algorithm (RSS between two subsequent iterations)
#' @param maxit maximal number of iterations
#' @return The relative subclone sizes of primary and relapse tumor (last value corresponds to normal tissue).
#' @export

EM.alg <- function(tree = "classic", read.count.matrix.ref, read.count.matrix.mut, parms, coverage.ratio, locus,
                   cutoff, maxit = 25, sex, sample.indicator, baf = baf, normal.copy.number.auto = 2, normal.copy.number.y = 1, normal.b.allele.norm = 1, 
                   joined.info = NA, extra.muts = NULL, false.negatives = F){
  
  
  ## initialization
  ## start with expectation
  prev.it <- Expectation(tree = tree, read.count.matrix.mut = read.count.matrix.mut, read.count.matrix.ref = read.count.matrix.ref, 
                         parms = parms, coverage.ratio = coverage.ratio, baf = baf, locus = locus, likelihood.space = F, sample.indicator = sample.indicator,
                         normal.copy.number.auto = 2, normal.copy.number.y = 1, sex = sex, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info,
                         false.negatives = false.negatives)
  
  
  if(is.na(prev.it[[1]][1])){return(NA)}
  if(prev.it[[1]][1]=="Impurity"){return(NA)}
  if(prev.it[[1]][1]=="Impurity"){return(NA)}
  
  if(prev.it$break.function==T){
    print("break.function = T")
    return(NA)}
  
  
  ## corresponding negative log likelihood
  prev.ll <- 0
  new.parms <- rep(0, length(parms))
  for(i in 1:ncol(coverage.ratio)){
    which.columns <- which(sample.indicator == LETTERS[i])*2
    which.columns. <- as.vector(rbind(which.columns-1, which.columns))
    
    copy.number. <- prev.it$copy.number[,which.columns,drop=F]
    copy.number. <- cbind(copy.number., 2)
    
    
    readcounts <- prev.it$readcounts[,which.columns.]
    
    ## add the normal clone
    readcounts <- cbind(readcounts, prev.it$normal.tissue[,i])
    readcounts <- cbind(readcounts, 0)
    
    copy.number. <- prev.it$copy.number[,which.columns.]
    copy.number. <- cbind(copy.number., 2)
    copy.number. <- cbind(copy.number., 2)
    if(sex=="male"){
      copy.number.[locus[,1]%in%c("X", "Y"),c(ncol(copy.number.)-1, ncol(copy.number.))] <- 1
    }
    norm.copy.number. <- prev.it$norm.copy.number[,which.columns.]
    norm.copy.number. <- cbind(norm.copy.number.,1)
    norm.copy.number. <- cbind(norm.copy.number.,1)
    
    prev.ll <- prev.ll + NegLogLikelihood(readcounts=readcounts, copy.number=norm.copy.number.*copy.number./
                                            colSums(c(parms[which.columns/2], 1-sum(parms[which.columns/2]))*t(copy.number.)), parms = rep(c(parms[which.columns/2], 1- sum(parms[which.columns/2])), each=2))
    
    
    ## maximization
    
    ## only maximize based on readcounts that have non-zero mutation counts
    readcounts <- readcounts[rowSums(readcounts[,seq(2, ncol(readcounts),2),drop=F])!=0 & !rownames(readcounts) %in% extra.muts,,drop=F]
    copy.number. <- copy.number.[rownames(readcounts),,drop=F]
    
    
    new.parms[which.columns/2] <- Maximization(readcounts = readcounts, copy.number = copy.number., 
                                               n.clones = ncol(readcounts)/2-1)
    
    
  }
  
  rm(prev.it)
  gc()
  
  prev.it <- new.parms
  
  ## iterate
  it.counter <- 0
  while(it.counter < maxit){
    
    ## set clones < 3 % to zero and remove them from the tree
    if(any(prev.it < 0.03)){
      prev.it[prev.it < 0.03] <- 0
      
    }
    
    it.clus <- Expectation(tree = tree, read.count.matrix.mut = read.count.matrix.mut, read.count.matrix.ref = read.count.matrix.ref, 
                           parms = prev.it, coverage.ratio = coverage.ratio, baf = baf, locus = locus, likelihood.space = F, sample.indicator = sample.indicator,
                           normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, sex = sex, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info,
                           false.negatives = false.negatives)
    
    
    if(is.na(it.clus[[1]][1])){return(NA)}
    
    if(it.clus[[1]][1]=="Impurity"){return(NA)}
    
    if(it.clus$break.function==T){
      print("break.function = T")
      return(NA)}
    
    ll <- 0
    new.parms <- rep(0, length(parms))
    for(i in 1:ncol(coverage.ratio)){
      which.columns <- which(sample.indicator == LETTERS[i])*2
      which.columns. <- as.vector(rbind(which.columns-1, which.columns))
      
      
      
      readcounts <- it.clus$readcounts[,which.columns.]
      ## add the normal clone
      readcounts <- cbind(readcounts, it.clus$normal.tissue[,i])
      readcounts <- cbind(readcounts, 0)
      
      copy.number. <- it.clus$copy.number[,which.columns.,drop=F]
      copy.number. <- cbind(copy.number., 2)
      copy.number. <- cbind(copy.number., 2)
      if(sex=="male"){
        copy.number.[locus[,1]%in%c("X", "Y"),c(ncol(copy.number.)-1, ncol(copy.number.))] <- 1
      }
      ## maximization
      
      ## only maximize based on readcounts that have non-zero mutation counts
      readcounts.for.max <- readcounts[rowSums(readcounts[,seq(2, ncol(readcounts),2),drop=F])!=0 & !rownames(readcounts) %in% extra.muts,,drop=F]
      copy.number.for.max <- copy.number.[rownames(readcounts.for.max),,drop=F]
      
      new.parms[which.columns/2] <- Maximization(readcounts = readcounts.for.max, coverage.ratio = copy.number.for.max, 
                                                 n.clones = ncol(readcounts)/2-1)
      
      
      
      norm.copy.number. <- it.clus$norm.copy.number[,which.columns.]
      norm.copy.number. <- cbind(norm.copy.number.,1)
      norm.copy.number. <- cbind(norm.copy.number.,1)
      
      
      
      ll <- ll + NegLogLikelihood(readcounts=readcounts, copy.number=norm.copy.number.*copy.number./
                                    colSums(c(parms[which.columns/2], 1-sum(parms[which.columns/2]))*t(copy.number.)), parms = rep(c(parms[which.columns/2], 1-sum(parms[which.columns/2])), each=2))
      
    }
    it <- new.parms
    
    if(sum(it[-length(it)]==0)==(length(it)-1)){return(NA)}
    
    
    
    rm(it.clus)
    gc()
    ## comparison between parameter sets is only possible for equally sized data frames. If we remove a cluster due to a mean of zero, we cannot compare the mean
    ## values to the former parameter set and thus, allow one more iteration before we do a comparison again
    if(length(it)==length(c(prev.it))){
      if(sum((it-prev.it)^2)<cutoff){break}
      
      if(ll > prev.ll){
        it <- prev.it
        break}
    }
    
    #if(ll > prev.ll){break}
    prev.it <- it
    
    prev.ll <- ll
    
    it.counter <- it.counter + 1
  }
  
  return(list(res=it))
  
  
  
}


#' Repeats the EM algorithm often with random starting parameters to find the global optimum (EM itself does only local optimization)
#' 
#' All possible trees for up to three primary and relapse subclones are automatically tested. More complex trees can be tested as well, but information must be provided on how these trees should look like
#' @param tree.list A list of candidate tree matrices (binary). each tree matrix conssists of columns corresponding to subclones; each row to possible relationships between mutations (1: mutation is present in clone, 0: mutation is absent in clone)
#' @param read.count.matrix.mut a matrix storing the number of mutant reads in each sample (row: location, column: sample)
#' @param read.count.matrix.ref a matrix storing the number of reference reads in each sample (row: location, column: sample)
#' @param coverage.ratio a matrix storing the measured coverage ratioper sample
#' @param baf a matrix storing the measured BAF per sample
#' @param locus matrix containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param cutoff the cutoff for the EM algorithm
#' @param it the number of iterations per candidate tree (random starting conditions)
#' @param sample.indicator.list list of indicator vectors; each vector indicates which subclone in the corresponding tree belongs to which sample
#' @param normal.copy.number.auto normal copy number on autosomes 
#' @param normal.copy.number.y normal copy number on the y chromosome
#' @param normal.b.allele.norm normal B allele count
#' @param trees.from.muts which trees are built based on the mutation structure in the tree and are sufficiently complex that they don't need to  be extended by the algorithm
#' @param joined.info List of lists containing prior information on copy number changes that occurred together. Length of list must correspond to number of loci. For example, if it is known that the copy number change at locus x is present in sample 'A' and 'B' set joined.info[[x]] <- c("A", "B). If there's no information, set joined.info[[x]] <- NA
#' @param extra.muts
#' @param false.negatives logical; should the possibility of false negative mutations be accounted for?
#' @param directory directory for output file (not saving directory)
#' @return A list containing for each candidate tree matrices which contain the relative subclone sizes for the primary and the relapse tumor for each fit (rows correspond to fits) and a vector containing the corresponding likelihood. In the case of more complex trees, an additional list element reports the candidate tree and vectors indicating which subclones are primary clones and which are relapse clones.
#' @export


model.fits <- function(tree.list, readcounts.mut, readcounts.ref, coverage.ratio, locus, cutoff=5e-4, it = 100, sex = "female",
                       directory = "", baf, sample.indicator.list, normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y,
                       normal.b.allele.norm = normal.b.allele.norm, trees.from.muts = 0, joined.info = NA, extra.muts = NULL, false.negatives = F,
                       clusterComputation = F, n.cores = 8){
  
  resultlist <- list(0)
  
  
  if(clusterComputation){
    library(doParallel)
    cl<-makeCluster(n.cores, type = "FORK", outfile=paste0(directory, "/Log.txt"))
    registerDoParallel(cores=n.cores)
  }
  
  
  ## test several phylogenetic trees and extensions thereof, stored in the tree list
  
  ii.start <- 1
  
  for(ii in ii.start:length(tree.list)){
    
    n.clones <- length(sample.indicator.list[[ii]])
    sim.res <- c()
    sim.like <- c()
    
    print(paste(ii, "out of ", length(tree.list), sep=" "))
    
    tree <- tree.list[[ii]]
    tree <- unique(tree)
    sample.indicator <- sample.indicator.list[[ii]]
    
    
    fits <- foreach(i=1:it, .combine=rbind) %do%{
      set.seed(i)
      
      parms <- rep(0, length(sample.indicator))
      
      for(j in 1:ncol(coverage.ratio)){
        parms. <- Rand.parms(sum(sample.indicator==LETTERS[j])+1)
        parms[which(sample.indicator==LETTERS[j])] <- parms.[-1]
      }
      
      res <- EM.alg(tree = tree, read.count.matrix.ref = readcounts.ref, read.count.matrix.mut = readcounts.mut, parms = parms, 
                    coverage.ratio = coverage.ratio, locus = locus, cutoff = cutoff, maxit = 50, sex = sex, sample.indicator = sample.indicator, baf = baf,
                    normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info, extra.muts = extra.muts,
                    false.negatives = false.negatives)
      
      if(!is.na(res[[1]])){
        
        parms <- res$res
        
        res <- EM.alg(tree = tree, read.count.matrix.ref = readcounts.ref, read.count.matrix.mut = readcounts.mut, parms = parms, 
                      coverage.ratio = coverage.ratio, locus = locus, cutoff = cutoff, maxit = 50, sex = sex, sample.indicator = sample.indicator, baf = baf,
                      normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info, extra.muts = extra.muts,
                      false.negatives = false.negatives)
      }
      
      if(!is.na(res[[1]])){
        
        if(res[[1]][1]!="Impurity"){
          
          res.clusters <- Expectation(tree = tree, read.count.matrix.mut = readcounts.mut, read.count.matrix.ref = readcounts.ref,
                                      parms = res$res, coverage.ratio = coverage.ratio, baf = baf, locus = locus, sex = sex, sample.indicator = sample.indicator,
                                      normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info,
                                      false.negatives = false.negatives)
          
          
          if(res.clusters[[1]][1]=="Impurity"){
            return.value <- c(rep(0, length(parms)), Inf)
          }else
            if( !any(res$res<0) & !any(res.clusters$readcounts<0)){
              
              
              like <- 0
              ## iterate through the samples and compute the likelihood
              for(jj in 1:ncol(coverage.ratio)){
                which.columns <- which(sample.indicator == LETTERS[jj])*2
                which.columns. <- as.vector(rbind(which.columns-1, which.columns))
                coverage.ratio. <- res.clusters$copy.number[,which.columns.,drop=F]
                coverage.ratio. <- cbind(coverage.ratio., 2)
                coverage.ratio. <- cbind(coverage.ratio., 2)
                if(sex=="male"){
                  coverage.ratio.[which(locus[,1] %in% c("X", "Y")),c(ncol(coverage.ratio.), ncol(coverage.ratio.)-1)] <- 1
                }
                readcounts <- res.clusters$readcounts[,which.columns.]
                
                ## add the normal clone
                readcounts <- cbind(readcounts, res.clusters$normal.tissue[,jj])
                readcounts <- cbind(readcounts, 0)
                
                
                norm.copy.number. <- res.clusters$norm.copy.number[,which.columns.]
                norm.copy.number. <- cbind(norm.copy.number.,1)
                norm.copy.number. <- cbind(norm.copy.number.,0)
                
                
                like <- like + NegLogLikelihood(readcounts=readcounts, copy.number=coverage.ratio.*norm.copy.number./
                                                  colSums(c(res$res[which.columns/2], 1-sum(res$res[which.columns/2]))*t(coverage.ratio.)), 
                                                parms = rep(c(res$res[which.columns/2],1-sum(res$res[which.columns/2])), each=2))
                
              }
              
              
              return.value <- c(res$res, like)
              
            }else{
              return.value <- c(rep(0, length(parms)), Inf)
            }
          rm(res.clusters)
          gc()
          
        }else{
          return.value <- c(rep(0, length(parms)), Inf)
        }
      }else{
        return.value <- c(rep(0, length(parms) ), Inf)
      }
      
      return.value
      
    }
    
    
    if(length(fits)>0){
      for(i.fits in 1:nrow(fits)){
        if(!is.infinite(fits[i.fits, ncol(fits)])){
          sim.like <- c(sim.like, fits[i.fits, ncol(fits)])
          sim.res <- rbind(sim.res, fits[i.fits, 1:(ncol(fits)-1)])
        }
      }
    }
    
    ## compute the expectation at the current best fit
    res.clusters <- Expectation(tree = tree, read.count.matrix.mut = readcounts.mut, read.count.matrix.ref = readcounts.ref,
                                parms = fits[which.min(fits[,ncol(fits)]),-ncol(fits)], 
                                coverage.ratio = coverage.ratio, baf = baf, locus = locus, sex = sex, 
                                sample.indicator = sample.indicator,
                                normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, joined.info = joined.info,
                                false.negatives = false.negatives)
    
    
    ## compute the adjusted BIC
    BIC.old <- BIC(-min(fits[,ncol(fits)]), n.parms = ncol(tree), n.data = sum(readcounts.mut + readcounts.ref),
                   Tau = res.clusters$number.of.possibilities, gamma = 0.9)
    
    resultlist[[length(resultlist)+1]] <- list(Parameters=sim.res,  Likelihood=sim.like, Tree = tree, indicator = sample.indicator)
    save(resultlist, file=paste0(directory,"/resultlist.RData"))
    
    ## try to extend the tree and see whether this gives a better fit
    extended.indicator <- sample.indicator.list[[ii]]
    extended.tree <- tree
    
    BIC.new <- BIC.old
    
    
    
    for(hh in 1:ncol(coverage.ratio)){
      print(paste("Sample Nr", hh))
      nodes <- which(extended.indicator==LETTERS[hh])
      if(ii %in% trees.from.muts){
        next
        nodes <- 1:length(extended.indicator)
      }
      
      ## if this sample already has three nodes, don't try to add more
      if(length(nodes)<3 | ii %in% trees.from.muts){
        try=3-length(nodes)
      }else{
        try=3
      }
      while(try < 3){
        try <- try + 1
        if(length(which(extended.indicator==LETTERS[hh]))>=3){
          next
        }
        ## try to add clones upstream of each node
        for(k in nrow(extended.tree):2){
          BIC.new. <- Inf
          ## which subclones arise from this node?
          subclones.from.node <- which(extended.tree[k,]==1)
          extended.tree. <- rbind(extended.tree, extended.tree[k,])
          extended.tree. <- rbind(extended.tree., 0)
          mutated.nodes.in.new.tree <-  c(which(rowSums(extended.tree.[,subclones.from.node,drop=F])==length(subclones.from.node)),
                                          nrow(extended.tree.))
          mutated.nodes.in.new.tree <- mutated.nodes.in.new.tree[-(length(mutated.nodes.in.new.tree)-1)]
          extended.tree. <- cbind(extended.tree., replace(rep(0, nrow(extended.tree.)), 
                                                          mutated.nodes.in.new.tree,1))
          colnames(extended.tree.)[ncol(extended.tree.)] <- LETTERS[hh]
          colnames(extended.tree.) <- make.unique(colnames(extended.tree.))
          extended.tree. <- extended.tree.[rowSums(extended.tree.)>0,,drop=F]
          extended.indicator. <- c(extended.indicator, LETTERS[hh])
          sim.res <- c()
          sim.like <- c()
          
          ## test the new tree
          fits <- foreach(i=1:100, .combine=rbind) %dopar%{
            set.seed(i)
            
            parms <- rep(0, length(extended.indicator.))
            
            for(j in 1:ncol(coverage.ratio)){
              parms. <- Rand.parms(sum(extended.indicator.==LETTERS[j])+1)
              parms[which(extended.indicator.==LETTERS[j])] <- parms.[-1]
            }
            
            
            res <- EM.alg(tree = extended.tree., read.count.matrix.ref = readcounts.ref, read.count.matrix.mut = readcounts.mut, parms = parms,
                          coverage.ratio = coverage.ratio, locus = locus, cutoff = cutoff, maxit = 50, sex = sex, sample.indicator = extended.indicator., baf = baf,
                          normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, 
                          joined.info = joined.info, extra.muts = extra.muts, false.negatives = false.negatives)
            
            if(!is.na(res[[1]])){
              parms <- res$res
              
              res <- EM.alg(tree = extended.tree., read.count.matrix.ref = readcounts.ref, read.count.matrix.mut = readcounts.mut, parms = parms,
                            coverage.ratio = coverage.ratio, locus = locus, cutoff = cutoff, maxit = 50, sex = sex, sample.indicator = extended.indicator., baf = baf,
                            normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, 
                            joined.info = joined.info, extra.muts = extra.muts, false.negatives = false.negatives)
            }
            
            if(!is.na(res[[1]])){
              if(res[[1]][1]!="Impurity"){
                
                res.clusters <- Expectation(tree = extended.tree., read.count.matrix.mut = readcounts.mut, read.count.matrix.ref = readcounts.ref,
                                            parms = res$res, coverage.ratio = coverage.ratio, baf = baf, locus = locus, sex = sex, sample.indicator = extended.indicator.,
                                            normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, 
                                            joined.info = joined.info, false.negatives = false.negatives)
                
                
                
                if(!any(res$res<0) & !any(res.clusters$readcounts<0)){
                  
                  
                  like <- 0
                  for(jj in 1:ncol(coverage.ratio)){
                    which.columns <- which(extended.indicator. == LETTERS[jj])*2
                    which.columns. <- as.vector(rbind(which.columns-1, which.columns))
                    coverage.ratio. <- res.clusters$copy.number[,which.columns.,drop=F]
                    coverage.ratio. <- cbind(coverage.ratio., 2)
                    coverage.ratio. <- cbind(coverage.ratio., 2)
                    if(sex=="male"){
                      coverage.ratio.[which(locus[,1] %in% c("X", "Y")),c(ncol(coverage.ratio.), ncol(coverage.ratio.)-1)] <- 1
                    }
                    readcounts <- res.clusters$readcounts[,which.columns.]
                    
                    ## add the normal clone
                    readcounts <- cbind(readcounts, readcounts.ref[,jj] - rowSums(readcounts[,seq(1,ncol(readcounts),2),drop=F]))
                    readcounts <- cbind(readcounts, 0)
                    
                    
                    norm.copy.number. <- res.clusters$norm.copy.number[,which.columns.]
                    norm.copy.number. <- cbind(norm.copy.number.,1)
                    norm.copy.number. <- cbind(norm.copy.number.,1)
                    
                    
                    like <- like + NegLogLikelihood(readcounts=readcounts, copy.number=coverage.ratio.*norm.copy.number./
                                                      colSums(c(res$res[which.columns/2], 1-sum(res$res[which.columns/2]))*t(coverage.ratio.)),
                                                    parms = rep(c(res$res[which.columns/2],1-sum(res$res[which.columns/2])), each=2))
                    
                  }
                  
                  
                  return.value <- c(res$res, like)
                  
                }else{
                  return.value <- c(rep(0, length(parms)), Inf)
                }
                rm(res.clusters)
                gc()
                
              }else{
                return.value <- c(rep(0, length(parms)), Inf)
              }
            }else{
              return.value <- c(rep(0, length(parms) ), Inf)
            }
            
            return.value
            
          }
          
          
          if(length(fits)>0){
            for(i.fits in 1:nrow(fits)){
              if(!is.infinite(fits[i.fits, ncol(fits)])){
                sim.like <- c(sim.like, fits[i.fits, ncol(fits)])
                sim.res <- rbind(sim.res, fits[i.fits, 1:(ncol(fits)-1)])
              }
            }
            resultlist[[length(resultlist)+1]] <- list(Parameters=sim.res,  Likelihood=sim.like, Tree = extended.tree., indicator = extended.indicator.)
            save(resultlist, file=paste0(directory,"/resultlist.RData"))
          }
          
          if(length(fits)>0){
            res.clusters <- Expectation(tree = extended.tree., read.count.matrix.mut = readcounts.mut, read.count.matrix.ref = readcounts.ref,
                                        parms = fits[which.min(fits[,ncol(fits)]),-ncol(fits)],
                                        coverage.ratio = coverage.ratio, baf = baf, locus = locus, sex = sex,
                                        sample.indicator = extended.indicator.,
                                        normal.copy.number.auto = normal.copy.number.auto, normal.copy.number.y = normal.copy.number.y, normal.b.allele.norm = normal.b.allele.norm, 
                                        joined.info = joined.info, false.negatives = false.negatives)
            
            if(!is.na(res.clusters[[1]])){
              BIC.new. <- BIC(-min(fits[,ncol(fits)]), n.parms = ncol(extended.tree.), n.data = sum(readcounts.mut + readcounts.ref), 
                              Tau = res.clusters$number.of.possibilities, gamma = 0.9)
            }else{
              BIC.new. <- Inf
            }
            
            
          }else{
            BIC.new. <- Inf
          }
          
          
          
          
          if(BIC.new. < BIC.new){
            extended.tree <- extended.tree.
            extended.indicator <- extended.indicator.
            BIC.new <- BIC.new.
          }
          
        }
        
        
      }
      
      
    }
    
    
    
    
    
  }
  
  
  
  
  stopImplicitCluster()
  stopCluster(cl)
  
  return(resultlist)
}



#' Function to generate all hypothetical distributions of a mutation among subclones.
#' 
#' This function generates matrices which indicate the combinations of subclones in a tree that can share a mutation and the number of alleles on which this mutation is present.
#' @param tree a binary tree matrix; each column corresponds to one subclone; each row to possible relationships between mutations (1: mutation is present in clone, 0: mutation is absent in clone)
#' @param copy.number Inferred copy number per sample
#' @param readcounts.mut Vector containing measured numbers of mutated read counts per sample
#' @param best.cn.sol.all binary matrix indicating whether a subclone carries a CNV or not
#' @param locus matrix containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param sample.identifier vector assigning each clade in the tree to a sample (e.g. c("A", "B", "C", "A", "B) for a tree with 5 clades)
#' @param normal.copy.number.auto normal copy number on autosomes 
#' @param normal.copy.number.y normal copy number on the y chromosome
#' @param false.negatives logical, whether the possbility of false negative calls should be accounted for
#' @return Two matrices containing the different possibilities for mutated alleles across the subclones as well as the possible copy numbers
#' @export
#

Index.Matrix <- function(tree, copy.number,  sample.identifier, readcounts.mut, locus, sex, best.cn.sol.all = NA, normal.copy.number.auto = 2,
                         normal.copy.number.y = 1, false.negatives = F){
  
  ## indicator matrix:which clones harbor an SNV?
  s <- tree
  ## copy number matrix: which copy number in which clone?
  p <- s
  
  mutated.samples <- names(readcounts.mut)[readcounts.mut!=0]
  
  non.mutated.samples <- names(readcounts.mut)[!names(readcounts.mut) %in% mutated.samples]
  
  ## Case I: we are looking at extra mutations (in case non-mutated loci containing a CNV were added in order to map their CNVs at the tree)
  
  if(all(sample.identifier %in% non.mutated.samples)){
    ## no SNVs
    s <- replace(s, s!=0, 0)
    
    p <- c()
    for(i in 1:length(best.cn.sol.all)){
      
      if(!is.na(best.cn.sol.all[[i]][1])){
        tmp <- replace(best.cn.sol.all[[i]],best.cn.sol.all[[i]]!=0, copy.number[[i]])
        if(length(p)>0){
          p[p==0] <-  tmp[p==0] + p[p==0]
        }else{
          p <- tmp
        }
      }
    }
    
    if(length(p)==0){
      p <- matrix(rep(ifelse(sex=="male" & locus[1] %in% c("X", "Y"), normal.copy.number.y, normal.copy.number.auto), ncol(s)), nrow=1)
      
    }
    
    p <- replace(p, p==0, ifelse(sex=="male" & locus[1] %in% c("X", "Y"), normal.copy.number.y, normal.copy.number.auto))
    
    ## now intermingle s with all ps
    p.new <- do.call(rbind, replicate(nrow(s), p, simplify=F))
    s.new <- matrix(rep(as.vector(s), each=nrow(p)), ncol=ncol(s))
    s <- s.new
    p <- p.new
    
    
    return(list(s=s, p=p))
  }
  
  ## Case II: we are looking at mutated sites
  ## make sure that all mutated samples are mutated in s and that no non-mutated sample is mutated in s
  
  to.keep <- rep(T, nrow(s))
  for(i in mutated.samples){
    to.keep <- to.keep * (rowSums(s[,sample.identifier %in% i, drop=F])!=0)
  }
  
  ## if false negatives are not considered, require that a mutation must be absent in a non mutated sample
  if(!false.negatives){
    for(i in non.mutated.samples){
      to.keep <- to.keep * (rowSums(s[,sample.identifier %in% i, drop=F])==0)
    }
  }
  
  s <- s[as.logical(to.keep),,drop=F]
  if(length(s)==0){return(list(s=s, p=s))}
  
  ## if we have a male patient and a locus on the x chromosome, change the normal copy.number to 1
  normal.copy.number <- 2
  if(locus[1] == "X" && sex == "male"){
    normal.copy.number <- 1
  }
  
  p <- c()
  for(i in 1:length(best.cn.sol.all)){
    
    
    if(!is.na(best.cn.sol.all[[i]][1])){
      tmp <- replace(best.cn.sol.all[[i]],best.cn.sol.all[[i]]!=0, copy.number[[i]])
      if(length(p)>0){
        p[p==0] <-  tmp[p==0] + p[p==0]
      }else{
        p <- tmp
      }
      
    }
  }
  
  if(length(p)==0){
    p <- matrix(rep(ifelse(sex=="male" & locus[1] %in% c("X", "Y"), normal.copy.number.y, normal.copy.number.auto), ncol(s)), nrow=1)
    
  }
  
  p <- replace(p, p==0, ifelse(sex=="male" & locus[1] %in% c("X", "Y"), normal.copy.number.y, normal.copy.number.auto))
  
  ## now intermingle s with all ps
  p.new <- do.call(rbind, replicate(nrow(s), p, simplify=F))
  s.new <- matrix(rep(as.vector(s), each=nrow(p)), ncol=ncol(s))
  s <- s.new
  p <- p.new
  rm(s.new)
  rm(p.new)
  
  
  ## if copy.number = 0 (full deletion), remove solutions where full deletion happened in all subclones (we also found a mutation here, so it needs to exist)
  ## also remove the solutions where the mutation is present on the deleted allele only
  
  for(i in 1:length(copy.number)){
    if(length(copy.number[[i]])>0){
      if(copy.number[[i]]==0){
        s <- s[rowSums(p[,which(best.cn.sol.all[[i]]!=0),drop=F])!=sum(best.cn.sol.all[[i]]!=0), ,drop=F]
        p <- p[rowSums(p[,which(best.cn.sol.all[[i]]!=0),drop=F])!=sum(best.cn.sol.all[[i]]!=0),,drop=F]
        
        mut.only.at.del <- apply(cbind(s[,which(best.cn.sol.all[[i]]!=0),drop=F],p[,which(best.cn.sol.all[[i]]!=0),drop=F]), 1, function(x){
          ifelse(sum((x[1:ncol(s[,which(best.cn.sol.all[[i]]!=0),drop=F])]!=0)==(x[(ncol(s[,which(best.cn.sol.all[[i]]!=0),drop=F])+1):(2*ncol(s[,which(best.cn.sol.all[[i]]!=0),drop=F]))]!=0))==ncol(s[,which(best.cn.sol.all[[i]]!=0),drop=F]),F,T)
        })
        s <- s[mut.only.at.del,,drop=F]
        p <- p[mut.only.at.del,,drop=F]
        p <- p[rowSums(s[,which(best.cn.sol.all[[i]]!=0),drop=F])!=0,,drop=F]
        s <- s[rowSums(s[,which(best.cn.sol.all[[i]]!=0),drop=F])!=0,,drop=F]
      }
    }
  }
  
  
  if(length(s)==0){return(list(s=s, p=s))}
  
  ## remove duplicates
  
  if(any(rowSums(s!=0)==ncol(s))){
    p <- rbind(p, p[rowSums(s!=0)==ncol(s),,drop=F])
    s <- rbind(s, replace(s[rowSums(s!=0)==ncol(s),,drop=F], s[rowSums(s!=0)==ncol(s),,drop=F]!=0, floor(normal.copy.number.auto/2)))
    
    s <- replace(s, s > p, p[s>p])
  }
  
  
  dupl.rows <- duplicated(cbind(s, p))
  
  s <- s[!dupl.rows,,drop=F]
  p <- p[!dupl.rows,,drop=F]
  
  
  ## returns the mutations and the copy numbers The number stored in mutations refers to the number of mutated alleles. 
  
  return(list(s = s, p = p))
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
#' This function infers the most likely combination of CNVs  based on a given tree structure and measured copy number ratios / B allele frequencies. Only one CNV per tumor sample can be inferred.
#' @param tree a binary tree matrix; each column corresponds to one subclone; each row to possible relationships between mutations (1: mutation is present in clone, 0: mutation is absent in clone)
#' @param coverage.ratios measured coverage ratio per sample
#' @param bafs measured B-allele frequencies per sample
#' @param locus matrix containing the chromosome and the coordinate of the mutation 
#' @param sex the gender of the patient ("male" / "female")
#' @param parms vector containing the relative subclone sizes of the individual clades in the phylogenetic tree
#' @param sample.indicator vector assigning each clade in the tree to a sample (e.g. c("A", "B", "C", "A", "B) for a tree with 5 clades)
#' @param normal.copy.number.auto normal copy number on autosomes 
#' @param normal.copy.number.y normal copy number on the y chromosome
#' @param normal.b.allele.norm normal B allele count
#' @param joined.info list of vectors containing prior information on copy number changes that occurred together.  For example, if it is known that the copy number change is present in sample 'A' and 'B' set joined.info <- list(c("A", "B")). If there's no information, set joined.info <- NA
#' @return A list containing the inferred copy numbers and B-alleles per sample as well as a matrix indicating in which clade of the subclone the copy number change is present. Also returns the squared error of the best solution. 
#' @export
#' @examples 

InferCopyNumber <- function(tree, coverage.ratios, bafs, sample.indicator, sex, locus, parms, normal.copy.number.auto = 2, normal.copy.number.y = 1,
                            normal.b.allele.norm = 1, joined.info = NA){
  
  if(!locus[1] %in% c("X", "Y") | sex == "female"){
    normal.baf <- 1/(normal.copy.number.auto/2)
    normal.copy.number <- normal.copy.number.auto
    normal.copy.number.normal <- 2
    normal.b.allele <- floor(normal.copy.number.auto/2)
  }else{
    normal.baf <- 0
    normal.copy.number <- normal.copy.number.y
    normal.copy.number.normal <- 1
    normal.b.allele <- 0
  }
  
  ## store the inferred B allele per subclone
  b.allele.list <- lapply(1:length(unique(sample.indicator)), function(x){normal.b.allele})
  ## store the inferred copy number per subclone
  copy.number.list <- lapply(1:length(unique(sample.indicator)), function(x){normal.copy.number})
  ## store the combination of subclones harboring the CNVs. A matrix with multiple rows in case there are multiple optimal solutions
  best.sol.list <- lapply(1:length(unique(sample.indicator)), function(x){NA})
  
  ## CNVs can be the same in different clades or separate events 
  indicator <- "separated"
  ## condensed tree (merging subclones per sample)
  pseudo.tree <- (sapply(unique(sample.indicator), function(x){
    rowSums(tree[,sample.indicator==x,drop=F])})!=0)*1
  
  ## list, containing combinations of clones harboring the CNV that is in agreement with the data
  tree.list <- list(0)
  
  ## if all samples have increased coverage ratios and increased or normal or decreased BAFs, and there is a solution in the tree supporting this as a single event, assume a single event
  if(any(rowSums(pseudo.tree[,which(coverage.ratios>1.05),drop=F])==sum(coverage.ratios>1.05) & (all(is.na(bafs)) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs>=normal.baf),drop=F]) == sum(bafs>=normal.baf) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs<=normal.baf),drop=F]) == sum(bafs<=normal.baf) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs==normal.baf),drop=F]) == sum(bafs==normal.baf)))){
    indicator <- "joined"
    
    ## add the combination of clones harboring the CNV
    tree.list[[length(tree.list)+1]] <- tree[which(rowSums(pseudo.tree[,which(coverage.ratios>1.05),drop=F])==sum(coverage.ratios>1.05)),,drop=F]
    
  }
  ## same for losses
  if(any(rowSums(pseudo.tree[,which(coverage.ratios<0.95),drop=F])==sum(coverage.ratios<0.95) & (all(is.na(bafs)) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs>=normal.baf),drop=F]) == sum(bafs>=normal.baf) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs<=normal.baf),drop=F]) == sum(bafs<=normal.baf) ||
                                                                                                 rowSums(pseudo.tree[,which(bafs==normal.baf),drop=F]) == sum(bafs==normal.baf) ))){
    
    indicator <- "joined"
    tree.list[[length(tree.list)+1]] <- tree[which(rowSums(pseudo.tree[,which(coverage.ratios<0.95),drop=F])==sum(coverage.ratios<0.95)),,drop=F]
  }
  
  ## Consider the cases in which we have a normal copy number but a change in BAF (LOH for example)
  if(all(!is.na(bafs))){
    bafs. <- replace(bafs, which(coverage.ratios < 0.95 | coverage.ratios > 1.05), normal.baf)
    
    if(any(bafs.>normal.baf) && 
       any(rowSums(pseudo.tree[,which(bafs.>normal.baf),drop=F])==sum(bafs.>normal.baf))){
      indicator <- "joined"
      
      tree.list[[length(tree.list)+1]] <- tree[which(rowSums(pseudo.tree[,which(bafs.>normal.baf),drop=F]==sum(bafs.>normal.baf)) &
                                                       colSums((coverage.ratios>1.05) == t(pseudo.tree))!=ncol(pseudo.tree) &
                                                       colSums((coverage.ratios<0.95) == t(pseudo.tree))!=ncol(pseudo.tree)),,drop=F]
      
    }
    
    if(any(bafs. < normal.baf) && 
       rowSums(pseudo.tree[,which(bafs.<normal.baf),drop=F]==sum(bafs.<normal.baf))){
      indicator <- "joined"
      
      tree.list[[length(tree.list)+1]] <- tree[which(rowSums(pseudo.tree[,which(bafs.<normal.baf),drop=F]==sum(bafs.<normal.baf)) &
                                                       colSums((coverage.ratios>1.05) == t(pseudo.tree))!=ncol(pseudo.tree) &
                                                       colSums((coverage.ratios<0.95) == t(pseudo.tree))!=ncol(pseudo.tree)),,drop=F]
      
    }
  }else{
    bafs. <- rep(NA, length(baf))
  }
  
  ## if we have both, higher and lower ratios, a joined solution does not work!
  if(any(coverage.ratios>1.05) & any(coverage.ratios < 0.95)){
    indicator <- "separated"
  }
  
  tree.list[[1]] <- NULL
  
  ## include single clades if they agree with the prior information. E.g., if we know that the same CNV is present in sample A and B 
  ## a solution assigning the CNV to a single clade only is not possible. 
  tree <- tree[rowSums(pseudo.tree)==1,,drop=F]
  if(!is.na(joined.info[[1]][1])){
    to.keep <- rep(F, nrow(tree))
    for(i in 1:length(joined.info)){
      to.keep. <- rep(T, nrow(tree))
      tmp <- tree[,sample.indicator %in% joined.info[[i]],drop=F]
      ind. <- sample.indicator[sample.indicator %in% joined.info[[i]]]
      for(j in joined.info[[i]]){
        to.keep.[rowSums(tmp[,which(ind.==j),drop=F])==0] <- F
      }
      to.keep <- to.keep + to.keep.
    }
    to.keep <- to.keep>0
    
    tree <- tree[to.keep,,drop=F]
  }
  
  ## I. test the combinations stored in "tree" (combinations of assigning the CNV to a single sample only); separated solutions
  if(indicator%in%c("joined", "separated") & nrow(tree)!=0){
    
    ## iterate through the samples
    for(i in 1:length(coverage.ratios)){
      
      if(locus[1]=="Y" | (locus[1]=="X" & sex == "male")){
        copynumber <- 1
        copynumber. <- 0
        b.allele <- 0
        b.allele. <- 0
        b.allele.norm <- 0
        copynumber.norm <- normal.copy.number.y
      }else{
        copynumber <- 1
        copynumber. <- 0
        b.allele <- 1
        b.allele. <- 1
        b.allele.norm <- floor(normal.copy.number.auto/2)
        normal.b.allele.norm <- 1
        copynumber.norm <- normal.copy.number.auto
      }
      
      ## initiate with copy number 1 in case of a loss and with the normal copy number else
      if(coverage.ratios[i] < 0.95){
        copynumber <- 1
      }else{
        copynumber <- copynumber.norm
      }
      
      ## copy numbe rindicator matrix
      p <- tree[,sample.indicator==LETTERS[i],drop=F]
      p <- p[rowSums(p)!=0,,drop=F]
      if(nrow(p)==0){
        squared.error <- Inf
        break
      }
      parms. <- parms[sample.indicator==LETTERS[i]]
      
      squared.error.cn <- Inf
      ## report the best solution of which subclones carry the CNV; initiate with the first possible combination
      best.sol.all <- p[1,]
      while(copynumber <= 100 ){
        squared.error <- Inf
        best.sol <- 1
        ## iterate through all possible combinations of subclones carrying the CNV
        for(j in 1:nrow(p)){
          ## squared error of the copy number change
          squared.error. <- ((sum(p[j,]*copynumber*parms.) + sum((1-p[j,])*parms.*copynumber.norm) + (1-sum(parms.))*normal.copy.number.normal)/
                               ((sum(parms.)*copynumber.norm) + sum((1-sum(parms.))*normal.copy.number.normal)) - coverage.ratios[i])^2
          
          
          ## test all possible b.allele numbers
          if(all(!is.na(bafs))){
            b.alleles <- seq(0, copynumber)
            baf <- rep(0, length(b.alleles))
            if(sum(c(p[j,],1)*c(parms., 1-sum(parms.)))!=0){
              baf <- sapply(b.alleles, function(x){
                (sum(p[j,]*x*parms.) + sum((1-p[j,])*b.allele.norm*parms.) +normal.b.allele.norm*(1-sum(parms.)))/
                  (sum(p[j,]*copynumber*parms.) + sum((1-p[j,])*copynumber.norm*parms.) + normal.copy.number.normal*(1-sum(parms.)))
              })
            }
            
            
            b.allele.sq.err <- (baf - bafs[i])^2 
            squared.error. <- squared.error. + min(b.allele.sq.err)
          }
          
          if(squared.error. < squared.error){
            best.sol <- j
            squared.error <- squared.error.
            if(all(!is.na(bafs))){
              b.allele.. <- b.alleles[which.min(b.allele.sq.err)]
            }else{
              b.allele.. <- b.allele.norm
              if(b.allele.. > copynumber){
                b.allele.. <- 0
              }
            }
            
          }
          
        }
        
        ## level of resolution: 10 % copy number change of +-1, corresponds to a squared error difference of 0.01 (default is no change)
        
        if(squared.error.cn > squared.error & abs(squared.error.cn - squared.error) >= 0.01){
          squared.error.cn <- squared.error
          best.sol.all <- p[best.sol,]
          copynumber. <- copynumber
          b.allele. <- b.allele..
        }else{
          break
        }
        
        copynumber <- copynumber + 1
        
      }
      
      copy.number.list[[i]] <- copynumber.
      
      if(all(is.na(bafs))){
        if(copynumber.%%copynumber.norm==0){
          b.allele.list[[i]] <- b.allele.*copynumber./copynumber.norm
          
          if(copynumber.norm %in% c(1,2)){
            b.allele.list[[i]] <- 1
          }
        }else{
          b.allele.list[[i]] <- b.allele.
        }
      }else{
        b.allele.list[[i]] <- b.allele.
      }
      
      
      
      if(copynumber. == copynumber.norm & (is.na(b.allele.) || b.allele. == b.allele.norm)){
        best.sol.all <- NA
        best.sol.list[[i]] <- NA
      }else{
        ## store the combinations of subclones in the tree that are in accordance with the best solution for each sample
        best.sol.all <- tree[which(colSums(t(tree[,sample.indicator==LETTERS[i],drop=F])==best.sol.all)==length(best.sol.all) & 
                                     rowSums(tree[,sample.indicator!=LETTERS[i],drop=F])==0),]
        
        best.sol.list[[i]] <- unique(matrix(best.sol.all, ncol=length(sample.indicator)))
      }
      
    }
  }else{
    squared.error.cn = Inf
  }
  
  solution0 <- list(copy.number.list = copy.number.list, b.allele.list=b.allele.list, best.sol.list = best.sol.list, squared.error.cn = squared.error.cn)
  
  ## II. test the combinations stored in "tree" (combinations of assigning the CNV to multiple samples); joined solutions
  if(indicator %in% c("joined", "separated") & length(tree.list)>0){
    
    ## iterate through tree.list
    for(k in 1:length(tree.list)){
      tree <- tree.list[[k]]
      
      ## select for solutions conformint to the prior information
      if(!is.na(joined.info[[1]][1])){
        to.keep <- rep(F, nrow(tree))
        for(i in 1:length(joined.info)){
          to.keep. <- rep(T, nrow(tree))
          tmp <- tree[,sample.indicator %in% joined.info[[i]],drop=F]
          ind. <- sample.indicator[sample.indicator %in% joined.info[[i]]]
          for(j in joined.info[[i]]){
            to.keep.[rowSums(tmp[,which(ind.==j),drop=F])==0] <- F
          }
          to.keep <- to.keep + to.keep.
        }
        to.keep <- to.keep>0
        
        tree <- tree[to.keep,,drop=F]
      }
      if(nrow(tree)==0){
        squared.error.cn = Inf
        break
      }
      
      ## collapse tree per sample
      pseudo.tree. <- matrix((sapply(unique(sample.indicator), function(x){
        rowSums(tree[,sample.indicator==x,drop=F])})!=0)*1, ncol=ncol(pseudo.tree))
      
      
      if(locus[1]=="Y" | (locus[1]=="X" & sex == "male")){
        copynumber <- 1
        copynumber. <- 0
        b.allele <- 0
        b.allele. <- 0
        b.allele.norm <- 0
        normal.b.allele.norm <- normal.b.allele.norm
        copynumber.norm <- normal.copy.number.y
      }else{
        copynumber <- 1
        copynumber. <- 0
        b.allele <- 1
        b.allele. <- 1
        b.allele.norm <- floor(normal.copy.number.auto/2)
        normal.b.allele.norm <- normal.b.allele.norm
        copynumber.norm <- normal.copy.number.auto
      }
      
      
      ## initiate with CN=1 in case of a loss and with the normal copy number else
      if(coverage.ratios[which(colSums(pseudo.tree.)!=0)[1]]>1.05){
        copynumber <- copynumber.norm
      }else{
        copynumber <- 1
      }
      
      ## find the best solutions by comparing expected and measured coverage ratios per sample combination, as before
      squared.error.cn <- Inf
      best.sol.all <- tree[1,]
      while(copynumber <= 100 ){
        squared.error <- Inf
        best.sol <- 1
        
        squared.error. <- 0
        for(j in 1:nrow(tree)){
          
          b.allele.sq.err. <- 0
          for(i in 1:length(coverage.ratios)){
            p <- tree[j,sample.indicator==LETTERS[i]]
            
            parms. <- parms[sample.indicator==LETTERS[i]]
            
            ## squared error of the copy number change
            squared.error. <- squared.error. + 
              ((sum(p*copynumber*parms.) + sum((1-p)*c(parms.)*copynumber.norm) + (1-sum(parms.))*normal.copy.number.normal)/
                 (copynumber.norm*sum(parms.) + normal.copy.number.normal*(1-sum(parms.))) - coverage.ratios[i])^2 
            
            ## test all possible b.allele numbers
            if(all(!is.na(bafs))){
              b.alleles <- seq(0, copynumber)
              baf <- rep(0, length(b.alleles))
              if(sum(c(p,1)*c(parms., 1-sum(parms.)))!=0){
                baf <- sapply(b.alleles, function(x){
                  (sum(p*x*parms.) + sum((1-p)*b.allele.norm*parms.) + normal.b.allele.norm*(1-sum(parms.)))/
                    (sum(p*copynumber*parms.) + sum(1-p)*copynumber.norm*parms. + normal.copy.number.normal*(1-sum(parms.)))
                })
              }
              
              
              b.allele.sq.err <- b.allele.sq.err. + sum((baf - bafs[i])^2)
            }
            
          }
          if(all(!is.na(bafs))){
            squared.error. <- squared.error. + min(b.allele.sq.err)
          }
          
          if(squared.error. < squared.error){
            best.sol <- j
            squared.error <- squared.error.
            if(all(!is.na(bafs))){
              b.allele.. <- b.alleles[which.min(b.allele.sq.err)]
            }else{
              b.allele.. <- b.allele.norm
              if(b.allele.. >= copynumber){
                b.allele.. <- 0
              }
            }
          }
          
        }
        
        ## level of resolution: 10 % copy number change of +-1, corresponds to a squared error difference of 0.01 (default is no change)
        
        if(squared.error.cn > squared.error & abs(squared.error.cn - squared.error) >= 0.01){
          squared.error.cn <- squared.error
          best.sol.all <- tree[best.sol,]
          copynumber. <- copynumber
          b.allele. <- b.allele..
        }else{
          break
        }
        
        copynumber <- copynumber + 1
        
      }
      
      pseudo.tree. <- matrix((sapply(unique(sample.indicator), function(x){
        rowSums(tree[,sample.indicator==x,drop=F])})!=0)*1, ncol=ncol(pseudo.tree))
      
      for(counter in which(colSums(pseudo.tree.)!=0)){
        copy.number.list[[counter]] <- copynumber.
        b.allele.list[[counter]] <- b.allele.
        
        if(all(is.na(bafs))){
          if(copynumber.%%copynumber.norm==0){
            b.allele.list[[counter]] <- b.allele.*copynumber./copynumber.norm
            if(copynumber.norm %in% c(1,2)){
              b.allele.list[[counter]] <- 1
            }
          }else{
            b.allele.list[[counter]] <- b.allele.
          }
        }else{
          b.allele.list[[counter]] <- b.allele.
        }
        
        if(copynumber. == copynumber.norm & (is.na(b.allele.) || b.allele. == b.allele.norm)){
          
          best.sol.all <- NA
        }
        
        best.sol.list[[counter]] <- unique(matrix(best.sol.all, ncol=length(sample.indicator)))
      }
      
      
    }
    
  }else{squared.error.cn = Inf}
  
  
  solution1 <- list(copy.number.list = copy.number.list, b.allele.list=b.allele.list, best.sol.list = best.sol.list, squared.error.cn = squared.error.cn)
  
  
  ## decide whether solution0 or solution1 is better 
  if(2*solution0$squared.error.cn < solution1$squared.error.cn ){
    return(solution0)
  }
  
  
  if(!is.infinite(solution1$squared.error.cn)){
    return(solution1)
    
  }else{
    return("no fit")
  }
  
  
}








