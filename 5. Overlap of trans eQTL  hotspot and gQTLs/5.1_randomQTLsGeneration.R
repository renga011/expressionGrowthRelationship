#!/usr/bin/env Rscript --vanilla

## script to generate 1000 sets of random-EQTLs and random-growthQTLs for whole genome QTL raindown 

#commands required for MSI 
setwd("~/ExpressionPhenotypeProject/")

library(dplyr)
library(stringr)
library(parallel)
library(GenomicRanges)

#GLOBAL VARIABLES ----

results_dir = "results/updates_qValue_justPhenoGeneticCorFilters_112421/"
RObj_dir = "RObjects_092522/"

nRandomSets = 1000
nIterations = 1000

#load datasets ----
load(paste0(results_dir, RObj_dir, "dataWithRPAndSEForExpressionAndGrowth.rda"))
load(paste0(results_dir, RObj_dir, "chromosomeSizes.rda"))

#variables from command line -------
batchIndex = commandArgs(trailingOnly = TRUE)
batchIndex = as.numeric(batchIndex)

#we run 40 random sets in one job
randomSet_start = (batchIndex-1) * 40 + 1
randomSet_end = batchIndex * 40

##### randomGrowthQTLs generation ---

#functions to generate random QTLs---

chooseChrStartEnd = function(chromosomeTable, width, forbiddenStartSpace){
  chromosomeTable = dplyr::filter(chromosomeTable, chromosome %in% 1:16) #only 16 chromosomes considered
  chromosomeTable$size_new = chromosomeTable$size - sapply(forbiddenStartSpace, function(x){return(length(x))}) #reduced chromosome length after deleting the forbidden spaces.
  
  #choose chromosome
  possibleChromosomeChoices = chromosomeTable$chromosome[which(chromosomeTable$size_new > 0)]
  
  if(!(rlang::is_empty(possibleChromosomeChoices))){ #conditional to check if any chromosomes are left to choose
  chromosomeChoice = ifelse(length(possibleChromosomeChoices) > 1, sample(possibleChromosomeChoices, 1), possibleChromosomeChoices) #have to include this conditional because sample does some weird thing when only one number is supplied for x
  
  #choose start point
  forbiddenStartSpace_i = forbiddenStartSpace[[chromosomeChoice]]
  
  availableStartSpace = setdiff(1:chromosomeTable$size[chromosomeTable$chromosome == chromosomeChoice], forbiddenStartSpace_i)
  
  startPointChoice = ifelse(length(availableStartSpace) > 1, sample(availableStartSpace, 1), availableStartSpace) 
  
  #choose end point
  endPointChoice = startPointChoice + width - 1
  
  genomeSpaceIndicator = 1
  } else {
    #when no more genome space is available to make QTL assignments, you return NAs and set the genome space indicator as FALSE => no more space available.
    chromosomeChoice = NA
    startPointChoice = NA
    endPointChoice = NA
    genomeSpaceIndicator = 0
  }
  
  return(c("chr_i" = chromosomeChoice,
           "start_i" = startPointChoice,
           "end_i" = endPointChoice,
           "genomeSpaceIndicator_i" = genomeSpaceIndicator))
}

checkOverlap = function(query_i, subject_i){
  queryRanges = GRanges(seqnames = query_i$chromosome,
                        IRanges(start = query_i$CILeft,
                                end = query_i$CIRight))
  
  subjectRanges = GRanges(seqnames = subject_i$chromosome,
                          IRanges(start = subject_i$CILeft,
                                  end = subject_i$CIRight))
  
  overlap = GenomicRanges::findOverlaps(queryRanges, subjectRanges)
  #x = ifelse(length(overlap) !=0, "TRUE", "FALSE")
  return(length(overlap))
}

randomQTLGenerator = function(variable_i, QTLTable){ #variable is either gene or trait here
  #get the QTLs corresponding to the variable_i from the QTL table
  QTLs_i = dplyr::filter(QTLTable, variable == variable_i)
  QTLs_i = QTLs_i[order(QTLs_i$width, decreasing = TRUE),]
  
  iteration = 1
  repeatTicker = 1 #initialize repeatTicker as 1 to enter the loop
  
  while(iteration < nIterations & repeatTicker == 1){
    rQTL = as.data.frame(matrix(0, nrow = nrow(QTLs_i), ncol = ncol(QTLs_i)))
    colnames(rQTL) = colnames(QTLs_i)
    rQTL$variable = QTLs_i$variable
    
    #do the random QTL assignment for the first QTL in the list 
    i = 1
    
    chr_i = QTLs_i$chromosome[i]
    start_i = QTLs_i$CILeft[i]
    end_i = QTLs_i$CIRight[i]
    width_i = QTLs_i$width[i]
    
  #forbidden start spaces for this QTL across chromosomes so that the QTL doesn't slide off the edge of the chromosome
    forbiddenStartSpace_i = lapply(1:16, function(chr_no){
        chr_length = chromosomeSizes$size[which(chromosomeSizes$chromosome == chr_no)]
          maxStartPoint = chr_length -width_i + 1
          maxStartPoint = ifelse(maxStartPoint < 1, 0, maxStartPoint) #if the chr_length is less than the width of the eQTL interval, entire chromosome is forbidden
          forbiddenSpace_i = (maxStartPoint+1):chr_length
          return(forbiddenSpace_i)
        })
    names(forbiddenStartSpace_i) = 1:16

  #generate the random QTL
  rQTL_i = chooseChrStartEnd(chromosomeTable = chromosomeSizes,
                            width = width_i,
                            forbiddenStartSpace = forbiddenStartSpace_i)

  rQTL[i,-1] = data.frame(chromosome = rQTL_i["chr_i"],
                          CILeft = rQTL_i["start_i"],
                          CIRight = rQTL_i["end_i"],
                          width = rQTL_i["end_i"] - rQTL_i["start_i"] + 1)

  iteration = 1
  genomeSpaceIndicator = rQTL_i["genomeSpaceIndicator_i"] #a variable to keep track whether sufficient genome space is left to generate the next non-overlapping random-QTL
    
  while(i < nrow(QTLs_i)){
      i = i+1 #increment the row number
      chr_i = QTLs_i$chromosome[i]
      start_i = QTLs_i$CILeft[i]
      end_i = QTLs_i$CIRight[i]
      width_i = QTLs_i$width[i]
      
      #forbidden start spaces for this QTL across chromosomes so that the QTL doesn't slide off the edge of the chromosome
      forbiddenStartSpace_i = lapply(1:16, function(chr_no){
        chr_length = chromosomeSizes$size[which(chromosomeSizes$chromosome == chr_no)]
        maxStartPoint = chr_length -width_i + 1
        maxStartPoint = ifelse(maxStartPoint < 1, 0, maxStartPoint) #if the chr_length is less than the width of the eQTL interval, entire chromosome is forbidden
        forbiddenSpace_i = (maxStartPoint+1):chr_length
        return(forbiddenSpace_i)
      })
      names(forbiddenStartSpace_i) = 1:16
      
      forbiddenStartSpace_i[[rQTL_i["chr_i"]]] = unique(c(forbiddenStartSpace_i[[rQTL_i["chr_i"]]], rQTL_i["start_i"]: rQTL_i["end_i"]))
      
      rQTL_i = chooseChrStartEnd(chromosomeTable = chromosomeSizes,
                                 width = width_i,
                                 forbiddenStartSpace = forbiddenStartSpace_i)
     
      
      #check overlap once for confirmation that things don't overlap with whats already been selected.
      if(!(is.na(rQTL_i["chr_i"]))){
        
          rQTL[i,-1] = c(chromosome = rQTL_i["chr_i"],
                                  CILeft = rQTL_i["start_i"],
                                  CIRight = rQTL_i["end_i"],
                                  width = rQTL_i["end_i"] - rQTL_i["start_i"] + 1)
      }
      if(is.na(rQTL_i["chr_i"])){
        break
      }
      }
      
  #if there is no more genome space for next QTL assignment, break out of the loop and go again
  genomeSpaceIndicator = rQTL_i["genomeSpaceIndicator_i"]
  if(genomeSpaceIndicator == 0){
      repeatTicker = 1
      break
      } else{
        repeatTicker = 0 #if there is genome space available then you don't need to go again. so this loop can finish and the bigger loop can end
      }
    
    iteration = iteration + 1
  }
  
  #print output
  if(repeatTicker == 0){
    print(paste(variable_i, ": randomQTLs assigned for all QTLs in", 
                iteration-1, "iterations"))
    #print(paste("overlaps within the randomQTL set generated:" , checkOverlap(rQTL, rQTL) - nrow(rQTL)))
    } else if(repeatTicker == 1){
                  print(paste(variable_i, ": randomQTLs couldn't be assigned for all QTLs in", 
                              iteration-1, "iterations")) 
                }
  
  return(rQTL)
}

#random growth-QTLs ------

growthQTLs = dplyr::select(BloomQTL, Trait, Chromosome, LOD_leftCI,
                           LOD_rightCI)
colnames(growthQTLs) = c("variable", "chromosome", "CILeft", "CIRight")

growthQTLs$width = growthQTLs$CIRight - growthQTLs$CILeft + 1

#apply QTL generator to growth-QTL table 

randomGQTLs = parallel::mclapply(randomSet_start:randomSet_end,
                                 FUN = function(randomSet_i){
                                   print(randomSet_i)
                                   rQTLs = mclapply(unique(growthQTLs$variable),
                                                    FUN = randomQTLGenerator,
                                                    QTLTable = growthQTLs,
                                                    mc.cores = parallelly::availableCores())
                                   rQTLs = do.call(rbind, rQTLs)
                                   return(rQTLs)
                                 },
                                 mc.cores = parallelly::availableCores())

#save individual batches of randomGQTLs

save(randomGQTLs,
     file = paste0(results_dir, RObj_dir,
                   "randomGQTLs_", batchIndex,".rda"))

#combine individual batches of random QTLs into one object

  randomGQTLs_all = list()
  
  for(batchIndex in 1:25){
    load(paste0(results_dir, RObj_dir,
                "randomGQTLs_", batchIndex,".rda"))
    randomGQTLs_all = c(randomGQTLs_all, randomGQTLs)
  }
  names(randomGQTLs_all) = 1:1000


#save ----

save(randomGQTLs_all,
     file = paste0(results_dir, RObj_dir, "randomGQTLs.rda"))
