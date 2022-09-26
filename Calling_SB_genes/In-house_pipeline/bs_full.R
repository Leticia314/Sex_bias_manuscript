# Script for identifying sex-biased genes (in-house pipeline)
# Author: Svetlana Ovchinnikova (s.ovchinnikova@zmbh.uni-heidelberg.de)

library(locfit)
library(future.apply)
plan(multisession, workers = 16) #number of R sessions running in parallel to speed up the calculations

tissues <- c("Brain", "Cerebellum", "Heart", "Kidney", "Liver")
species <- c("Human", "Mouse", "Rat", "Rabbit", "Chicken", "Opossum")

for(sp in species) {
  print(sp)
  
  #loading the data 
  data <- logTransform(loadData(sp))
  #adding letters in front of timepoints' name. They are later used as names for a vector and thus shouldn't be just numbers
  if(sp == "Rabbit")
    data$group <- paste0("r", data$group)
  if(sp == "Opossum")
    data$group <- paste0("o", data$group)
  
  #timewarping (it's already precalculated and stored in files)
  #adds new column to the sample table
  data <- newTimePoints(data, sp)
  
  for(t in tissues) {
    print(t)
    #sample table for a specific tissue
    st <- subset(as.data.frame(colData(data)), tissue == t)
    #get only the time points where we have both male and female samples
    st %>%
      group_by(timePoint, group) %>%
      summarise(femCount = sum(sex == "Female"), maleCount = sum(sex == "Male")) %>%
      filter(femCount != 0 & maleCount != 0) %>%
      select(timePoint, group) -> tps
    
    st <- subset(st, timePoint %in% tps$timePoint)
    
    #get expression values and remove genes that are expressed in less than 3 samples
    exprs <- assay(data, "logcounts")[, rownames(st)]
    exprs <- exprs[rowSums(exprs != 0) > 2, ]
    
    #define threshold
    thr <- 0.1
    
    #simple function that returns an element with the largest absolute value
    signedmax <- function(x) x[ which.max(abs(x)) ] 
    
    l <- 0 #number of resamplings
    # this variable stores the results of bootstraping
    # 1 - number of times we got above thr
    # 2 - number of times we got below -thr
    # 3 - total number of resamplings for this gene
    res <- matrix(0, nrow = nrow(exprs), ncol = 3) 
    rownames(res) <- rownames(exprs)
    stop <- rep(F, nrow(exprs)) # for each gene: whether we should stop resampling this gene
    names(stop) <- rownames(exprs)
    i <- 0
    tm <- Sys.time()
    tm_all <- Sys.time()
    estPerm <- 500 # estimated number of resamplings needed for reasonable multiple testing correction
    
    # we will stop either when the estimated number of resamplings is reached
    # or when for each gene we have an estimate for p-value (never happens actually)
    while(sum(!stop) > 0 & l < estPerm) { 
      i <- i + 1
      
      print(paste0("Round ", i))
      print(paste0("Number of permutations: ", l))
      l0 <- l
      # the smaller number of genes have significant difference
      # between sexes, the more number of resamplings is required to reach
      # FDR value of roughly 0.1 (actually a smaller value is used)
      # estPerm is always calculated based on the current number of promising 
      # genes, however this number will decrease (more resamplings will show that some of
      # this promising genes are not in fact significant). Therefore, on each step we do more
      # resamplings, than estimated
      l <- round(1.1 * estPerm)
      
      # We will do resaplings in batches. This is the size of the batch. It is somewhere
      # between 25 and 1000  and we aim to do about 5% of planned resamplings in one batch
      resampl <- min(1000, max(round((l - l0) * 0.05), 25))
      
      print(paste0("Genes left: ", sum(!stop)))
      print(paste0("Estimated number of permutations: ", estPerm))
      print(Sys.time() - tm)
      print("-------------------------------------------")
      tm <- Sys.time()
      
      tps <- sort(unique(st$newTimePoint))
      # band width for locfit regression
      # the bigger is this number, the smoother will be the curve
      # it is defined so that at any point we use at least 3 time points for curve fitting
      # it looks complicated because of the time warping - without it it would be just about 2.5-3
      h <- max(diff(range(st$newTimePoint))/4, (tps[3] - tps[1]) * 1.25, (tps[length(tps)] - tps[length(tps) - 2] ) * 1.25)
      
      # bootstraping
      # each gene for which we still don't have an estimate for p-value will be processed in
      # parallel
      newRes <- t(future_apply(cbind(res, exprs)[!stop, ], 1, function(expr) {
        perms <- expr[3]
        hitsUp <- expr[1]
        hitsDown <- expr[2]
        #expression of the gene
        expr <- expr[-(1:3)]
        
        dataTable <- st
        dataTable$expr <- expr
        
        dataTableM <- subset(dataTable, sex == "Male")
        dataTableF <- subset(dataTable, sex == "Female")
        
        #startTime <- Sys.time()
        # "design matrix" (sort of) for both curves (male and female separately) 
        xM <- lp(dataTableM$newTimePoint, deg = 1, h = h)
        xF <- lp(dataTableF$newTimePoint, deg = 1, h = h)
        
        # timepoints where we are going to get the fitted values for each curve
        # we are using unexported and direct calls instead of 'locfit' function to improve performance time
        grid <- lfgrid(mg = 30, ll = min(st$newTimePoint), ur = max(st$newTimePoint))
        # we will keep resampling this gene until we reach a required number of resamplings
        # until we step across the threshold at least five times or
        # unitl we reach the planned number of permutations
        while((hitsUp < 5 | hitsDown < 5) & perms < l) { 
          
          res <- t( replicate( resampl, {
            # fitting the curve with weighted regression
            # weights are randomly drawn from a gamma distribution
            fitM <- locfit.raw(xM, dataTableM$expr, weights = rgamma( nrow(dataTableM), shape=1, scale=1 ), ev = grid )
            fitF <- locfit.raw(xF, dataTableF$expr, weights = rgamma( nrow(dataTableF), shape=1, scale=1 ), ev = grid )
            
            # evaluate models at predefined points (that is an ugly but the most fast and direct way I've found to get these values)
            male <- (locfit:::preplot.locfit.raw(fitM, NULL, "fitp", "coef", "none"))$y
            female <- (locfit:::preplot.locfit.raw(fitF, NULL, "fitp", "coef", "none"))$y
            
            # average distance between curves
            dist <- mean(male - female)
            #keep only the values that have the same sign as dist
            vals <- male - female
            vals <- abs(vals[vals * dist >= 0])
            
            #maximal and average distance are returned
            c( dist = dist,
               maxDist = max(vals))
          }))
          
          #calculating the score
          score <- sqrt(abs(res[, 1] * res[, 2])) * sign(res[, 1])
          
          # check, how many times we step across the threshold and encrease the number of performed resamplings
          hitsUp <- hitsUp + sum(score < thr)
          hitsDown <- hitsDown + sum(score > -thr)
          
          perms <- perms + resampl
        }
        
        #print(Sys.time() - start)
        #c(as.numeric(startTime), as.numeric(Sys.time()), as.numeric(Sys.time() - startTime), Sys.getpid())
        #c((hitsUp + 1)/(perms + 1), (hitsDown + 1)/(perms + 1), perms)
        c(hitsUp, hitsDown, perms)
      }))
      # update results for the genes we've tested at this step
      res[!stop, ] <- newRes
      # check, which genes still need more resamplings
      stop <- pmin(res[, 1], res[, 2]) >= 5
      # based on this number, update estimated number of required permutations
      estPerm <- round(500 * nrow(exprs)/max(sum(!stop), 1))
    }
    print(Sys.time() - tm_all)
    write_rds(res, paste0(resPath, sp, "_", t, "_bootstrap_0.1_notw.rds"), compress = "gz")
  }
}
