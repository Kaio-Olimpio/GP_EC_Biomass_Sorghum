##' Function assignation
##'
##' @title Preparing the data for cross-validation
##' 
##' @description
##' This function divides the data into k folds and repeat this process according
##' to the number of folds and repetitions provided by the user. It also stores 
##' all the necessary information to perform the cross-validation using [BGLR::BGLR()]. 
##' 
##' @param data A data frame containing the observed data
##' @param y The name of the column that contains the trait values
##' @param gen The name of the column that contains the genotypes information
##' @param env The name of the column that contains the environments information
##' @param niter An integer indicating the number of iteration
##' @param burnin An integer indicating the number of burn-in iterations
##' @param seed An integer indicating the seed to be used. Mandatory, for the sake
##' of reproducibility
##' @param nfolds An integer indicating the number of folds for the cross-validation
##' @param nrept An integer indicating the number of times the cross-validation 
##' should repeat
##' @param cv A CHARACTER: 'cv1' for CV1, 'cv2' for CV2, 'loo_cv0' for classic CV0 using a
##' leave-one-out scheme, and 'loo_cv00' for CV00 using a leave-one-out scheme
##' 
##' @return The function will return a dataframe named `cvinfo` with the 
##' information necessary to perform the cross-validation using [BGLR::BGLR()], 
##' and a list called `cvdata` containing the data that will be used in the 
##' cross-validations. 
##' 
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##'

assignation = function(data, y, gen, env, seed, nfolds,
                       nrept, cv, niter, burnin){
  
  data = data
  data$gen = data[,gen]
  data$env = data[,env]
  data$trait = data[,y]
  
  if (cv == 'cv1' | cv == 'cv00'){
    set.seed(seed)
    sets = split(
      rep(
        1:nfolds, length(unique(data$gen)) * nrept
      )[order(runif(length(unique(data$gen)) * nrept))],
      f = 1:nrept
    )
    cvdata = lapply(sets, function(x){
      cvdata = data[,c('gen','env','trait')]
      cvdata = merge(cvdata, data.frame(
        gen = unique(data$gen),
        set = x
      ), by = gen)
    })
    names(cvdata) = paste0('R',1:nrept)
    
  }else if (cv == 'cv2' | cv == 'cv0'){
    cvdata = list()
    for (j in 1:nrept) {
      set.seed(seed+j)
      cvdata[[j]] = data[,c('gen','env','trait')]
      cvdata[[j]]$set = NA
      for (i in unique(data$gen)) {
        cvdata[[j]][cvdata[[j]]$gen == i,'set'] = sample(
          1:nfolds, 
          size = dim(cvdata[[j]][cvdata[[j]]$gen == i,])[1],
          replace = dim(cvdata[[j]][cvdata[[j]]$gen == i,])[1] > nfolds
        )
      }
    }
    names(cvdata) = paste0('R',1:nrept)
  }else if (cv == 'loo_cv0' | cv == 'loo_cv00'){
    cvdata = data[, c('gen', 'env', 'trait')]
  }
  cvinfo = data.frame(nfolds = nfolds, nrept = nrept, cv = cv, niter = niter, burnin = burnin)
  
  return(list(cvdata = cvdata, cvinfo = cvinfo))
}

##' Function crossvalid
##'
##' @title Cross-validation schemes
##' 
##' @description
##' This functions performs parallelized cross-validations for obtaining the predictive ability
##' of a genomic prediction model. The function can perform CV1, CV2, CV0 and CV00.
##' The models are fitted using [BGLR::BGLR()]. Before using `crossvalid`, users may
##' prepare the data using the `assignation()` function.
##' 
##' @details
##' Cross-validations are useful for assessing how good the model can be in predicting
##' the performance of genotypes in distinct scenarios. Each cross-validation 
##' scheme (CV1, CV2, CV0 and CV00) mimics a different situation. At CV1, the goal
##' is to predict the performance of untested genotypes. At CV2, genotypes are 
##' evaluated at some environments, but not all of them, so the goal is to predict
##' the performance of genotypes in environment where it was not tested. CV0 assess 
##' the model's ability to predict the performance of tested genotypes in untested 
##' environments. Finally, CV00 evaluates how good is the model's predictions for 
##' untested genotypes performance in untested environments.
##' 
##' @param cvdata A list obtained from the `assignation()` function
##' @param ETA A two-level list used to specify the regression function (or linear 
##' predictor). See [BGLR::BGLR()] manual for more details.
##' @param results If NULL, the results will be provided per repetition. 
##' If `results = 'mean'`, the function will return the mean across repetitions.
##' @param ... passed to [BGLR::BGLR()].
##' 
##' @return The function will return a list with:
##' \itemize{
##' \item \code{yhat} : A data frame containing the observed and predicted values
##' for each genotype, in each environment.
##' \item \code{corr} : The correlation between observed and predicted values.
##' \item \code{corr_env} : The correlation between observed and predicted values
##'  per environment
##' \item \code{wapa} : The weighted average predictive ability of the model
##' \item \code{mspe} : The mean square prediction error of the model
##' \item \code{slope} : The value of the slope of a linear regression between the 
##' observed and predicted values
##' } 
##' 
##' @import BGLR parallel doParallel foreach
##' 
##' @author Saulo F. S. Chaves (saulo.chaves at ufv.br)
##'

crossvalid = function(cvdata,  ETA, results = NULL,...){
  
  cvdata = cvdata

  if (cvdata$cvinfo$cv == 'cv1'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo$nrept), 
                             f = 1:ceiling(cvdata$cvinfo$nrept/ parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo$nrept), 
              f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F)))[[i]]
        ], function(z){
        lapply(1:cvdata$cvinfo$nfolds, function(x){
          
          yNA = z$trait
          yNA[z$set == x] = NA
          
          crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter, 
                                burnIn = cvdata$cvinfo$burnin, verbose = F)
          unlink("*.dat")
          
          cbind(
            z[z$set == x,], 
            yhat = crossval$yHat[z$set == x]
          )
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    )
    
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    corr_env = lapply(corr_env, function(x) setNames(x, unique(cvdata$cvdata[[1]]$env)))
    corr_env = lapply(corr_env, function(x) x[order(names(x))])
    Vr = lapply(corr_env, function(x) (1-x^2)/(tapply(cvdata$cvdata[[1]]$gen, cvdata$cvdata[[1]]$env, length) - 2))
    wapa = list()
    for (i in names(corr_env)) {
      wapa[[i]] = sum(corr_env[[i]]/Vr[[i]]) / sum(1/Vr[[i]])
    }
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cvdata$cvinfo$cv == 'cv2'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo$nrept), 
                             f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo$nrept), 
              f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F)))[[i]]
        ], function(z){
        lapply(1:cvdata$cvinfo$nfolds, function(x){
          
          yNA = z$trait
          yNA[z$set == x] = NA
          
          crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter,
                                burnIn = cvdata$cvinfo$burnin, verbose = F, ...)
          unlink("*.dat")
          
          cbind(
            z[z$set == x,], 
            yhat = crossval$yHat[z$set == x]
          )
        })
      })
    }

    parallel::stopCluster(cl)
    
    reps = unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    )

    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    corr_env = lapply(corr_env, function(x) setNames(x, unique(cvdata$cvdata[[1]]$env)))
    corr_env = lapply(corr_env, function(x) x[order(names(x))])
    Vr = lapply(corr_env, function(x) (1-x^2)/(tapply(cvdata$cvdata[[1]]$gen, cvdata$cvdata[[1]]$env, length) - 2))
    wapa = list()
    for (i in names(corr_env)) {
      wapa[[i]] = sum(corr_env[[i]]/Vr[[i]]) / sum(1/Vr[[i]])
    }
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if (cvdata$cvinfo$cv == 'cv0'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo$nrept),
                             f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, X = cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo$nrept),
              f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F)))[[i]]
        ], fun = function(z){
        lapply(unique(cvdata$cvdata[[1]]$env), function(w){
          lapply(1:cvdata$cvinfo$nfolds, function(x){
            
            yNA = z$trait
            yNA[z$env == w] = NA
            
            crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter,
                                  burnIn = cvdata$cvinfo$burnin, verbose = T, ...)
            unlink("*.dat")
            
            cbind(
              z[z$set == x & z$env == w,], 
              yhat = crossval$yHat[z$set == x & z$env == w]
            )
          })
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = lapply(unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    ), function(y){do.call(rbind,y)})
    
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    corr_env = lapply(corr_env, function(x) setNames(x, unique(cvdata$cvdata[[1]]$env)))
    corr_env = lapply(corr_env, function(x) x[order(names(x))])
    Vr = lapply(corr_env, function(x) (1-x^2)/(tapply(cvdata$cvdata[[1]]$gen, cvdata$cvdata[[1]]$env, length) - 2))
    wapa = list()
    for (i in names(corr_env)) {
      wapa[[i]] = sum(corr_env[[i]]/Vr[[i]]) / sum(1/Vr[[i]])
    }
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    
    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
    
  }else if (cvdata$cvinfo$cv == 'cv00'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F))
    parallel::clusterEvalQ(cl, library(BGLR))
    parallel::clusterExport(cl, varlist = c('cvdata', 'ETA'))
    
    a = list()
    for (i in 1:length(split(paste0('R',1:cvdata$cvinfo$nrept), 
                             f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F))))) {
      a[[i]] = parallel::parLapply(cl = cl, X = cvdata$cvdata[
        split(paste0('R',1:cvdata$cvinfo$nrept),
              f = 1:ceiling(cvdata$cvinfo$nrept/parallel::detectCores(logical = F)))[[i]]
        ], fun = function(z){
        lapply(unique(cvdata$cvdata[[1]]$env), function(w){
          lapply(1:cvdata$cvinfo$nfolds, function(x){
            
            yNA = z$trait
            yNA[z$env == w & z$set == x] = NA
            
            crossval = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter,
                                  burnIn = cvdata$cvinfo$burnin, verbose = F, ...)
            unlink("*.dat")
            
            cbind(
              z[z$set == x & z$env == w,], 
              yhat = crossval$yHat[z$set == x & z$env == w]
            )
          })
        })
      })
    }
    
    parallel::stopCluster(cl)
    
    reps = lapply(unlist(
      lapply(
        lapply(a, function(x) do.call(rbind, x)),
        function(w){
          apply(w, 1, function(z){do.call(rbind, z)})
        }
      ),recursive = F
    ), function(y){do.call(rbind,y)})
 
    corr = lapply(reps, function(x) cor(x[, 'trait'], x[,'yhat']))
    corr_env = lapply(reps, function(x){
      as.vector(by(x[,c('trait', 'yhat')], x[,'env'], 
                   function(z) cor(z[,'trait'], z[,'yhat'])))
    })
    corr_env = lapply(corr_env, function(x) setNames(x, unique(cvdata$cvdata[[1]]$env)))
    corr_env = lapply(corr_env, function(x) x[order(names(x))])
    Vr = lapply(corr_env, function(x) (1-x^2)/(tapply(cvdata$cvdata[[1]]$gen, cvdata$cvdata[[1]]$env, length) - 2))
    wapa = list()
    for (i in names(corr_env)) {
      wapa[[i]] = sum(corr_env[[i]]/Vr[[i]]) / sum(1/Vr[[i]])
    }
    mspe = lapply(reps, function(x) mean((x[,'yhat'] - x[, 'trait'])^2))
    slope = lapply(reps, function(x) unname(lm(x[,'trait'] ~ x[,'yhat'])$coefficients[2]))
    

    if (is.null(results)){
      res = list(yhat = reps, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = reps,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }

  }else if(cvdata$cvinfo$cv == 'loo_cv0'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F), type = 'PSOCK')
    doParallel::registerDoParallel()
    yhat = foreach::foreach(i = unique(cvdata$cvdata$env)) %dopar% {
      yNA = cvdata$cvdata$trait
      yNA[cvdata$cvdata$env == i] = NA
      
      cvloo = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter,
                         burnIn = cvdata$cvinfo$burnin, verbose = F, ...)
      unlink("*.dat")
      
      cbind(
        cvdata$cvdata[cvdata$cvdata$env == i,],
        yhat = cvloo$yHat[cvdata$cvdata$env == i] 
      )
    }
    parallel::stopCluster(cl)
    
    yhat = do.call(rbind, yhat)
    corr = cor(yhat$trait, yhat$yhat)
    corr_env = as.vector(by(yhat[,c('trait', 'yhat')], yhat[,'env'],
                            function(z) cor(z[,'trait'], z[,'yhat'])))
    names(corr_env) = unique(cvdata$cvdata$env)
    corr_env = corr_env[order(names(corr_env))]
    Vr = (1-corr_env^2)/(tapply(yhat$gen, yhat$env, length) - 2)
    wapa = sum(corr_env/Vr) / sum(1/Vr)
    mspe = mean((yhat[,'yhat'] - yhat[, 'trait'])^2)
    slope = unname(lm(yhat[,'trait'] ~ yhat[,'yhat'])$coefficients[2])
    if (is.null(results)){
      res = list(yhat = yhat, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
      return(res)
    }else{
      res = list(
        yhat = yhat,
        corr = mean(do.call(rbind, corr)),
        corr_env = apply(do.call(cbind,lapply(corr_env,as.matrix)), 1, mean),
        wapa = mean(do.call(rbind, wapa)),
        mspe = mean(do.call(rbind, mspe)),
        slope = mean(do.call(rbind, slope))
      )
      return(res)
    }
  }else if(cvdata$cvinfo$cv == 'loo_cv00'){
    cl = parallel::makeCluster(parallel::detectCores(logical = F), type = 'PSOCK')
    doParallel::registerDoParallel()
    yhat = foreach::foreach(i = unique(cvdata$cvdata$env)) %dopar% {
      yNA = cvdata$cvdata$trait
      yNA[cvdata$cvdata$env == i | cvdata$cvdata$gen %in% cvdata$cvdata[cvdata$cvdata$env == i, 'gen']] = NA
      
      cvloo = BGLR::BGLR(y = yNA, ETA = ETA, nIter = cvdata$cvinfo$niter,
                         burnIn = cvdata$cvinfo$burnin, verbose = F, ...)
      unlink("*.dat")
      
      cbind(
        cvdata$cvdata[cvdata$cvdata$env == i,],
        yhat = cvloo$yHat[cvdata$cvdata$env == i] 
      )
    }
    parallel::stopCluster(cl)
    
    yhat = do.call(rbind, yhat)
    corr = cor(yhat$trait, yhat$yhat)
    corr_env = as.vector(by(yhat[,c('trait', 'yhat')], yhat[,'env'],
                            function(z) cor(z[,'trait'], z[,'yhat'])))
    names(corr_env) = unique(cvdata$cvdata$env)
    corr_env = corr_env[order(names(corr_env))]
    Vr = (1-corr_env^2)/(tapply(yhat$gen, yhat$env, length) - 2)
    wapa = sum(corr_env/Vr) / sum(1/Vr)
    mspe = mean((yhat[,'yhat'] - yhat[, 'trait'])^2)
    slope = unname(lm(yhat[,'trait'] ~ yhat[,'yhat'])$coefficients[2])

    res = list(yhat = yhat, corr = corr, corr_env = corr_env, wapa = wapa, 
                 mspe = mspe, slope = slope)
    return(res)

  }
}

