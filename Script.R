rm(list = ls())

# Library -----------------------------------------------------------------
library(tidyverse)
source('Supporting/crossvalid.R')
library(ggpubr)
library(ggtext)
library(FactoMineR)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(doParallel)
library(foreach)
library(BGLR)
library(ggnewscale)
library(gghighlight)
library(EnvRtype)

# Data --------------------------------------------------------------------
pheno = readRDS('Data/pheno.rds')
geno = load('Data/geno.rda')

# Diversity ---------------------------------------------------------------
G = tcrossprod(apply(genoF, 2, scale))/ncol(genoF)
rownames(G) = colnames(G) = rownames(genoF)
pca = PCA(G, scale.unit = F)

data.frame(
  line = names(pca$ind$coord[,1]),
  pc1 = pca$ind$coord[,1],
  pc2 = pca$ind$coord[,2],
  cat = ifelse(rownames(pca$ind$coord) %in% unique(pheno$Aline), 'A', 'R'),
  row.names = NULL
) %>% ggplot(aes(x = pc1, y = pc2, color = cat, shape = cat)) + 
  geom_point(size = 2) + 
  labs(x = paste0('PC1 (', round(pca$eig[1,2],2), '%)'),
       y = paste0('PC2 (', round(pca$eig[2,2],2), '%)'),
       color = 'Parental lines', shape = 'Parental lines', 
       fill = 'Parental lines') + 
  theme(legend.position = 'top', legend.text = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold')) + 
  scale_colour_manual(values = c('#1CB80F','#386cb0')) +
  scale_fill_manual(values = c('#1CB80F','#386cb0')) #+ ylim(-6,7.5) + xlim(-6,7.5)

Heatmap(
  G[order(rownames(G)), order(colnames(G))], 
  column_order = sort(colnames(G)), row_order = sort(rownames(G)),
  #col = colorRamp2(c(-.5,0,2), c('#5e3c99', '#f7f7f7', '#e66101')),
  row_split = c(rep("A lines", length(unique(pheno$Aline))), 
                rep("R Lines", length(unique(pheno$Rline)))),
  column_split = c(rep("A lines", length(unique(pheno$Aline))), 
                   rep("R Lines", length(unique(pheno$Rline)))),
  heatmap_legend_param = list(title = 'Kin', legend_height = unit(5, 'cm'),
                              border = 'black'),
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9)
)

rm(pca, G, genoF)

# Genomic kinship matrices -----------------------------------------------
GA = tcrossprod(apply(genoA, 2, scale))/ncol(genoA)
rownames(GA) = colnames(GA) = rownames(genoA)
GA = GA[order(rownames(GA)), order(colnames(GA))]
GR = tcrossprod(apply(genoR, 2, scale))/ncol(genoR)
rownames(GR) = colnames(GR) = rownames(genoR)
GR = GR[order(rownames(GR)), order(colnames(GR))]

rm(genoA, genoR, genoF)

# Environmental data ------------------------------------------------------
envcov = read.csv('Data/envcov.csv')
id.vars = names(envcov)[-which(names(envcov) %in% 
                                 c('env','YEAR', 'MM','DD','DOY','YYYYMMDD',
                                   'daysFromStart','soc','silt','sand','phh2o',
                                   'bdod','alt','clay','LON','LAT','nit'))]

## Omega matrix: mean ------------------------------------------------------
W_mean = W_matrix(env.data = envcov, env.id = 'env', var.id = id.vars)
W_mean = cbind(W_mean, scale(unique(envcov[,c('env', 'soc', 'silt', 'sand', 
                                              'phh2o', 'bdod', 'alt', 'clay', 
                                              'LON', 'LAT', 'nit')])[,-1]))
O_mean = tcrossprod(W_mean)/ncol(W_mean)

## Omega matrix: 4 stages ------------------------------------------------------
W_stages = EnvRtype::W_matrix(env.data = envcov, env.id = 'env', var.id = id.vars,
                              by.interval = T, time.window = c(20,60,95,130), 
                              names.window = c('E0_E1','E2_E5','E6_E7','E8_E9'))
W_stages = W_stages[,-grep('NA', colnames(W_stages))]
W_stages = cbind(W_stages, scale(unique(envcov[,c('env', 'soc', 'silt', 'sand', 
                                                  'phh2o', 'bdod', 'alt', 'clay', 
                                                  'LON', 'LAT', 'nit')])[,-1]))
O_stages = tcrossprod(W_stages)/ncol(W_stages)

## Omega matrix: daily ------------------------------------------------------
W_daily  = EnvRtype::W_matrix(env.data = envcov, env.id = 'env', var.id = id.vars,
                              by.interval = T, time.window = 1:130)
W_daily = cbind(W_daily, scale(unique(envcov[,c('env', 'soc', 'silt', 'sand', 
                                                'phh2o', 'bdod', 'alt', 'clay', 
                                                'LON', 'LAT', 'nit')])[,-1]))
O_daily = tcrossprod(W_daily)/ncol(W_daily)

## Critical time window ------------------------------------------------------
pb = winProgressBar(
  min = 0, max = length(unique(envcov$daysFromStart)), 
  label = 'Completed', width = 300L, initial = 0,
  title = 'Critical time window'
)
corr = list()
fit = list()
for (i in 1:length(unique(envcov$daysFromStart))) {
  
  for (j in 1:length(unique(envcov$daysFromStart))) {
    
    if (length(unique(envcov$daysFromStart)) - i < 6) break
    if (j-i < 6) next
    
    test = list()
    for (k in unique(pheno$env)) {
      test[[k]] = apply(
        as.matrix(envcov[
          which(envcov$daysFromStart %in% i:j & envcov$env == k),
          -which(names(envcov) %in% 
                   c('env','YEAR', 'MM','DD','DOY','YYYYMMDD',
                     'daysFromStart', 'alt', 'LON', 'LAT',
                     'soc','silt','sand','phh2o','nit','clay','bdod'))
        ]), 2, mean
      )  
    }
    
    test = cbind(do.call(rbind,test)[order(rownames(do.call(rbind,test))),],
                 y = tapply(pheno$EBLUE, pheno$env, mean))
    correl = as.matrix(cor(test)[,'y'])
    R2 = as.matrix(apply(test[,-which(colnames(test) == 'y')], 2, function(x){
      summary(lm(test[,'y'] ~ x))$r.squared}))
    colnames(correl) = colnames(R2) = paste(i, j, sep = '-')
    corr[[paste(i, j, sep = '-')]] = correl
    fit[[paste(i, j, sep = '-')]] = R2
    rm(test, correl, R2)
  }
  pctg = paste0(round(i/length(unique(envcov$daysFromStart)) * 100, 0), '% done')
  setWinProgressBar(pb, i, label = pctg)
}
close(pb)

cpt.corr = do.call(cbind, corr)
cpt.corr = as.matrix(apply(cpt.corr[-which(rownames(cpt.corr) == 'y'),], 1, 
                           function(x) names(which.max(x))))

cpt.fit = do.call(cbind, fit)
cpt.fit = as.matrix(apply(cpt.fit, 1, function(x) names(which.max(x))))

## Leave-One-Out cross-validation
pb = winProgressBar(
  min = 0, max = length(unique(pheno$env)), 
  label = 'Progress', width = 300L, initial = 0,
  title = 'Critical time window - LOO'
)
results.fit = list()
results.corr = list()
for (f in unique(pheno$env)) {
  
  loo = envcov[-which(envcov$env == f),]
  means = tapply(pheno[-which(pheno$env == f),'EBLUE'],
                 pheno[-which(pheno$env == f),'env'],
                 mean)
  corr = list()
  fit = list()
  for (i in 1:length(unique(envcov$daysFromStart))) {
    
    for (j in 1:length(unique(envcov$daysFromStart))) {
      
      if (length(unique(envcov$daysFromStart)) - i < 6) break
      if (j-i < 6) next
      
      test = list()
      for (k in unique(pheno$env)[-which(unique(pheno$env) == f)]) {
        test[[k]] = apply(
          as.matrix(loo[
            which(loo$daysFromStart %in% i:j & loo$env == k),
            -which(names(loo) %in% 
                     c('env','YEAR', 'MM','DD','DOY','YYYYMMDD',
                       'daysFromStart', 'alt', 'LON', 'LAT',
                       'soc','silt','sand','phh2o','nit','clay','bdod'))
          ]), 2, mean
        )  
      }
      
      test = cbind(do.call(rbind,test)[order(rownames(do.call(rbind,test))),],
                   y = means)
      correl = as.matrix(cor(test)[,'y'])
      R2 = as.matrix(apply(test[,-which(colnames(test) == 'y')], 2, function(x){
        summary(lm(test[,'y'] ~ x))$r.squared}))
      colnames(correl) = colnames(R2) = paste(i, j, sep = '-')
      corr[[paste(i, j, sep = '-')]] = correl
      fit[[paste(i, j, sep = '-')]] = R2
      rm(test,correl)
    }
  }
  results.corr[[f]] = do.call(cbind, corr)
  results.fit[[f]] = do.call(cbind, fit)
  rm(corr, fit)
  pctg = paste0(round(which(unique(pheno$env) == f)/length(unique(pheno$env)) * 100, 0), '% done')
  setWinProgressBar(pb, f, label = pctg)
}
close(pb)

cpt_loo.corr = do.call(cbind, lapply(results.corr, function(x){
  as.matrix(apply(x[-which(rownames(x) == 'y'),], 1, 
                  function(y) names(which.max(y))))
}))

cpt_loo.fit = do.call(cbind, lapply(results.fit, function(x){
  as.matrix(apply(x, 1, function(y) names(which.max(y))))
}))

ctw.corr = as.matrix(apply(cbind(cpt_loo.corr, cpt.corr), 1, 
                           function(x) names(which.max(table(x)))))
ctw.fit = as.matrix(apply(cbind(cpt_loo.fit, cpt.fit), 1, 
                          function(x) names(which.max(table(x)))))
ctw = cbind(corr = ctw.corr, R2 = ctw.fit)
cpt = data.frame(
  var = rownames(ctw),
  start = as.numeric(sub('-', '',str_extract(ctw[,1], pattern = '\\d+\\-'))),
  end = as.numeric(sub('-', '',str_extract(ctw[,1], pattern = '\\-\\d+'))) 
)

## Optimized omega matrix: mean ------------------------------------------------------
W_opt_mean = matrix(NA, nrow = length(unique(envcov$env)), 
                    ncol = length(unique(cpt$var)),
                    dimnames = list(unique(envcov$env), unique(cpt$var)))

for (i in 1:nrow(cpt)) {
  W_opt_mean[,i] = tapply(
    as.matrix(envcov[which(envcov$daysFromStart >= cpt[i, 'start'] &
                             envcov$daysFromStart <= cpt[i, 'end']), 
                     c('env','daysFromStart',cpt[i, 'var'])][,cpt[i, 'var']]),
    as.matrix(envcov[which(envcov$daysFromStart >= cpt[i, 'start'] &
                             envcov$daysFromStart <= cpt[i, 'end']), 
                     c('env','daysFromStart',cpt[i, 'var'])][,'env']),
    mean
  )
}

W_opt_mean = as.matrix(cbind(W_opt_mean, unique(envcov[
  ,c('env','LON','LAT','alt','soc','silt','sand','phh2o','nit','clay','bdod')
])[,c('LON','LAT','alt','soc','silt','sand','phh2o','nit','clay','bdod')]))
W_opt_mean = apply(W_opt_mean, 2, scale)
O_opt_mean = tcrossprod(W_opt_mean)/ncol(W_opt_mean)

## Most influential environmental covariates --------------------------------
colnames(W_mean) = c('iihs', "hiri", 'etp', 'frue', 'gdd', 'dh', 'sh', 'petp', 'prec',
                     'rh', 'rta', 'spv', 'tmax', 'tmean', 'tmin', 'trange', 'tdew',
                     'vpd', 'ws', 'soc', 'silt', 'sand', 'phh2o', 'bdod', 'alt', 'clay',
                     'lon', 'lat', 'nit')

blupsgei = predict(modgei, classify = 'gen:env', pworkspace = '2gb')$pvals
blupsgei = left_join(blupsgei %>% select(-std.error, -status), 
                     as.data.frame(W_mean) %>% rownames_to_column('env'), by = 'env')

d = lm(
  formula = reformulate(paste(colnames(W_mean)[-c(15,16,9)], collapse = '+'), 
                        response = 'predicted.value'),
  data = blupsgei
)
summary(d)

GEIrel = relaimpo::calc.relimp(d, type = 'lmg', rela = F, rank = T)
GEIrel.df = data.frame(GEIrel@lmg) %>% rownames_to_column('covamb') %>%
  mutate(beta = scale(coef(d))[-1]) 

ggarrange(
  GEIrel.df %>% ggplot(aes(x = reorder(covamb,beta), y = beta, fill = beta)) + 
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 90), legend.position = 'none') + 
    labs(x = '', y = 'Regression coefficient')+
    scale_fill_viridis_c(option = 'D'),
  GEIrel.df %>% ggplot(aes(x = reorder(covamb,GEIrel.lmg), 
                           y = (GEIrel.lmg/GEIrel@R2.decomp)*100, 
                           fill = (GEIrel.lmg/GEIrel@R2.decomp)*100)) + 
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 90), legend.position = 'none') + 
    labs(x = 'Environmental covariates', y = 'Part in R² (%)', 
         caption = paste('R² =', round(GEIrel@R2.decomp, 2))) + 
    scale_fill_viridis_c(option = 'H'),
  nrow = 2,labels = 'auto'
)

# GP Models ---------------------------------------------------------------
## Design matrices ------------------------------------------------------
ZA = model.matrix(~-1 + Aline, data = pheno)
ZR = model.matrix(~-1 + Rline, data = pheno)
ZE = model.matrix(~-1 + env, data = pheno)
ZH = model.matrix(~-1 + gen, data = pheno)

rownames(ZA) = rownames(ZR) = rownames(ZE) = rownames(ZH) = pheno$gen

## Connectivity ------------------------------------------------------------
con = crossprod(ZE,ZH)
colnames(con) = sub('gen','',colnames(con))
rownames(con) = sub('env','',rownames(con))
con = ifelse(con == 1, 'P', 'A')
Heatmap(
  con, col = c('white', '#386cb0'),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 7),
  show_column_names = F, show_row_names = F,
  column_title = 'Hybrids', row_title = 'Environments'
)
table(con)['P']/sum(table(con)) * 100
rm(con)

## Full covariances matrices ------------------------------------------------
ZGZ_A = tcrossprod(tcrossprod(ZA, GA),ZA)
ZGZ_R = tcrossprod(tcrossprod(ZR, GR),ZR)
ZGZ_H = ZGZ_A * ZGZ_R
ZZ_E = tcrossprod(ZE)
GEI_A = ZGZ_A * ZZ_E
GEI_R = ZGZ_R * ZZ_E
GEI_H = ZGZ_H * ZZ_E

ZOZ_E_mean = tcrossprod(tcrossprod(ZE, O_mean), ZE)
GEIO_A_mean = ZGZ_A * ZOZ_E_mean
GEIO_R_mean = ZGZ_R * ZOZ_E_mean
GEIO_H_mean = ZGZ_H * ZOZ_E_mean

ZOZ_E_stages = tcrossprod(tcrossprod(ZE, O_stages), ZE)
GEIO_A_stages = ZGZ_A * ZOZ_E_stages
GEIO_R_stages = ZGZ_R * ZOZ_E_stages
GEIO_H_stages = ZGZ_H * ZOZ_E_stages

ZOZ_E_daily = tcrossprod(tcrossprod(ZE, O_daily), ZE)
GEIO_A_daily = ZGZ_A * ZOZ_E_daily
GEIO_R_daily = ZGZ_R * ZOZ_E_daily
GEIO_H_daily = ZGZ_H * ZOZ_E_daily

ZOZ_E_opt_daily = tcrossprod(tcrossprod(ZE, O_opt_daily), ZE)
GEIO_A_opt_daily = ZGZ_A * ZOZ_E_opt_daily
GEIO_R_opt_daily = ZGZ_R * ZOZ_E_opt_daily
GEIO_H_opt_daily = ZGZ_H * ZOZ_E_opt_daily

ZOZ_E_opt_mean = tcrossprod(tcrossprod(ZE, O_opt_mean), ZE)
GEIO_A_opt_mean = ZGZ_A * ZOZ_E_opt_mean
GEIO_R_opt_mean = ZGZ_R * ZOZ_E_opt_mean
GEIO_H_opt_mean = ZGZ_H * ZOZ_E_opt_mean

## Model 1 - Only GCA -------------------------------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZGZ_A)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS')
)

M1 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m1 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M1$yHat
)

varcomp_m1 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'e'),
  varcomp = c(M1$ETA[[1]]$varU, M1$ETA[[2]]$varU, M1$ETA[[3]]$varU, M1$varE),
  perc = c(M1$ETA[[1]]$varU/sum(c(M1$ETA[[1]]$varU, M1$ETA[[2]]$varU, M1$ETA[[3]]$varU, M1$varE)),
           M1$ETA[[2]]$varU/sum(c(M1$ETA[[1]]$varU, M1$ETA[[2]]$varU, M1$ETA[[3]]$varU, M1$varE)),
           M1$ETA[[3]]$varU/sum(c(M1$ETA[[1]]$varU, M1$ETA[[2]]$varU, M1$ETA[[3]]$varU, M1$varE)),
           M1$varE/sum(c(M1$ETA[[1]]$varU, M1$ETA[[2]]$varU, M1$ETA[[3]]$varU, M1$varE)))
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M1_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M1_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M1_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M1_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 2 - GCA and SCA --------------------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS')
)

M2 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m2 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M2$yHat
)

varcomp_m2 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid','e'),
  varcomp = c(
    M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU, M2$ETA[[4]]$varU,
    M2$varE
  ),
  perc = c(
    M2$ETA[[1]]$varU/sum(c(M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU, 
                           M2$ETA[[4]]$varU, M2$varE)),
    M2$ETA[[2]]$varU/sum(c(M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU, 
                           M2$ETA[[4]]$varU, M2$varE)),
    M2$ETA[[3]]$varU/sum(c(M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU, 
                           M2$ETA[[4]]$varU, M2$varE)),
    M2$ETA[[4]]$varU/sum(c(M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU, 
                           M2$ETA[[4]]$varU, M2$varE)),
    M2$varE/sum(c(M2$ETA[[1]]$varU, M2$ETA[[2]]$varU, M2$ETA[[3]]$varU,
                  M2$ETA[[4]]$varU, M2$varE))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M2_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M2_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M2_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M2_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 3 - GCA, SCA and Interactions -------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_A)$values, model = 'RKHS')
)

M3 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m3 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M3$yHat
)

varcomp_m3 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, M3$ETA[[4]]$varU,
    M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, M3$ETA[[7]]$varU, M3$varE
  ),
  perc = c(
    M3$ETA[[1]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[2]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU,
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[3]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[4]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[5]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[6]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$ETA[[7]]$varU/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU, 
                           M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                           M3$ETA[[7]]$varU, M3$varE)),
    M3$varE/sum(c(M3$ETA[[1]]$varU, M3$ETA[[2]]$varU, M3$ETA[[3]]$varU,
                  M3$ETA[[4]]$varU, M3$ETA[[5]]$varU, M3$ETA[[6]]$varU, 
                  M3$ETA[[7]]$varU, M3$varE))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M3_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M3_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M3_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M3_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 4.1 - GCA, SCA and Interactions w/ EnvCov - mean -------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_mean)$vectors, d = eigen(GEIO_A_mean)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_mean)$vectors, d = eigen(GEIO_R_mean)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_mean)$vectors, d = eigen(GEIO_H_mean)$values, 
            model = 'RKHS')
)

M4.1 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m4.1 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M4.1$yHat
)

varcomp_m4.1 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW','e'),
  varcomp = c(
    M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, M4.1$ETA[[4]]$varU,
    M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, M4.1$ETA[[7]]$varU, M4.1$varE
  ),
  perc = c(
    M4.1$ETA[[1]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[2]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU,
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[3]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[4]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[5]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[6]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$ETA[[7]]$varU/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU, 
                             M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                             M4.1$ETA[[7]]$varU, M4.1$varE)),
    M4.1$varE/sum(c(M4.1$ETA[[1]]$varU, M4.1$ETA[[2]]$varU, M4.1$ETA[[3]]$varU,
                    M4.1$ETA[[4]]$varU, M4.1$ETA[[5]]$varU, M4.1$ETA[[6]]$varU, 
                    M4.1$ETA[[7]]$varU, M4.1$varE))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M4.1_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M4.1_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M4.1_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M4.1_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 4.2 - GCA, SCA and Interactions w/ EnvCov - 4 stages -------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_stages)$vectors, d = eigen(GEIO_A_stages)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_stages)$vectors, d = eigen(GEIO_R_stages)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_stages)$vectors, d = eigen(GEIO_H_stages)$values, 
            model = 'RKHS')
)

M4.2 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m4.2 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M4.2$yHat
)

varcomp_m4.2 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW','e'),
  varcomp = c(
    M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, M4.2$ETA[[4]]$varU,
    M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, M4.2$ETA[[7]]$varU, M4.2$varE
  ),
  perc = c(
    M4.2$ETA[[1]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[2]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU,
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[3]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[4]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[5]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[6]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$ETA[[7]]$varU/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU, 
                             M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                             M4.2$ETA[[7]]$varU, M4.2$varE)),
    M4.2$varE/sum(c(M4.2$ETA[[1]]$varU, M4.2$ETA[[2]]$varU, M4.2$ETA[[3]]$varU,
                    M4.2$ETA[[4]]$varU, M4.2$ETA[[5]]$varU, M4.2$ETA[[6]]$varU, 
                    M4.2$ETA[[7]]$varU, M4.2$varE))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M4.2_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M4.2_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M4.2_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M4.2_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 4.3 - GCA, SCA and Interactions w/ EnvCov - daily -------------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_daily)$vectors, d = eigen(GEIO_A_daily)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_daily)$vectors, d = eigen(GEIO_R_daily)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_daily)$vectors, d = eigen(GEIO_H_daily)$values, 
            model = 'RKHS')
)

M4.3 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m4.3 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M4.3$yHat
)

varcomp_m4.3 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW','e'),
  varcomp = c(
    M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, M4.3$ETA[[4]]$varU,
    M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, M4.3$ETA[[7]]$varU, M4.3$varE
  ),
  perc = c(
    M4.3$ETA[[1]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[2]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU,
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[3]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[4]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[5]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[6]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$ETA[[7]]$varU/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU, 
                             M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                             M4.3$ETA[[7]]$varU, M4.3$varE)),
    M4.3$varE/sum(c(M4.3$ETA[[1]]$varU, M4.3$ETA[[2]]$varU, M4.3$ETA[[3]]$varU,
                    M4.3$ETA[[4]]$varU, M4.3$ETA[[5]]$varU, M4.3$ETA[[6]]$varU, 
                    M4.3$ETA[[7]]$varU, M4.3$varE))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M4.3_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M4.3_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M4.3_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M4.3_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 5.1 - GCA, SCA and Interactions w/ EnvCov and without Envcov - mean ----
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_mean)$vectors, d = eigen(GEIO_A_mean)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_mean)$vectors, d = eigen(GEIO_R_mean)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_mean)$vectors, d = eigen(GEIO_H_mean)$values, 
            model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_R)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_R)$values, model = 'RKHS')
)

M5.1 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m5.1 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M5.1$yHat
)

varcomp_m5.1 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW',
             'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, M5.1$ETA[[4]]$varU,
    M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU,
    M5.1$ETA[[9]]$varU, M5.1$ETA[[10]]$varU, M5.1$varE
  ),
  perc = c(
    M5.1$ETA[[1]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[2]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[3]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[4]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[5]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[6]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[7]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[8]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[9]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                             M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                             M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                             M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$ETA[[10]]$varU/sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                              M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                              M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                              M5.1$ETA[[10]]$varU, M5.1$varE)),
    M5.1$varE/sum(sum(c(M5.1$ETA[[1]]$varU, M5.1$ETA[[2]]$varU, M5.1$ETA[[3]]$varU, 
                        M5.1$ETA[[4]]$varU, M5.1$ETA[[5]]$varU, M5.1$ETA[[6]]$varU, 
                        M5.1$ETA[[7]]$varU, M5.1$ETA[[8]]$varU, M5.1$ETA[[9]]$varU,
                        M5.1$ETA[[10]]$varU, M5.1$varE)))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M5.1_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M5.1_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M5.1_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M5.1_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 5.2 - GCA, SCA and Interactions w/ EnvCov and without Envcov - 4 stages ----
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_stages)$vectors, d = eigen(GEIO_A_stages)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_stages)$vectors, d = eigen(GEIO_R_stages)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_stages)$vectors, d = eigen(GEIO_H_stages)$values, 
            model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_R)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_R)$values, model = 'RKHS')
)

M5.2 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m5.2 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M5.2$yHat
)

varcomp_m5.2 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW',
             'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, M5.2$ETA[[4]]$varU,
    M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU,
    M5.2$ETA[[9]]$varU, M5.2$ETA[[10]]$varU, M5.2$varE
  ),
  perc = c(
    M5.2$ETA[[1]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[2]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[3]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[4]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[5]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[6]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[7]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[8]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[9]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                             M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                             M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                             M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$ETA[[10]]$varU/sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                              M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                              M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                              M5.2$ETA[[10]]$varU, M5.2$varE)),
    M5.2$varE/sum(sum(c(M5.2$ETA[[1]]$varU, M5.2$ETA[[2]]$varU, M5.2$ETA[[3]]$varU, 
                        M5.2$ETA[[4]]$varU, M5.2$ETA[[5]]$varU, M5.2$ETA[[6]]$varU, 
                        M5.2$ETA[[7]]$varU, M5.2$ETA[[8]]$varU, M5.2$ETA[[9]]$varU,
                        M5.2$ETA[[10]]$varU, M5.2$varE)))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M5.2_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M5.2_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M5.2_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M5.2_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 5.3 - GCA, SCA and Interactions w/ EnvCov and without Envcov - daily ----
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_daily)$vectors, d = eigen(GEIO_A_daily)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_daily)$vectors, d = eigen(GEIO_R_daily)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_daily)$vectors, d = eigen(GEIO_H_daily)$values, 
            model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_R)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_R)$values, model = 'RKHS')
)

M5.3 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m5.3 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M5.3$yHat
)

varcomp_m5.3 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW',
             'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, M5.3$ETA[[4]]$varU,
    M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU,
    M5.3$ETA[[9]]$varU, M5.3$ETA[[10]]$varU, M5.3$varE
  ),
  perc = c(
    M5.3$ETA[[1]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[2]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[3]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[4]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[5]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[6]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[7]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[8]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[9]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                             M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                             M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                             M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$ETA[[10]]$varU/sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                              M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                              M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                              M5.3$ETA[[10]]$varU, M5.3$varE)),
    M5.3$varE/sum(sum(c(M5.3$ETA[[1]]$varU, M5.3$ETA[[2]]$varU, M5.3$ETA[[3]]$varU, 
                        M5.3$ETA[[4]]$varU, M5.3$ETA[[5]]$varU, M5.3$ETA[[6]]$varU, 
                        M5.3$ETA[[7]]$varU, M5.3$ETA[[8]]$varU, M5.3$ETA[[9]]$varU,
                        M5.3$ETA[[10]]$varU, M5.3$varE)))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M5.3_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M5.3_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M5.3_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M5.3_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 6 - Same as model 5 with crucial periods of time - daily ----------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_opt_daily)$vectors, d = eigen(GEIO_A_opt_daily)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_opt_daily)$vectors, d = eigen(GEIO_R_opt_daily)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_opt_daily)$vectors, d = eigen(GEIO_H_opt_daily)$values, 
            model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_R)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_R)$values, model = 'RKHS')
)

M6 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m6 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M6$yHat
)

varcomp_m6 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW',
             'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, M6$ETA[[4]]$varU,
    M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, M6$ETA[[7]]$varU, M6$ETA[[8]]$varU,
    M6$ETA[[9]]$varU, M6$ETA[[10]]$varU, M6$varE
  ),
  perc = c(
    M6$ETA[[1]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[2]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[3]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[4]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[5]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[6]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[7]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[8]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[9]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                           M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                           M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                           M6$ETA[[10]]$varU, M6$varE)),
    M6$ETA[[10]]$varU/sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                            M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                            M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                            M6$ETA[[10]]$varU, M6$varE)),
    M6$varE/sum(sum(c(M6$ETA[[1]]$varU, M6$ETA[[2]]$varU, M6$ETA[[3]]$varU, 
                      M6$ETA[[4]]$varU, M6$ETA[[5]]$varU, M6$ETA[[6]]$varU, 
                      M6$ETA[[7]]$varU, M6$ETA[[8]]$varU, M6$ETA[[9]]$varU,
                      M6$ETA[[10]]$varU, M6$varE)))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M6_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M6_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M6_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M6_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

## Model 7 - Same as model 5 with crucial periods of time - mean ----------
ETA = list(
  E = list(V = eigen(ZZ_E)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  A = list(V = eigen(ZGZ_A)$vectors, d = eigen(ZZ_E)$values, model = 'RKHS'),
  R = list(V = eigen(ZGZ_R)$vectors, d = eigen(ZGZ_R)$values, model = 'RKHS'),
  H = list(V = eigen(ZGZ_H)$vectors, d = eigen(ZGZ_H)$values, model = 'RKHS'),
  AW = list(V = eigen(GEIO_A_opt_mean)$vectors, d = eigen(GEIO_A_opt_mean)$values, 
            model = 'RKHS'),
  RW = list(V = eigen(GEIO_R_opt_mean)$vectors, d = eigen(GEIO_R_opt_mean)$values, 
            model = 'RKHS'),
  HW = list(V = eigen(GEIO_H_opt_mean)$vectors, d = eigen(GEIO_H_opt_mean)$values, 
            model = 'RKHS'),
  AE = list(V = eigen(GEI_A)$vectors, d = eigen(GEI_A)$values, model = 'RKHS'),
  RE = list(V = eigen(GEI_R)$vectors, d = eigen(GEI_R)$values, model = 'RKHS'),
  HE = list(V = eigen(GEI_H)$vectors, d = eigen(GEI_R)$values, model = 'RKHS')
)

M7 = BGLR(y = pheno$EBLUE, ETA = ETA, nIter = 12000, burnIn = 2000)
unlink("*.dat")

pred_m7 = data.frame(
  A_line = pheno$Aline,
  R_line = pheno$Rline, 
  gen = pheno$gen,
  env = pheno$env,
  EBLUE = pheno$EBLUE,
  GEBV = M7$yHat
)

varcomp_m7 = data.frame(
  effect = c("E", 'Aline', 'Rline', 'Hybrid', 'AxW', 'RxW', 'HxW',
             'AxE', 'RxE', 'HxE','e'),
  varcomp = c(
    M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, M7$ETA[[4]]$varU,
    M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, M7$ETA[[7]]$varU, M7$ETA[[8]]$varU,
    M7$ETA[[9]]$varU, M7$ETA[[10]]$varU, M7$varE
  ),
  perc = c(
    M7$ETA[[1]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[2]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[3]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[4]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[5]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[6]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[7]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[8]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[9]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                           M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                           M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                           M7$ETA[[10]]$varU, M7$varE)),
    M7$ETA[[10]]$varU/sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                            M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                            M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                            M7$ETA[[10]]$varU, M7$varE)),
    M7$varE/sum(sum(c(M7$ETA[[1]]$varU, M7$ETA[[2]]$varU, M7$ETA[[3]]$varU, 
                      M7$ETA[[4]]$varU, M7$ETA[[5]]$varU, M7$ETA[[6]]$varU, 
                      M7$ETA[[7]]$varU, M7$ETA[[8]]$varU, M7$ETA[[9]]$varU,
                      M7$ETA[[10]]$varU, M7$varE)))
  )
)

### CV 1 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 5, nrept = 10, cv = 'cv1', niter = 12000,
                     burnin = 2000)

M7_CV1 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 2 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = 2, nrept = 10, cv = 'cv2', niter = 12000,
                     burnin = 2000)

M7_CV2 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 0 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv0', 
                     niter = 12000, burnin = 2000)

M7_CV0 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)

### CV 00 --------------------------------------------------------------------
cvdata = assignation(data = pheno, y = 'EBLUE', gen = 'gen', env = 'env', 
                     seed = 1000, nfolds = NA, nrept = NA, cv = 'loo_cv00', 
                     niter = 12000, burnin = 2000)

M7_CV00 = crossvalid(cvdata = cvdata, ETA = ETA, results = NULL)


# Variance components -----------------------------------------------------
VC = list(
  M1 = varcomp_m1,
  M2 = varcomp_m2,
  M3 = varcomp_m3,
  M4.1 = varcomp_m4.1,
  M4.2 = varcomp_m4.2,
  M4.3 = varcomp_m4.3,
  M5.1 = varcomp_m5.1,
  M5.2 = varcomp_m5.2,
  M5.3 = varcomp_m5.3,
  M6 = varcomp_m6,
  M7 = varcomp_m7
)

## Plot --------------------------------------------------------------------

names(VC)[names(VC) == 'M6'] = 'M6.1'
names(VC)[names(VC) == 'M7'] = 'M6.2'

for (i in names(VC)) {
  VC[[i]]$Model = i
}

do.call(rbind, VC) |> 
  ggplot(aes(x = Model, y = perc,
             fill = factor(effect, levels = c(
               'AxW', 'RxW', 'HxW', 'AxE', 'RxE','HxE',
               'Aline', 'Rline', 'Hybrid', 'E', 'e'
             )))) + 
  geom_bar(stat = 'identity') + 
  labs(y = 'Percentage of total variance', fill = 'Variance \n component') + 
  scale_fill_manual(
    values = c(
      'E' = '#1f78b4',
      'Aline' = '#a6cee3', 
      'Rline' = '#b2df8a',
      'Hybrid' = '#33a02c',
      'AxW' = '#fb9a99',
      'RxW' = '#e31a1c',
      'HxW' = '#fdbf6f',
      'AxE' = '#ff7f00', 
      'RxE' = '#cab2d6',
      'HxE' = '#6a3d9a',
      'e' = '#ffff99'
    ),
    label = c(
      'E' = expression(sigma[e]^2),
      'Aline' = expression(sigma[a]^2), 
      'Rline' = expression(sigma[r]^2),
      'Hybrid' = expression(sigma[h]^2),
      'AxW' = expression(sigma[a~omega]^2),
      'RxW' = expression(sigma[r~omega]^2),
      'HxW' = expression(sigma[h~omega]^2),
      'AxE' = expression(sigma[ae]^2), 
      'RxE' = expression(sigma[re]^2),
      'HxE' = expression(sigma[he]^2),
      'e' = expression(sigma[epsilon]^2)
    )
  )+ 
  theme(legend.position = 'right')

# CV results --------------------------------------------------------------

## CV1 ---------------------------------------------------------------------
CV1 = do.call(list,lapply(ls()[grep('CV1', ls())], as.name))
names(CV1) = sub('_CV1','',ls()[grep('CV1', ls())])[-1]

CV1 = unlist(lapply(CV1, function(x){
  lapply(x, function(y){
    do.call(rbind, y)
  })
}), recursive = F)
CV1.yhat = CV1[grep('yhat', names(CV1))]
CV1 = CV1[-grep('yhat', names(CV1))]

a = do.call(rbind,lapply(CV1[grep('corr_env', names(CV1))], reshape2::melt))
colnames(a) = c("rept", 'env', 'corr_env')
a$Model = rep(sub('_CV1','',ls()[grep('CV1$', ls())])[-1], 
              each = length(unique(a$env)) * 10)
rownames(a) = NULL

CV1 = left_join(a, data.frame(
  Model = rep(sub('_CV1','',ls()[grep('CV1$', ls())])[-1], each = 10),
  rept = rep(paste0('R',1:10), times = 11),
  corr = c(do.call(cbind, CV1[grep('corr\\b', names(CV1))])),
  wapa = c(do.call(cbind, CV1[grep('wapa', names(CV1))])),
  mspe = c(do.call(cbind, CV1[grep('mspe', names(CV1))])),
  slope = c(do.call(cbind, CV1[grep('slope', names(CV1))]))
), by = c('Model', 'rept')
) |> relocate(corr_env, .after = corr)

## CV2 ---------------------------------------------------------------------
CV2 = do.call(list,lapply(ls()[grep('CV2', ls())], as.name))
names(CV2) = sub('_CV2','',ls()[grep('CV2', ls())])[-1]

CV2 = unlist(lapply(CV2, function(x){
  lapply(x, function(y){
    do.call(rbind, y)
  })
}), recursive = F)
CV2.yhat = CV2[grep('yhat', names(CV2))]
CV2 = CV2[-grep('yhat', names(CV2))]

a = do.call(rbind,lapply(CV2[grep('corr_env', names(CV2))], reshape2::melt))
colnames(a) = c("rept", 'env', 'corr_env')
a$Model = rep(sub('_CV2','',ls()[grep('CV2$', ls())])[-1], 
              each = length(unique(a$env)) * 10)
rownames(a) = NULL

CV2 = left_join(a, data.frame(
  Model = rep(sub('_CV2','',ls()[grep('CV2$', ls())])[-1], each = 10),
  rept = rep(paste0('R',1:10), times = 11),
  corr = c(do.call(cbind, CV2[grep('corr\\b', names(CV2))])),
  wapa = c(do.call(cbind, CV2[grep('wapa', names(CV2))])),
  mspe = c(do.call(cbind, CV2[grep('mspe', names(CV2))])),
  slope = c(do.call(cbind, CV2[grep('slope', names(CV2))]))
), by = c('Model', 'rept')
) |> relocate(corr_env, .after = corr)

## CV0 ---------------------------------------------------------------------
CV0 = do.call(list,lapply(ls()[grep('CV0\\b', ls())], as.name))
names(CV0) = sub('_CV0','',ls()[grep('CV0\\b', ls())])[-1]

CV0 = unlist(CV0, recursive = F)
CV0.yhat = CV0[grep('yhat', names(CV0))]
CV0 = CV0[-grep('yhat', names(CV0))]

a = do.call(rbind, lapply(
  lapply(CV0[grep('corr_env', names(CV0))], reshape2::melt),
  function(x) cbind(x, env = rownames(x))
))
a$Model = rep(sub('_CV0','',ls()[grep('CV0$', ls())])[-1], 
              each = length(unique(a$env)))
rownames(a) = NULL
colnames(a)[colnames(a) == 'value'] = 'corr_env'

CV0 = left_join(a, data.frame(
  Model = sub('_CV0','',ls()[grep('CV0$', ls())])[-1],
  corr = c(do.call(cbind, CV0[grep('corr\\b', names(CV0))])),
  wapa = c(do.call(cbind, CV0[grep('wapa', names(CV0))])),
  mspe = c(do.call(cbind, CV0[grep('mspe', names(CV0))])),
  slope = c(do.call(cbind, CV0[grep('slope', names(CV0))]))
), by = c('Model')
) |> mutate(rept = 'R1') |>  
  relocate(corr_env, .after = corr) |> relocate(rept, .before = env)

## CV00 ---------------------------------------------------------------------
CV00 = do.call(list,lapply(ls()[grep('CV00\\b', ls())], as.name))
names(CV00) = sub('_CV00','',ls()[grep('CV00\\b', ls())])[-1]

CV00 = unlist(CV00, recursive = F)
CV00.yhat = CV00[grep('yhat', names(CV00))]
CV00 = CV00[-grep('yhat', names(CV00))]

a = do.call(rbind, lapply(
  lapply(CV00[grep('corr_env', names(CV00))], reshape2::melt),
  function(x) cbind(x, env = rownames(x))
))
a$Model = rep(sub('_CV00','',ls()[grep('CV00$', ls())])[-1], 
              each = length(unique(a$env)))
rownames(a) = NULL
colnames(a)[colnames(a) == 'value'] = 'corr_env'

CV00 = left_join(a, data.frame(
  Model = sub('_CV00','',ls()[grep('CV00$', ls())])[-1],
  corr = c(do.call(cbind, CV00[grep('corr\\b', names(CV00))])),
  wapa = c(do.call(cbind, CV00[grep('wapa', names(CV00))])),
  mspe = c(do.call(cbind, CV00[grep('mspe', names(CV00))])),
  slope = c(do.call(cbind, CV00[grep('slope', names(CV00))]))
), by = c('Model')
) |> mutate(rept = 'R1') |>  
  relocate(corr_env, .after = corr) |> relocate(rept, .before = env)

## Compilation -------------------------------------------------------------
CV = rbind(CV1, CV2, CV0, CV00)
CV$cv = c(rep(c('CV1', 'CV2'), each = dim(CV1)[1]), 
          rep(c("CV0", 'CV00'), each = dim(CV0)[1]))
CV[CV$Model == 'M6','Model'] = 'M6.1'
CV[CV$Model == 'M7','Model'] = 'M6.2'

rm(CV1, CV2, CV0, CV00)

## Plots -------------------------------------------------------------------

# Box-plots: CV1 and CV2
ggarrange(subset(CV, subset = cv %in% c("CV1", 'CV2')) %>%
            ggplot(aes(x = Model, y = wapa, fill = cv)) +
            geom_boxplot() + 
            theme(legend.position = 'top') + 
            labs(y = 'Mean weighted average preditive ability') +
            scale_fill_manual("Cross-validation scheme", 
                              values = c('CV1' = "#ffff99", 'CV2' = '#386cb0')),
          subset(CV, subset = cv %in% c("CV1", 'CV2')) %>%
            ggplot(aes(x = Model, y = mspe, fill = cv)) +
            geom_boxplot() + 
            theme(legend.position = 'top') + 
            labs(y = 'Mean square prediction error') +
            scale_fill_manual("Cross-validation scheme", 
                              values = c('CV1' = "#ffff99", 'CV2' = '#386cb0')),
          common.legend = T, labels = 'auto', ncol = 2)

# Bar plots: CV1, CV2, CV0 and CV00
ggarrange(
  CV |> reframe(CI = 1.96 * sd(wapa), wapa = mean(wapa), .by = c(Model, cv)) |>
    mutate(cv = factor(cv, levels = c("CV1", 'CV2', 'CV0', 'CV00'))) |> 
    ggplot(aes(x = Model, y = wapa)) +
    geom_bar(stat = 'identity', position = position_dodge(), aes(fill = cv)) +
    geom_errorbar(aes(ymin = wapa - CI, ymax = wapa +CI), 
                  width = .2) + 
    geom_text(aes(label = round(wapa,2), y = 0), vjust = -1, fontface = 'bold',
              size = 3) + 
    labs(y = 'Mean weighted average preditive ability') + 
    facet_wrap(~cv) + 
    scale_fill_manual("Cross-validation scheme", 
                      values = c('CV1' = "#ffff99", 'CV2' = '#386cb0',
                                 'CV0' = '#d7191c', 'CV00' = '#abdda4')) + 
    theme(legend.position = 'none',axis.text.x = element_text(angle = 90)) + 
    ylim(0,1),
  
  CV |> reframe(CI = 1.96 * sd(mspe), mspe = mean(mspe), .by = c(Model, cv)) |>  
    mutate(cv = factor(cv, levels = c("CV1", 'CV2', 'CV0', 'CV00'))) |> 
    ggplot(aes(x = Model, y = mspe)) +
    geom_bar(stat = 'identity', position = position_dodge(), aes(fill = cv)) +
    geom_errorbar(aes(ymin = mspe - CI, ymax = mspe +CI), 
                  width = .2) + 
    geom_text(aes(label = round(mspe,2), y = 0), vjust = -1, fontface = 'bold',
              size = 2) + 
    labs(y = 'Mean square prediction error') + 
    facet_wrap(~cv) + 
    scale_fill_manual("Cross-validation scheme", 
                      values = c('CV1' = "#ffff99", 'CV2' = '#386cb0',
                                 'CV0' = '#d7191c', 'CV00' = '#abdda4')) + 
    theme(legend.position = 'none', axis.text.x = element_text(angle = 90)),
  labels = 'auto', ncol = 1
)

# Within-environment predictive ability

heat = as.matrix(CV |> reframe(corr_env = mean(corr_env),
                               .by = c(cv, env, Model)) |> 
                   pivot_wider(names_from = env, values_from = corr_env) |> 
                   mutate(cv = factor(cv, levels = c('CV1', 'CV2', 'CV0', 'CV00'))) |> 
                   arrange(cv) |> 
                   mutate(mod_cv = paste(Model, cv, sep = '_')) |> 
                   column_to_rownames('mod_cv') |> select(-cv, -Model))

Heatmap(
  t(heat), 
  show_row_dend = F, show_column_dend = F, 
  col = colorRamp2(c(-1,-.5,0,.5,1), c("red4","red1","white",
                                       "darkolivegreen2","darkolivegreen4")),
  column_order = rownames(heat), row_order = unique(CV$env),
  column_labels = do.call(rbind, str_split(rownames(heat),'_'))[,1],
  row_names_gp = gpar(fontsize = 7), column_names_rot = 90,
  column_split = rep(c('CV1', 'CV2', 'CV0', 'CV00'), each = 11)
)

## Top 5 coincidence within environments ------------------------
### CV1 ---------------------------------------------------------
ranks.CV1 = lapply(CV1.yhat, function(x){
  x |> reframe(trait = mean(trait), yhat = mean(yhat),
                 .by = c(gen, env)) |> 
    arrange(env) |> slice_max(trait, n = 5, by = env) |> 
    rename(gen.tr = gen) |> 
    mutate(
      gen.yhat = x |> reframe(trait = mean(trait), yhat = mean(yhat),
                                .by = c(gen, env)) |> 
        arrange(env) |> slice_max(yhat, n = 5, by = env) |> select(gen)
    ) |> select(env, gen.tr, gen.yhat)
})

ranks.CV1 = as.matrix(cbind(ranks.CV1$M1.yhat[,-3], 
                            do.call(cbind, lapply(ranks.CV1, function(x) x[3]))))

colnames(ranks.CV1) = c('env', 'gen.tr', c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                                           'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                                           'M6.1', 'M6.2'))


coin.mat.CV1 = matrix(nrow = length(unique(ranks.CV1[,'env'])), 
                      ncol = length(names(CV1.yhat)), 
                      dimnames = list(unique(ranks.CV1[,'env']), 
                                      colnames(ranks.CV1)[grep('M', colnames(ranks.CV1))]))

for (i in rownames(coin.mat.CV1)) {
  for (j in colnames(coin.mat.CV1)) {
    coin.mat.CV1[i,j] = mean(ranks.CV1[which(ranks.CV1[,'env'] == i),'gen.tr'] %in%
                               ranks.CV1[which(ranks.CV1[,'env'] == i), j]) * 100
  }
}

colnames(coin.mat.CV1) = paste(colnames(coin.mat.CV1), 'CV1', sep='_')


### CV2 -------------------------------------------------------
ranks.CV2 = lapply(CV2.yhat, function(x){
  x |> reframe(trait = mean(trait), yhat = mean(yhat),
                 .by = c(gen, env)) |> 
    arrange(env) |> slice_max(trait, n = 5, by = env) |> 
    rename(gen.tr = gen) |> 
    mutate(
      gen.yhat = x |> reframe(trait = mean(trait), yhat = mean(yhat),
                                .by = c(gen, env)) |> 
        arrange(env) |> slice_max(yhat, n = 5, by = env) |> select(gen)
    ) |> select(env, gen.tr, gen.yhat)
})

ranks.CV2 = as.matrix(cbind(ranks.CV2$M1.yhat[,-3], 
                            do.call(cbind, lapply(ranks.CV2, function(x) x[3]))))

colnames(ranks.CV2) = c('env', 'gen.tr', c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                                           'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                                           'M6.1', 'M6.2'))


coin.mat.CV2 = matrix(nrow = length(unique(ranks.CV2[,'env'])), 
                      ncol = length(names(CV2.yhat)), 
                      dimnames = list(unique(ranks.CV2[,'env']), 
                                      colnames(ranks.CV2)[grep('M', colnames(ranks.CV2))]))

for (i in rownames(coin.mat.CV2)) {
  for (j in colnames(coin.mat.CV2)) {
    coin.mat.CV2[i,j] = mean(ranks.CV2[which(ranks.CV2[,'env'] == i),'gen.tr'] %in%
                               ranks.CV2[which(ranks.CV2[,'env'] == i), j]) * 100
  }
}

colnames(coin.mat.CV2) = paste(colnames(coin.mat.CV2), 'CV2', sep='_')


### CV0 -------------------------------------------------------
ranks.CV0 = lapply(CV0.yhat, function(x){
  x |> arrange(env) |> slice_max(trait, n = 5, by = env) |> 
    rename(gen.tr = gen) |> 
    mutate(
      gen.yhat = x |> arrange(env) |> slice_max(yhat, n = 5, by = env) |> select(gen)
    ) |> select(env, gen.tr, gen.yhat)
})

ranks.CV0 = as.matrix(cbind(ranks.CV0$M1.yhat[,-3], 
                            do.call(cbind, lapply(ranks.CV0, function(x) x[3]))))

colnames(ranks.CV0) = c('env', 'gen.tr', c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                                           'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                                           'M6.1', 'M6.2'))


coin.mat.CV0 = matrix(nrow = length(unique(ranks.CV0[,'env'])), 
                      ncol = length(names(CV0.yhat)), 
                      dimnames = list(unique(ranks.CV0[,'env']), 
                                      colnames(ranks.CV0)[grep('M', colnames(ranks.CV0))]))

for (i in rownames(coin.mat.CV0)) {
  for (j in colnames(coin.mat.CV0)) {
    coin.mat.CV0[i,j] = mean(ranks.CV0[which(ranks.CV0[,'env'] == i),'gen.tr'] %in%
                               ranks.CV0[which(ranks.CV0[,'env'] == i), j]) * 100
  }
}

colnames(coin.mat.CV0) = paste(colnames(coin.mat.CV0), 'CV0', sep='_')

### CV00 -------------------------------------------------------
ranks.CV00 = lapply(CV00.yhat, function(x){
  x |> arrange(env) |> slice_max(trait, n = 5, by = env) |> 
    rename(gen.tr = gen) |> 
    mutate(
      gen.yhat = x |> arrange(env) |> slice_max(yhat, n = 5, by = env) |> select(gen)
    ) |> select(env, gen.tr, gen.yhat)
})

ranks.CV00 = as.matrix(cbind(ranks.CV00$M1.yhat[,-3], 
                             do.call(cbind, lapply(ranks.CV00, function(x) x[3]))))

colnames(ranks.CV00) = c('env', 'gen.tr', c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                                            'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                                            'M6.1', 'M6.2'))


coin.mat.CV00 = matrix(nrow = length(unique(ranks.CV00[,'env'])), 
                       ncol = length(names(CV00.yhat)), 
                       dimnames = list(unique(ranks.CV00[,'env']), 
                                       colnames(ranks.CV00)[grep('M', colnames(ranks.CV00))]))

for (i in rownames(coin.mat.CV00)) {
  for (j in colnames(coin.mat.CV00)) {
    coin.mat.CV00[i,j] = mean(ranks.CV00[which(ranks.CV00[,'env'] == i),'gen.tr'] %in%
                                ranks.CV00[which(ranks.CV00[,'env'] == i), j]) * 100
  }
}

colnames(coin.mat.CV00) = paste(colnames(coin.mat.CV00), 'CV00', sep='_')

### Plot -------------------------------------------------------
heat = cbind(coin.mat.CV1, coin.mat.CV2, coin.mat.CV0, coin.mat.CV00)

Heatmap(
  heat, 
  show_row_dend = F, show_column_dend = F, 
  col = colorRamp2(c(0,20,60,100), c("red4","red1", "darkolivegreen2","darkolivegreen4")),
  column_order = colnames(heat), row_order = rownames(heat),
  column_labels = do.call(rbind, str_split(colnames(heat),'_'))[,1],
  row_names_gp = gpar(fontsize = 7), column_names_rot = 90,
  column_split = rep(c('CV1', 'CV2', 'CV0', 'CV00'), each = 11),
  heatmap_legend_param = list(title = 'Coincidence (%)')
)

## Coincidence across environments -----------------------------------------
CV1.yhat = lapply(CV1.yhat, function(x){
  x |> reframe(yhat = mean(yhat), .by = c(gen, env, trait))
})
CV1.yhat = data.frame(do.call(rbind, CV1.yhat)) |> 
  mutate(model = rep(c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                       'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                       'M6.1', 'M6.2'), each = nrow(pheno)), CV = 'CV1')

CV2.yhat = lapply(CV2.yhat, function(x){
  x |> reframe(yhat = mean(yhat), .by = c(gen, env, trait))
})

CV2.yhat = data.frame(do.call(rbind, CV2.yhat)) |> 
  mutate(model = rep(c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                       'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                       'M6.1', 'M6.2'), each = nrow(pheno)),CV = 'CV2')


CV0.yhat = data.frame(do.call(rbind, CV0.yhat)) |> 
  mutate(model = rep(c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                       'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                       'M6.1', 'M6.2'), each = nrow(pheno)), CV = 'CV0')

CV00.yhat = data.frame(do.call(rbind, CV00.yhat)) |> 
  mutate(model = rep(c('M1', 'M2', 'M3', 'M4.1', 'M4.2', 
                       'M4.3', 'M5.1', 'M5.2', 'M5.3', 
                       'M6.1', 'M6.2'), each = nrow(pheno)), CV = 'CV00')

CV.yhat = rbind(CV1.yhat, CV2.yhat, CV0.yhat, CV00.yhat)
rm(CV1.yhat, CV2.yhat, CV0.yhat, CV00.yhat)

plot.list = list()
for (i in apply(expand.grid(unique(CV.yhat$CV), c('M1', 'M5.1', 'M6.2')), 1, 
                function(x){
                  paste(x, collapse = '_')
                })) {
  test = CV.yhat |> filter(CV == str_extract(i, 'CV\\d+') & 
                             model == gsub('CV\\d+_', '', i)) 
  
  plot.list[[i]] = test |> 
    ggplot(aes(x = yhat, y = trait)) + 
    geom_rug() + 
    geom_hline(
      data = test |> reframe(qq = quantile(trait, probs = c(.2, .5, .75 ,.90))),
      aes(yintercept = qq, color = c('20%', '50%', '75%', '90%')), 
      linetype = 'dashed',
      linewidth = 1.2
    ) + 
    geom_vline(
      data = test |> reframe(qq = quantile(yhat, probs = c(.2, .5, .75 ,.90))),
      aes(xintercept = qq, color = c('20%', '50%', '75%', '90%')),
      linetype = 'dashed', 
      linewidth = 1.2
    ) + 
    scale_colour_manual('Quantiles',
                        values = c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a')) + 
    new_scale_colour()+
    geom_point(aes(colour = gen), size = 1.7, show.legend = F) + 
    scale_colour_viridis_d('Coincident', option = 'inferno')  +
    gghighlight(trait >= quantile(test$trait, .9) & yhat >= quantile(test$yhat, .9),
                unhighlighted_params = list(colour = 'darkgrey', size = 1), 
                use_direct_label = F) +
    labs(x = 'Predicted', y = 'Observed', 
         title = paste(paste(gsub('CV\\d+_', '', i), str_extract(i, 'CV\\d+'), 
                             sep = ': '),
                       paste0(round(sum(test$trait >= quantile(test$trait, .9) & 
                                          test$yhat >= quantile(test$yhat, .9)) / 
                                      sum(test$yhat >= quantile(test$yhat, .9)) * 
                                      100, 2), '%'),
                       sep = ' - '))
}

do.call(ggarrange, c(plot.list, common.legend = T, nrow = 3, ncol = 4))


