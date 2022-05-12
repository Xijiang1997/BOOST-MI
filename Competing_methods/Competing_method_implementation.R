

# load the data
load('bc_pattern_zero_10_replicate_1.RData')


# BOOST-MI
##############
source("functions/Boost_Ising_function.R")
    
re_BOOST_MI <- Boost_Ising(count, loc)

# BOOST-GP
#############
# Need to download R file from https://github.com/Minzhe/BOOST-GP/tree/main/R
source('R/boost.gp.R')
re_BOOST_GP <- boost.gp(count, loc, 1600, 1000, size.factor = rowSums(count)/mean(rowSums(count)))
    

# SPARK
#############
library(SPARK)

loc_spark <- data.frame(loc)
si <- apply(count, 1, sum)
count_spark <- t(count)
count_spark <- data.frame(count_spark)
colnames(count_spark) <- rownames(loc_spark)
    
spark <- CreateSPARKObject(counts=count_spark,
                               location=loc_spark,
                               percentage = 0,
                               min_total_counts = 0)
spark@lib_size <- si
spark <- spark.vc(spark,
                      covariates = NULL,
                      lib_size = spark@lib_size,
                      num_core = 5,
                      verbose = F)
spark <- spark.test(spark,
                        check_positive = T,
                        verbose = T)
re_SPARK <- spark@res_mtest
   

# BinSpect
###############
library(Giotto)

# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

n <- dim(loc)[1]
rownames(count) <- as.character(1:n)
colnames(count) <- paste0('X', 1:dim(count)[2])

loc <- data.frame(loc)
loc$ID <- as.character(1:n)

# run Binspect
SS_simu <- createGiottoObject(raw_exprs = t(count), spatial_locs = loc)
SS_simu <- normalizeGiotto(gobject = SS_simu, norm_method = 'standard', scalefactor = 6000, verbose = T)
SS_simu <- createSpatialNetwork(gobject = SS_simu, method = 'kNN', k = 5, name = 'spatial_network')
re_binSpect_rank <- binSpect(SS_simu, bin_method = "rank", expression_values = "normalized",subset_genes = NULL, spatial_network_name = "spatial_network")

# run Binspect
SS_simu <- createGiottoObject(raw_exprs = t(count), spatial_locs = loc)
SS_simu <- normalizeGiotto(gobject = SS_simu, norm_method = 'standard', scalefactor = 6000, verbose = T)
SS_simu <- createSpatialNetwork(gobject = SS_simu, method = 'kNN', k = 5, name = 'spatial_network')
re_binSpect_kmeans <- binSpect(SS_simu, bin_method = "kmeans", expression_values = "normalized",subset_genes = NULL, spatial_network_name = "spatial_network")


# SpatialDE
##################
library(spatialDE)

n <- dim(loc)[1]
rownames(count) <- as.character(1:n)
colnames(count) <- paste0('X', 1:dim(count)[2])
count <- t(count)
loc <- data.frame(loc)
loc$total_counts <- colSums(count)
X <- loc[, c("x", "y")]
norm_expr <- stabilize(count)
resid_expr <- regress_out(norm_expr, sample_info = loc)

re_SpatialDE <- spatialDE::run(resid_expr, coordinates = X)
    

# calculate AUC
#######################

library(pROC)

# AUC of BOOST-MI
auc(roc(gamma, re_BOOST_MI$BF_neg))

# AUC of BOOST-GP: need to tranfrom PPI to Bayes factor
auc(roc(gamma, 19*(re_BOOST_GP$PPI/(1-re_BOOST_GP$PPI + 10^(-8)))))   

# AUC of SPARK
auc(roc(gamma, re_SPARK$combined_pvalue))

# AUC of SpatialDE: need to reorder the result
re_SpatialDE_pvalue <- numeric(100)
for (i in 1:100){
  re_SpatialDE_pvalue[i] <- re_SpatialDE$pval[re_SpatialDE$g == paste0('X', i)]
}
auc(roc(gamma, re_SpatialDE_pvalue))
    
# AUC of BinSpect: need to reorder the result
re_binSpect_rank_pvalue <- numeric(100)
re_binSpect_kmeans_pvalue <- numeric(100)
    for (i in 1:100){
      re_binSpect_rank_pvalue[i] <- re_binSpect_rank$p.value[re_binSpect_rank$genes == paste0('X', i)]
      re_binSpect_kmeans_pvalue[i] <- re_binSpect_kmeans$p.value[re_binSpect_kmeans$genes == paste0('X', i)]
    }
auc(roc(gamma, re_binSpect_rank_pvalue))
auc(roc(gamma, re_binSpect_kmeans_pvalue))


