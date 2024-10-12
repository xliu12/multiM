library(glue)
library(tidyverse)
library(mvtnorm)
library(SuperLearner)
library(origami)
library(fastDummies)

devtools::load_all("MediationClustered")

# load data
load("_example/data_example.RData")

Sname <- "school"

Xnames <- colnames(data_example)[grep("^X", colnames(data_example))]
Wnames <- colnames(data_example)[grep("^W", colnames(data_example))] # W_nj, the cluster size scaled

Aname <- colnames(data_example)[grep("^A_", colnames(data_example))]
Mnames <- colnames(data_example)[grep("^M_", colnames(data_example))]
Yname <- colnames(data_example)[grep("^Y_", colnames(data_example))]

learners_a <- learners_m <- learners_y <- c("SL.gam", "SL.glm","SL.nnet") 

# "cluster mean + cluster dummies" (referred to as cwc.FE)
methd <- data.frame(expand.grid(cluster_opt=c("cwc.FE"), Morder = c("12","21")))
i <- 1


estimates_list <- map(.x = 1:nrow(methd), .f = function(.x){
    i <- .x
    out <- ClusterMed(data = data_example,
                     Sname = Sname,
                     Wnames = Wnames, 
                     Xnames = Xnames,
                     Aname = Aname,
                     Mnames = Mnames, 
                     Yname = Yname, 
                     cluster_opt = methd$cluster_opt[i],
                     Morder = methd$Morder[i],
                     num_folds = 5,
                     learners_a = learners_a,  
                     learners_m = learners_m,
                     learners_y = learners_y 
    )
    
    data.frame(out, methd[.x, ], row.names = NULL)
})

# cluster means ----
methd <- data.frame(expand.grid(cluster_opt=c("cwc"), Morder = c("12","21")))
i <- 1



estimates_cwc_list <- map(.x = 1:nrow(methd), .f = function(.x){
    i <- .x
    out <- ClusterMed(data = data_example,
                     Sname = Sname,
                     Wnames = Wnames, 
                     Xnames = Xnames,
                     Aname = Aname,
                     Mnames = Mnames, 
                     Yname = Yname, 
                     cluster_opt = methd$cluster_opt[i],
                     Morder = methd$Morder[i],
                     num_folds = 5,
                     learners_a = learners_a,  
                     learners_m = learners_m,
                     learners_y = learners_y 
    )
    
    data.frame(out, methd[.x, ], row.names = NULL)
})

results <- c(estimates_list,estimates_cwc_list) %>% reduce(rbind)
results
