## Random Forest for multiclass classification (tidymodels inside)

library("vip")
library("ggplot2")
library("tidyverse")
library("tidymodels")
library("data.table")
library("randomForest")
library("themis")
library('ape')
library('phytools')
library('vita')
#library('modelStudio')
#library('DALEX')

#set up the input data and read in the files
{
#setwd("/Users/cotinga/jsb439@cornell.edu/Code/avian_molecular_shifts")

LHT<-as.data.frame(readRDS(file="./RDS/RF_test_data.RDS"))
megaLHT<-readRDS(file="./RDS/megaLHT.RDS")
models<-as.factor(megaLHT$exon_models)
megaLHT<-megaLHT[!colnames(megaLHT) %in% colnames(LHT)]
megaLHT<- as.data.frame(lapply(megaLHT[,c(2,3, 4, 7, 8, 9, 10)], as.factor))

#getting phylo distance pcoas
phylo.dist.gene <-readRDS('./RDS/pcoa.dist.gene.RDS')
phylo.dist.gene <- phylo.dist.gene$vectors[,c(1)]

phylo.dist <- readRDS('./RDS/pcoa.dist.RDS')
phylo.dist <- phylo.dist$vectors[,c(1)]

# #generate time distance predictors using MEMs (from spatialRF)
tree<-readRDS(file="./RDS/RF_timetree.RDS")
distance_matrix<-cophenetic.phylo(as.phylo(tree))
# mems <- spatialRF::mem(distance.matrix = distance_matrix)
# 
# mem.rank <- spatialRF::rank_spatial_predictors(
#   distance.matrix = distance_matrix,
#   spatial.predictors.df = mems,
#   ranking.method = "moran"
# )

#mems <- mems[, mem.rank$ranking]
# #plot(mems$mem_2~ mems$mem_3)

merged<-cbind(exon_models=models, LHT, megaLHT, phylo.dist)

merged <- merged %>%
  select(c(exon_models, mass, mean.clutch, gen_length,
           breeding, longevity, chickPC1, latitude, freshwater,
           marine, migrant, diet, nocturnal, forest#, 
           # mem_1,
           # mem_2,
           # mem_3,
           # mem_4,
           # mem_5,
           # mem_6,
           # mem_7,
           # mem_8,
           # mem_9,
           # mem_10
           )) #phylo.dist

}

#how many of each class do we have in the given partition set?
merged %>% group_by(exon_models) %>%
  summarise(N=n())

#### Data splitting

#split into training and test sets
{
set.seed(5)
merged_dt <- merged#select(merged, -c(`exon_models`))
merged_split <- initial_split(merged_dt, strata = exon_models, prop = 0.7)
merged_train <- training(merged_split)
merged_test <- testing(merged_split)
}

#### Preprocessing
#We use tidymodels to build a recipe for data preprocessing:
{
merged_recipe <- merged_train %>%
  recipe(exon_models ~ .) %>%
  step_corr(all_numeric(), threshold = 0.95) %>%
  step_normalize(all_numeric(), -all_outcomes()) %>%
  step_dummy(all_nominal_predictors(), one_hot=T) %>%
  #step_upsample(exon_models, over_ratio=0.5) #%>%
  #step_rose(exon_models, over_ratio=0.5, seed=5)
  step_adasyn(exon_models, over_ratio=0.5, neighbors=2, seed=1)
#class_weights()

#step_interact( ~ all_predictors())
#%>%
#step_zv(all_numeric(), -all_outcomes()) %>%
#step_normalize(all_numeric(), -all_outcomes()) %>%
#step_impute_knn(all_numeric(), neighbors = 5) ## there are no missing data here, but in case!
}

#prep the model
{
prep_merged <- prep(merged_recipe)
print(prep_merged)
}

#juice the model
{
training_set <- juice(prep_merged)
head(training_set)
}

#### Model building
# We now specify the structure of our model:
#   
# - hyperparameters to tune: `mtry` (number of features to sample for each tree) and `min_n` (minimum number of data points in a node to allow further splitting)
# - number of trees in the forest
# - the problem at hand (classification)
# - the engine (R package)
# 
# Then we put this in a workflow together with the preprocessing recipe
{
tune_spec <- rand_forest(
  mtry = tune(),
  trees = 1000,
  min_n = tune(),
  #proxmity=distance_matrix,
) %>%
  set_mode("classification") %>%
  set_engine("randomForest", importance=TRUE) #tree # proxmity=distance_matrix


tune_wf <- workflow() %>%
  add_formula(exon_models ~ .) %>%
  add_model(tune_spec)
}

#### Tuning the hyperparameters
#We use k-fold cross-validation to tune the hyperparameters in the training set
#10 folds repeated 10 times
{
trees_folds <- vfold_cv(training_set, v = 10, repeats = 10, strata=exon_models)
print(trees_folds)
}

#tune the parameters
{
doParallel::registerDoParallel(cores = 10)

tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 20, control=control_grid(save_pred=T) ## n. of tuning combinations
)

print(tune_res)
}

#plot min_n vs mtry
{
library("repr")
options(repr.plot.width=14, repr.plot.height=8)

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

}

#set up the grid for parameter tuning
{
m <- round(sqrt(ncol(training_set)-1),0)
print(m)
rf_grid <- grid_regular(
  mtry(range = c(m-3, m+6)),
  min_n(range = c(1, 15)),
  levels = 4
)
print(rf_grid)

}

#tune the parameters
{
regular_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = rf_grid, control=control_grid(save_pred=T)
)

print(regular_res)
}

#look at model performance across different model parameters
{
regular_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")
}

#### Final model

#We now select the best model from the hyperparameters tuning, and fit it to the training set:
#selecting the best model based on AUC:
{
best_auc <- select_best(tune_res, "roc_auc")
print(best_auc)
}

#finalise the model:
{
final_rf <- finalize_model(
  tune_spec,
  best_auc
)
  
print(final_rf)

}

#finalise the workflow and fit it to the initial split (training and test data):
{
final_wf <- workflow() %>%
  add_recipe(merged_recipe) %>%
  add_model(final_rf)

final_res <- final_wf %>%
  last_fit(merged_split)
}

#evaluate the fine-tuned RF model:
{
print(final_res)
final_res %>%
  collect_metrics()
}

#get variable importance:
{
imp<-final_res %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_engine() %>% 
  #vip(num_features = 20, geom = "point")
  #varimp(conditional=T)
  #vip(num_features = 25, geom="point", all_permutations=F, include_type=T, method='model')+theme_bw() 
  # 
  #PIMP(X=training_set[,-c(8)], y= training_set$exon_models, rForest=., parallel=T, ncores=10) %>%
  #PimpTest() %>% 
  #summary()
  
  vip(method='permute', target="exon_models", metric="error", pred_wrapper = predict, train=training_set, nsim=500, parallel=F,
      geom = "violin", include_type=T, all_permutations=F, num_features=25, mapping = aes_string(fill = "Variable"),
      aesthetics=list(scale="width", draw_quantiles=0.5, adjust=2.5, color = "black", size = 0.1)) +
  #ggplot2::scale_fill_manual(values = fill.color) + 
  ggplot2::ylab("") + ggplot2::ylab("Increase in error when permuted") + 
  ggplot2::theme_minimal() + theme(panel.grid.major.x = element_blank(),
                                   panel.grid.minor.x = element_blank()) +
  ggplot2::theme(legend.position = "none", axis.text.y = element_text(size = 10))

}

pdf(file="RF_importance.pdf", height=6, width=3)
imp+theme(axis.text.y = element_text(size = 7), axis.text.x= element_text(size = 8),
          axis.title.x = element_text(size = 8)) + scale_fill_manual(values=rep("grey", 24))
dev.off()

#check randomForestExplainer -- tree depth 
{
require(randomForestExplainer)
final_res %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_engine() %>% 
  plot_min_depth_distribution()
}


