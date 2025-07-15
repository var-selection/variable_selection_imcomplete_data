# generate missing

library(mice)
library(micemd)
library(miceadds)
library(dplyr)
library(missMethods) #generating missing values for simulation
library(ggplot2)
library(multilevel) #calculating icc
library(lme4)
library(lmerTest)
library(broom.mixed)
#library(mitml)
library(readr)


fixed_effects <- c(paste0("x", 1:10))
rand_effects  <- c(paste0("z", 1:10))


#within-group size=30, group size=15
#p=q=10
ds_test <- readRDS("~/data/data_15_30_small/cs_data_15_30_small_2.rds")

#generating missing rate=10%
#missingness on x2, x4, x6
ds_test_use <- ds_test %>%
  dplyr::select(y, all_of(fixed_effects), within_sub_id, group_id)


complete_data <- ds_test_use


#cor(complete_data)

# ---- generate_missing_data_function ------------------------------------------

#1. Generate missing data 2024 version --> much easy this time
#generate_miss <- function(data_long, cmis, contr, varnum, pnum, id, time, predictor){#func start
#https://cran.r-project.org/web/packages/missMethods/vignettes/Generating-missing-values.html

generate_miss2024 <- function(data, prob){
  ds_miss_generate <- data %>%
    #delete_MAR_censoring(., p=prob, cols_mis  = "y",  cols_ctrl ="x2", where="both", sorting=F ) %>% #generating y
    delete_MAR_censoring(., p=prob, cols_mis  = "x2",  cols_ctrl ="x3", where="both", sorting=F ) %>%
    delete_MAR_censoring(., p=prob, cols_mis  = "x4",  cols_ctrl ="x5", where="both", sorting=F ) %>%
    delete_MAR_censoring(., p=prob, cols_mis  = "x6",  cols_ctrl ="x7", where="both", sorting=F )

}#end generate_miss2024


incomplete_data <- generate_miss2024(ds_test_use, 0.1)

count_NA(incomplete_data$x2, type = "default")
count_NA(incomplete_data$x4, type = "default")
count_NA(incomplete_data$x6, type = "default")


#plot testing
test_com <- ds_test_use %>% dplyr::select(x2, x3)
md.pattern(incomplete_data[ , c("y", fixed_effects)])
multilevel::ICC1(aov(y ~ as.factor(within_sub_id), data = incomplete_data))







