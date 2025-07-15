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



# ---- imputation-test --------------------------------------------------------------
#true model
fit_test <- lmer(y ~ x1+x2+x3+x4+x5-1  + (x1+x2+x3-1 | group_id), data = ds_test_use)
summary(fit_test)


skip_var_index <- c("x1", "within_sub_id")
mi_meth <- "2l.pan"
save_path <- "~/data/imputed_data/"
prob <- 0.1
coeff_size <- "big"




vs_imputation2024a <- function(incomplete_data, skip_var_index, mi_meth, prob,
                               imputed_data_save_path, coeff_size){

  replicates <- incomplete_data %>% dplyr::select(replicat) %>% unique() %>% dplyr::pull()
  gs         <- length(unique(incomplete_data$group_id)) #assuming the incomplete and complete data dimensions are the same
  wgs        <- length(unique(incomplete_data$within_sub_id))
  mr         <- paste0(prob*100, "%")

  file_name1 <- paste0(paste("c", gs, wgs, mr, coeff_size, sep ="_"), "/", paste("r", replicates, sep = "_"))
  file_name2 <- paste0(paste("c", gs, wgs, mr, coeff_size, sep ="_"), "/", paste("MAIE_MSIE", replicates, sep = "_"))



  #assuming for the same condition
  complete_data <- dt_list[replicates] %>% do.call(rbind, .)

  ds_com_use <- complete_data %>%
    dplyr::select(com_y = y, com_x2=x2, com_x4=x4, com_x6=x6,
                  within_sub_id, group_id)


  pred_matrix <- make.predictorMatrix(incomplete_data)
  pred_matrix[skip_var_index,] <- 0
  pred_matrix[,skip_var_index] <- 0
  pred_matrix["y","x2"] <- 2
  pred_matrix["y","x3"] <- 2

  pred_matrix["x2","x3"] <- 2
  pred_matrix["x3","x2"] <- 2

  pred_matrix["x2", "y"] <- 2
  pred_matrix["x3", "y"] <- 2
  pred_matrix[,"group_id"]<- -2

  ds_mids <- mice::mice(incomplete_data, predictorMatrix = pred_matrix, method = mi_meth, m=10, seed=12345 , maxit=5)
  ds_mi_long <- complete(ds_mids,action = "long",include = TRUE) %>%
    dplyr::group_by(.imp) %>%
    dplyr::mutate(
      z1 = x1,
      z2 = x2,
      z3 = x3,
      z4 = x4,
      z5 = x5,
      z6 = x6,
      z7 = x7,
      z8 = x8,
      z9 = x9,
      z10= x10,
      # `.imp` = as.numeric(`.imp`),
    ) %>%
    dplyr::ungroup()


  ds_mi_long1 <- ds_mi_long %>%
    dplyr::left_join(ds_com_use, by=c("group_id", "within_sub_id")) %>%
    dplyr::filter(`.imp` != 0  ) %>%
    dplyr::mutate(
      ve_y     = (y-com_y)^2,
      ve_x2    = (x2-com_x2)^2,
      ve_x4    = (x4-com_x4)^2,
      ve_x6    = (x6-com_x6)^2,

      abs_ve_y = abs(y-com_y),
      abs_ve_x2= abs(x2-com_x2),
      abs_ve_x4= abs(x4-com_x4),
      abs_ve_x6= abs(x6-com_x6),
    ) %>%
    dplyr::select(
      `.imp`, `.id`, within_sub_id, group_id,
      ve_y, ve_x2,ve_x4, ve_x6,
      abs_ve_y, abs_ve_x2, abs_ve_x4, abs_ve_x6,
    )

  write_rds(ds_mi_long, paste0(imputed_data_save_path, file_name1, ".rds"))
  write_rds(ds_mi_long1, paste0(imputed_data_save_path, file_name2, ".rds"))


  ds_mi_mids <- mice::as.mids(ds_mi_long)
  ds_mi_list <- miceadds::mids2datlist(ds_mi_mids)

  model_fit <- with(ds_mi_list, exp=lme4::lmer(y ~ x1+x2+x3+x4+x5 -1 + (z1+z2+z3-1| group_id)))

  model_res  <- miceadds::lmer_pool(model_fit)
  res <- summary(model_res)
  mod_table <- cbind(rownames(res), res) %>% dplyr::rename("Variables" = `rownames(res)` )
  mod_result <- mod_table %>%
    dplyr::mutate(
      Variables = dplyr::recode(Variables, "group_id.z1"="z1",
                                "group_id.z2"="z2",
                                "group_id.z3"="z3" ),

      results = sprintf("%0.4f", est),
      se      = sprintf("%0.3f", se),
      t       = sprintf("%0.3f", t),
      p       = sprintf("%0.3f", p),
      lwr     = sprintf("%0.3f", `lo 95`),
      upr     = sprintf("%0.3f", `hi 95`),
    ) %>%
    dplyr::filter(
      Variables %in% fixed_effects |
        Variables %in% rand_effects )


}#end of imputation fun

test_output <- vs_imputation2024a(ds_miss_generate1, skip_var_index, mi_meth, prob = prob,
                                  imputed_data_save_path=save_path, coeff_size=coeff_size)



