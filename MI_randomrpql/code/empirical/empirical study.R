#-------empirical-------------------------------------

library(rpql)
library(mice)
library(micemd)
library(lme4)
library(dplyr)
ds_combine_final_use <- readRDS("~/uk_national_child_development_study/final_cleaned/ds_combine_final_use1.rds")

ds=ds_combined_use

ds=ds%>%
  filter(all_no_na=="TRUE")%>%
  filter(wave!="G")%>%
  mutate(wave=case_when(
    wave=="B"~1,
    wave=="C"~2,
    wave=="D"~3,
    wave=="E"~4,
    wave=="F"~5
  ))%>%
  group_by(MCSID)%>%
  mutate(n=dplyr::n())%>%
  filter(n>4)%>%
  dplyr::select(-n)%>%
  # recode data
  mutate(mwork=case_when(
    DCWRK00==1|DCWRK00==2|DCWRK00==5|DCWRK00==7~1,
    DCWRK00==3|DCWRK00==4|DCWRK00==6|DCWRK00==10~0
  ),
  fwork=case_when(
    DCWRK00==1|DCWRK00==3|DCWRK00==7~1,
    DCWRK00==2|DCWRK00==4|DCWRK00==8~0

  ),
  hh=case_when(
    DROOW00==1|DROOW00==2~1,
    DROOW00==3|DROOW00==4|DROOW00==5|DROOW00==6~2,
    DROOW00==7|DROOW00==8|DROOW00==9~3
  ),
  DMINH00=case_when(
    DMINH00==1~1,
    DMINH00==2~0,
    DMINH00==3~NA_real_

  ),
  DFINH00=case_when(
    DFINH00==1~1,
    DFINH00==2~0,
    DFINH00==3~NA_real_

  ))

m=lmer(sdq_total~1+(1|MCSID),data=ds)
jtools::summ(m)




skip_vars=names(ds)[c(2:3,5:9,11:12,15,18,19,25)]
ds1=ds%>%
  dplyr::select(-skip_vars)

data=ds1


x1=rep(1,nrow(data))
z1=rep(1,nrow(data))

data=cbind(x1,z1,data)
num=ncol(data)-4

random_effect <-paste0("z",seq(1,num,1))
fix_effect <- paste0("x",seq(1,num,1))

random_effect_name <-cbind(old=names(data)[-c(1,3,4,5)],new=random_effect)

fix_effect_name <-cbind(old=names(data)[-c(2:5)],new=fix_effect)

data=cbind(data[,-c(2:5)],data[,-c(1,3,4,5)],data[,c(3,5,4)])
names(data)=c(fix_effect,random_effect,"participant_id","y","wave")

#------training data testing data----------------------------------
set.seed(1234)

ind=sample(unique(data$participant_id),300)

data_f=
  data%>%
  filter(participant_id%in% ind)

#--------stage 1--------------------------------

mi_meth        <- "2l.pan"
incomplete_data=data_f[,c(1,27,2:13,28:29)]
incomplete_data$participant_id=as.numeric(as.factor(incomplete_data$participant_id))
skip_var_index=c("participant_id")
gs         <- length(unique(data_f$participant_id)) #assuming the incomplete and complete data dimensions are the same
wgs        <- length(unique(data_f$wave))

pred_matrix <- make.predictorMatrix(incomplete_data)
pred_matrix[, c("y","x1")] <- 2 #The main difference is here!! only random intercept!!!!
pred_matrix[skip_var_index,] <- 0
pred_matrix[,skip_var_index] <- 0

pred_matrix[, "participant_id"]<- -2
pred_matrix["participant_id", ]<- -2
#c(1,1,1,1,1,1,1,1,1,1,1,0,-2)
diag(pred_matrix) <-0

ds_mids <- mice::mice(incomplete_data, predictorMatrix = pred_matrix,
                      method = mi_meth, m=10, seed=1234567 , maxit=5)





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
    z11= x11,
    z12= x12,
    z13= x13
    # `.imp` = as.numeric(`.imp`),
  ) %>%
  dplyr::ungroup()

set.seed(1234567)


ds_mi_long1 <- ds_mi_long %>% #we directly use the generated missing data
  dplyr::filter(`.imp` != 0  )%>%
  filter(participant_id%in% ind)#excluding the missing data part







#-------stage2--------------------
base::source(file="~/code/model_selection_functions.R")
data=ds_mi_long1
data=data%>%
  mutate(x14=wave,
         z14=wave)


id_vector<- c("participant_id","wave")
no_scale <- c("participant_id", "z1","x1","wave","x6","x7","x8","x9","x10","x11","x12","x14",
              "z6","z7","z8","z9","z10","z11","z12","z14")
random_effect <-c(paste0("z",seq(1,14,1)))
fix_effect <- c(paste0("x",seq(1,14,1)))

# random_effect=random_effect[-6]
# fix_effect =fix_effect [-6]
fixed_effect_name=data.frame(new=paste0("x",seq(1,14,1)))%>%mutate(old=fix_effect)
random_effect_name=data.frame(new=paste0("z",seq(1,14,1)))%>%mutate(old=random_effect)


ycol <- "y"
cluster_id="participant_id"
im=1
verbose<- TRUE
family<-"gaussian"



id_vector<- c("participant_id", "wave")
no_scale <- c("participant_id", "z1","x1","wave","x6","x7","x8","x9","x10","x11","x12","x14",
              "z6","z7","z8","z9","z10","z11","z12","z14")
ycol <- "y"
cluster_id="participant_id"
im=1
verbose<- TRUE
family<-"gaussian"

result1<- rpql_selection1(q1_fix = 7,
                          q1_random=4,
                          bt=10,
                          im=1,
                          data=data,
                          fixed_effect=fix_effect,
                          random_effect=random_effect,
                          ycol,
                          id_vector,
                          no_scale,
                          cluster_id="participant_id",
                          time="wave",
                          family = "gaussian",
                          ci_criteria=5,
                          pen.type="adl",
                          lam=exp(seq(from=log(0.5e-3),to=log(1e-14),length.out=100)),
                          verbose = TRUE,
                          random_draw_sample=0.7

)
rlist::list.save(result1,paste0("~/empirical analysis/","result1.rds"))

#------------stage3--------------------------------------------------
lam=exp(seq(from=log(1e-2),to=log(1e-5),length.out=100))
result2<- rpql_selection2(q2_fix =4,q2_random=3,bt=10, #q2 will be a mix issue in this function
                          data=data,
                          fixed_effect=fix_effect,
                          random_effect=random_effect,
                          ycol=ycol,
                          id_vector=id_vector,
                          no_scale=no_scale,
                          cluster_id="participant_id",
                          time="wave",
                          family = family,
                          pen.type="adl",
                          result=result1,
                          lam=lam,
                          verbose = TRUE,
                          random_draw_sample = 1,
                          ci_criteria=5,
                          use_fre_fun1 =T,
                          use_lmer = F)


rlist::list.save(result2,paste0("~/empirical analysis/","result2.rds"))
#-------------stage4--------------------------------------------------------
base::source(file="~/code/rwe.R")
base::source(file="~/code/hidden_functions.R")

tmp_co_t=NULL



rank=change_rank(res2=result2,fix_effect,random_effect)

rp_change_fun(rank,plot_fixed_effect=T,plot_random_effect=T,fix_effect,random_effect,fixed_effect_name,random_effect_name)





#--------testing-------------------------------------------------------------------------------------
ds_testing=
  ds%>%
  filter(!participant_id%in%ind)

# imputation again
mi_meth        <- "2l.pan"
incomplete_data=data_testing[,c(1,27,2:13,28:29)]
incomplete_data$participant_id=as.numeric(as.factor(incomplete_data$participant_id))
skip_var_index=c("participant_id")
gs         <- length(unique(data_f$participant_id)) #assuming the incomplete and complete data dimensions are the same
wgs        <- length(unique(data_f$wave))

pred_matrix <- make.predictorMatrix(incomplete_data)
pred_matrix[, c("y","x1")] <- 2 #The main difference is here!! only random intercept!!!!
pred_matrix[skip_var_index,] <- 0
pred_matrix[,skip_var_index] <- 0

pred_matrix[, "participant_id"]<- -2
pred_matrix["participant_id", ]<- -2
#c(1,1,1,1,1,1,1,1,1,1,1,0,-2)
diag(pred_matrix) <-0

ds_mids <- mice::mice(incomplete_data, predictorMatrix = pred_matrix,
                      method = mi_meth, m=10, seed=1234567 , maxit=5)





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
    z11= x11,
    z12= x12,
    z13= x13
    # `.imp` = as.numeric(`.imp`),
  ) %>%
  dplyr::ungroup()

set.seed(1234567)




ds_mi_mids <- mice::as.mids(ds_mi_long)
ds_mi_list <- miceadds::mids2datlist(ds_mi_mids)

model_fit <- with(ds_mi_list, exp=lme4::lmer(y ~ x1+wave-1 + (1| participant_id))) #true model

model_res  <- miceadds::lmer_pool(model_fit)
res <- summary(model_res)
mod_table <- cbind(rownames(res), res) %>% dplyr::rename("Variables" = `rownames(res)` )
mod_result <- mod_table %>%
  dplyr::mutate(
    Variables = dplyr::recode(Variables, "participant_id.z1"="z1", ),

    results = sprintf("%0.4f", est),
    se      = sprintf("%0.3f", se),
    t       = sprintf("%0.3f", t),
    p       = sprintf("%0.3f", p),
    lwr     = sprintf("%0.3f", `lo 95`),
    upr     = sprintf("%0.3f", `hi 95`),

  )
