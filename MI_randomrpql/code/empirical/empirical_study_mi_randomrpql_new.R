library(readr)
library(dplyr)
library(mice)
library(miceadds)
library(micemd)
#---chanmice#---change-column-names----------------------------------
child_data_dictionary_structured <- read_csv("~/Model selection2/empirical data/child_data_dictionary_structured.csv")

new_data=read_rds("~/Model selection2/empirical data/new_data.rds")

df=new_data
nm_data=
names(df)%>%
  as.data.frame()
names(nm_data)="col_names"

nm_data1=
nm_data%>%
  left_join(child_data_dictionary_structured%>%
              dplyr::select(-short_description,-value_summary),by=c("col_names"="column_name"))%>%
  mutate(col_name_new=paste0(series_index,"_",year_from_label))

names(df)=nm_data1$col_name_new

# homework from 2004 need to divded by 100

df1=
df%>%
tidyr::pivot_longer(
  cols = matches("_(19\\d{2}|200\\d|201[0-8])$"),
  names_to = c(".value", "year"),
  names_pattern = "^(.*)_(\\d{4})$",
  names_transform = list(year = as.integer)
  )%>%
  dplyr::select(child_id=child_id_2025,mother_id=Mother_id_2025,year,MATHZ,race=Race_2025,
                birth_of_child_of_mom=birht_of_child_of_mom_2025,
                age_of_child, gender,DOB,age_of_mom,education,child_number,
                weight,COGNZ,EMOTZ,BPI,height,SPPCS,SPPCG,DIGITZ,
                RECOGZ,COMPZ,PPVTZ,homework_in_school,homework_after_school,
                income
                )%>%
  ungroup()%>%
  mutate(homework_in_school=if_else(year>=2004,round(homework_in_school/100),homework_in_school),
         homework_after_school=if_else(year>=2004,round(homework_after_school/100),homework_after_school),
         # bmi=703 * weight / (height^2),
          homework_hours=(homework_in_school+homework_after_school),
          age_of_child=round(age_of_child/12,2))%>%
     dplyr::select(-homework_in_school,-homework_after_school)%>%
group_by(child_id)%>%
  mutate(DOB=max(DOB,na.rm=T), gender=max(gender,na.rm=T),
         age_of_child=if_else(is.na(age_of_child),year-DOB,age_of_child),
         age_of_child=if_else(age_of_child<0,NA,age_of_child),
         age_of_child_new=(year-DOB))%>%
  arrange(child_id,year)%>%
  ungroup()
#------------choose-which-year-need-to-use-------------------------
# using the age of child selected the outcome
df2=
  df1%>%
  filter(age_of_child_new>=5 & age_of_child_new <=13)%>%
  group_by(child_id)%>%
  mutate(miss=any(is.na(MATHZ)))%>%
  filter(!miss)%>%
  mutate(n=n())%>%
  filter(n>4)%>%
  mutate(wave=case_when(age_of_child_new==5~1,
                        age_of_child_new==7~2,
                        age_of_child_new==9~3,
                        age_of_child_new==11~4,
                        age_of_child_new==13~5
  ))#n_distinct(df2$child_id) #1181


child_id_keep=
df2%>%
  group_by(mother_id)%>%
  filter(year==min(year))%>%
  summarize(keep_id=unique(child_id))%>%
  mutate(n=n())%>% # 970 month has the same age child 26 month has two same age child random select one
  filter(n==1)%>%
  ungroup()

df3=
df2%>%
  filter(child_id%in%child_id_keep$keep_id)%>%
  group_by(mother_id)%>%
  mutate(n_m=n_distinct(child_id))%>% #n_distinct(df3$child_id) 970
  dplyr::select(child_id,mother_id,MATHZ,race,birth_of_child_of_mom,gender,education,
                child_number,weight,height,COGNZ,EMOTZ,BPI,SPPCS,SPPCG,DIGITZ,RECOGZ,COMPZ,
                PPVTZ,homework_hours,wave,age_of_child_new)%>%
  ungroup()


# delete the child with too much missing data
missing_rate=
df3%>%
  group_by(child_id,wave)%>%
  summarize(
    miss=mean(is.na(c_across(SPPCS:homework_hours)))
  )%>%
  group_by(child_id)%>%
  summarize(n_miss=mean(miss))



df4=
df3

m1=lme4::lmer(MATHZ~1+(1|child_id),data=df4)
jtools::summ(m1)
# icc=0.56
names(df4)
cor_pre=round(cor(df4%>%dplyr::select(-child_id,-mother_id,-race,-age_of_child_new),use="pairwise.complete.obs"),2)
diag(cor_pre)=NA
cor_pre=abs(cor_pre)
min(cor_pre,na.rm=T) #0
max(cor_pre,na.rm=T) #0.77

# missing rate
miss_rate_each=
df4%>%
  summarize(across(MATHZ:homework_hours,~mean(is.na(.x))))

a=
  df4%>%
  group_by(child_id,wave)%>%
  summarize(
    miss=mean(is.na(c_across(MATHZ:homework_hours)))
  )%>%
  group_by(child_id)%>%
  summarize(n_miss=mean(miss))
#----------------imputation----------------------------------------
names(df4)

mi_meth        <- "2l.pan"
incomplete_data=df4%>%dplyr::select(-mother_id)%>%mutate(x1=1)

skip_var_index=c("child_id")
gs         <- length(unique(df4$child_id)) #assuming the incomplete and complete data dimensions are the same
wgs        <- length(unique(df4$wave))

pred_matrix <- make.predictorMatrix(incomplete_data)
pred_matrix[, c("MATHZ","x1")] <- 2 #The main difference is here!! only random intercept!!!!
pred_matrix[skip_var_index,] <- 0
pred_matrix[,skip_var_index] <- 0

pred_matrix[, "child_id"]<- -2
pred_matrix["child_id", ]<- -2

diag(pred_matrix) <-0

ds_mids <- mice::mice(incomplete_data, predictorMatrix = pred_matrix,
                      method = mi_meth, m=10, seed=1234567 , maxit=20)



plot(ds_mids)
ds_mi_long <- complete(ds_mids,action = "long",include = TRUE)

#-------split---400(train):150(testing)-----------------------
set.seed(1234567)
n_distinct(ds_mi_long$child_id)
# 488 person
# random select 300
ind=sample(unique(ds_mi_long$child_id),400)
ind_test=sample(unique(ds_mi_long$child_id),150)



ds_mi_long1 <- ds_mi_long %>%
  dplyr::filter(`.imp` != 0  )%>%
  filter(child_id%in% ind)

ds_test <- df4 %>%
  filter(child_id%in% ind_test)



# readr::write_rds(ds_mi_long1, paste0("~/model-selection/empirical_analysis/train_data_round2", ".rds"))
# readr::write_rds(ds_test, paste0("~/model-selection/empirical_analysis/test_data_round2", ".rds"))

#-------------------run-selection-------------------------------

#--------load-data------------------------------------------------------
training_data=ds_mi_long1
testing_data=ds_test

#-------change name of fix and random effect--------------------
base::source(file="~/empirical/code/model_selection_model_functions.R")

id_vector<- c("child_id","wave")
no_scale <- c("child_id","x1","x2","x4","z1","z2","z4","wave")

training_data$gender=training_data$gender-1

random_effect_name=names(training_data)[c(5:22,24)]
fixed_effect_name=names(training_data)[c(5:22,24)]

train_ds=cbind(training_data[,c(".imp","child_id","wave","MATHZ","x1",fixed_effect_name[-19])],
               training_data[,c("x1",random_effect_name[-19])])

names(train_ds)=c(".imp","child_id","wave","MATHZ",paste0("x",seq(1,length(random_effect_name),1)),
                  paste0("z",seq(1,length(random_effect_name),1)))

names_index=cbind(names(train_ds)[2:23],
                  c("child_id","wave","MATHZ","x1",fixed_effect_name[-19]))
fixed_effect=paste0("x",seq(1,length(random_effect_name),1))
random_effect=paste0("z",seq(1,length(random_effect_name),1))

ycol <- "MATHZ"
cluster_id="child_id"
im=10
verbose<- TRUE
family<-"gaussian"


fixed_effect=fixed_effect[-2]
random_effect=random_effect[-c(2,3,4,5)]
fixed_effect_name=fixed_effect_name[-1]
random_effect_name=random_effect_name[-c(1,2,3,4)]

time_res1=system.time({
  result1<- rpql_selection1(q1_fix = 6,
                            q1_random=4,
                            bt=400,
                            im=10,
                            data=train_ds,
                            fixed_effect=fixed_effect,
                            random_effect=random_effect,
                            ycol=ycol,
                            id_vector=id_vector,
                            no_scale=no_scale,
                            cluster_id="child_id",
                            time="wave",
                            family = "gaussian",
                            ci_criteria=5,
                            pen.type="adl",
                            lam=exp(seq(from=log(0.5e-2),to=log(1e-10),length.out=100)),
                            verbose = TRUE,
                            random_draw_sample=1

  )
})



# }
time_res2=system.time({
  lam=exp(seq(from=log(0.5e-1),to=log(1e-4),length.out=100))
  result2<- rpql_selection2(q2_fix =4,q2_random=4,bt=400, #q2 will be a mix issue in this function
                            data=train_ds,
                            im=10,
                            fixed_effect=fixed_effect,
                            random_effect=random_effect,
                            ycol=ycol,
                            id_vector=id_vector,
                            no_scale=no_scale,
                            cluster_id="child_id",
                            time="wave",
                            family = family,
                            pen.type="adl",
                            result=result1,
                            lam=lam,
                            verbose = TRUE,
                            random_draw_sample = 1,
                            ci_criteria=5,
                            use_fre_fun1 =F,
                            use_lmer = F)
})

#----------------rank--------------------------------
#-----------ranking---------------------------------------------------------------------
library(foreach)
library(doParallel)
library(rlist)
library(parallel)
library(tidyr)
library(lares)
library(dplyr)
library(MASS)
library(missMethods) #generating missing values for simulation
library(miceadds)
library(mice)
library(rpql)
library(PlackettLuce)

base::source(file="~/MI_randomrpql/code/rwe.R")
base::source(file="~/MI_randomrpql/code/hidden_functions.R")

tmp_co_t=NULL

fixed_effect_list=names_index[-c(1,2,3,5),]%>%as.data.frame()
names(fixed_effect_list)=c("new","old")

random_effect_list=names_index[-c(1,2,3,5,6,7,8),]%>%as.data.frame()
random_effect_list[,1]=random_effect
names(random_effect_list)=c("new","old")

time_res3=system.time({
  rank=change_rank(res2=result2,fixed_effect,random_effect)

  res_f=rp_change_fun(rank,plot_fixed_effect=T,plot_random_effect=T,fixed_effect,random_effect,fixed_effect_list,random_effect_list)
})



# t1=time_res1[3]/60
t1=time_res1[3]/60
t2=time_res2[3]/60
t3=time_res3[3]/60
t4=t1+t2+t3
# time in minutes
t_total=cbind(time_res1=t1,time_res2=t2,rank=t3,total_time=t4)

#--------run-glmm-------------------------------------------------------------------------------------

#----imputate---------------
mi_meth        <- "2l.pan"
incomplete_data=testing_data%>%dplyr::select(-mother_id)%>%mutate(x1=1)

skip_var_index=c("child_id")
gs         <- length(unique(testing_data$child_id)) #assuming the incomplete and complete data dimensions are the same
wgs        <- length(unique(testing_data$wave))

pred_matrix <- make.predictorMatrix(incomplete_data)
pred_matrix[, c("MATHZ","x1","x13")] <- 2 #The main difference is here!! only random intercept!!!!
pred_matrix[skip_var_index,] <- 0
pred_matrix[,skip_var_index] <- 0

pred_matrix[, "child_id"]<- -2
pred_matrix["child_id", ]<- -2

diag(pred_matrix) <-0

ds_mids <- mice::mice(incomplete_data, predictorMatrix = pred_matrix,
                      method = mi_meth, m=10, seed=1234567 , maxit=20)



plot(ds_mids)

ds_mi_list <- miceadds::mids2datlist(ds_mids)
# #
#
#
#
model_fit <- with(ds_mi_list, exp=lme4::lmer(MATHZ ~ PPVTZ+RECOGZ+ (1+RECOGZ|child_id))) #true model
#
#
#
model_res  <- miceadds::lmer_pool(model_fit)
res <- summary(model_res)
mod_table <- cbind(rownames(res), res) %>% dplyr::rename("Variables" = `rownames(res)` )
mod_result <- mod_table %>%
  dplyr::mutate(
    Variables = dplyr::recode(Variables, "child_id.z1"="z1", ),
#
    results = sprintf("%0.4f", est),
    se      = sprintf("%0.3f", se),
    t       = sprintf("%0.3f", t),
    p       = sprintf("%0.3f", p),
    lwr     = sprintf("%0.3f", `lo 95`),
    upr     = sprintf("%0.3f", `hi 95`),

  )
mod_result






