#--------------model_selection---------------------------
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

base::source(file="~/code/model_selection_functions.R")

beta="small"
# beta="big"
n_group         <- 50 #group size
n_time          <- 5 #within group: cluster-size
#n_time          <- 10 #within group: cluster-size
im=10 # number of imputated dataset
missing_rate="5%"
#missing_rate="20%"
path_ds="simulated_data/"

random_effect <-paste0("z",seq(1,10,1))
fixed_effect <- paste0("x",seq(1,10,1))
id_vector<- c("group_id", "within_sub_id")
no_scale <- c("group_id", "within_sub_id", "z1","x1")
ycol <- "y"
cluster_id="group_id"
verbose<- TRUE
family<-"gaussian"


id=opt$id
#id=204
a=Sys.time()
# idlist=c(7,8)


data_name=paste0(path_ds,"c_",n_group,"_",n_time,"_",missing_rate,"_",beta,"/r_",id,".rds")
data=readr::read_rds(data_name)

#----------stage2---------------------------
# result1_name=paste0(path_res,"result1_",n_group,"_",n_time,"_",missing_rate,"_",beta,"/res1_",id,".rds")
#result1=tryCatch(readr::read_rds(result1_name),error=function(e){return(NA)})
#if(is.na(result1)){
lam<-exp(seq(from=log(1e-1),to=log(1e-10),length.out=100))
result1 <-rpql_selection1(q1_fix = 4,
                          q1_random=4,
                          bt=4,
                          im=10,
                          ds_mids=data,
                          fixed_effect=fixed_effect,
                          random_effect=random_effect,
                          ycol=ycol,
                          id_vector=id_vector,
                          no_scale=no_scale,
                          cluster_id="group_id",
                          time="within_sub_id",
                          family = family,
                          ci_criteria=5,
                          pen.type="adl",
                          lam=lam,
                          verbose = TRUE,
                          random_draw_sample=1)


rlist::list.save(result1,paste0(path_res,"result1_",n_group,"_",n_time,"_",missing_rate,"_",beta,"_",id,".rds"))
}




#---------------stage3-------------------------------------------------------------------

lam=exp(seq(from=log(2),to=log(1e-4),length.out=100))
result2<- rpql_selection2(q2_fix =3,q2_random=3,bt=4, #q2 will be a mix issue in this function
                          data=data,
                          fixed_effect=fixed_effect,
                          random_effect=random_effect,
                          ycol=ycol,
                          id_vector=id_vector,
                          no_scale=no_scale,
                          cluster_id="group_id",
                          time="within_sub_id",
                          family = family,
                          pen.type="adl",
                          result=result1,
                          lam=lam,
                          verbose = TRUE,
                          random_draw_sample = 1,
                          ci_criteria=5,
                          use_fre_fun1 =T,
                          use_lmer = F)

rlist::list.save(result2,paste0(path_res2,"result2_",n_group,"_",n_time,"_",missing_rate,"_",beta,"_",id,".rds"))

#-------------stage 4-----------------------------------------------------------------------------------------------



base::source(file="~/code/rwe.R")
base::source(file="~/code/hidden_functions.R")




fixed_effect_name=data.frame(new=paste0("x",seq(1,10,1)))%>%mutate(old=new)
random_effect_name=data.frame(new=paste0("z",seq(1,10,1)))%>%mutate(old=new)
rank=change_rank(res2=result2,fixed_effect,random_effect)
res_f=rp_change_fun(rank,plot_fixed_effect=F,plot_random_effect=F,fixed_effect,random_effect,fixed_effect_name,random_effect_name)


rlist::list.save(res_f,paste0(path_rank,"selection_results_",n_group,"_",n_time,"_",missing_rate,"_",beta,"_",id,".rds"))
