library(dplyr)
library(mice)
library(miceadds)
library(rpql)

header.true <- function(data) {
  names(data) <- as.character(unlist(data[1,]))
  data[-1,]
}



#----------------step1-function----------------------------------------
rpql_selection1 <- function(q1_fix = 10, # number of fixed effect used in each bootscrap step
                            q1_random=10,# number of random effect used in each bootscrap step
                            bt=10, # time of bootscrapping
                            im=im,data=data,
                            random_effect, fixed_effect,
                            ycol, id_vector, no_scale,
                            cluster_id="group_id",
                            time="within_sub_id", family = "gaussian",
                            ci_criteria,# select whether ci used to select the optimal model
                            pen.type="adl", lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70)),
                            verbose = TRUE,
                            random_draw_sample = 0.5){
  if(is.null(q1_fix)){
    q1_fix=floor(length(fixed_effect)/2)
  }
  
  if(is.null(q1_random)){
    q1_random=length(t(unique(data%>%dplyr::select({{time}}))))-1
  }
  num_w= length(t(unique(data%>%dplyr::select(time))))
  if(num_w<=q1_random){
    stop("number of random effect used in each bootscrapping step need to be less than the within group size")
  } else{
    
    
    var_list_x <-c(fixed_effect, random_effect, id_vector)
    
    var_list              <- c(var_list_x, ycol)
    var_model_selection_x <- setdiff(var_list_x, id_vector)
    fl                    <- length(fixed_effect)
    rl                    <- length(random_effect)
    rlcov                 <- (rl*(rl+1))/2-rl
    
    p  <-dim(data[,var_list_x])[2]
    v  <-dim(data[,var_model_selection_x])[2]
    n  <-dim(data[,var_list_x])[1]
    
    x_name<-c(fixed_effect, random_effect)
    check_name <- x_name %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = gsub("\\d{1}$","\\" ,name),)
    
    set.seed(1234)
    
    
    
    
    
    #M  <-matrix(0,10*im,v+cate_level)
    # full variable lists
    name_res=as.data.frame(c(fixed_effect,random_effect))
    names(name_res)="name"
    fre_list=name_res # list of number of variables selected
    
    
    rep=floor(q1_fix/q1_random)
    bt_rep=ceiling(bt/rep)
    
    M  <-matrix(0,ncol=bt*im*rep,nrow=v)
    #s=1
    for (s in 1:im){#imputation
      if (verbose) message(sprintf("Starting imputation % 3i of %i", s, im))
      MM           <- data%>%filter(.imp==s)%>%dplyr::select(var_list)
      ### test nonmissing
      #MM           <- ds_mids[,var_list]
      #MM[,"group_id"]=as.numeric(unlist(MM[,"group_id"]))
      #####################################3
      #x            <- as.matrix(MM[,var_list_x])
      y            <- as.matrix(MM[,ycol])
      scale_x_name <- setdiff(var_list_x, no_scale)
      xs <- MM %>%
        dplyr::select(all_of(var_list_x),  )%>%
        dplyr::mutate_at(., vars(scale_x_name), as.numeric) %>%
        dplyr::mutate_at(., vars(scale_x_name), scale)
      
      
      fc=0
      for (j in 1:bt){ #set the bootstrapping as bt.
        if (verbose) message(sprintf("Starting bootstrap % 3i of %i", j, bt))
        unique_id2 <- as.data.frame(unique(data[, cluster_id]))##!!changed
        names(unique_id2)="group_id"
        #IMPORTANT: randomly sample half of the people rather than half of the data.
        random_unique <- as.numeric(sample(unique_id2$group_id, round(dim(unique_id2)[1])*random_draw_sample,  replace = T))
        # non-missing change
        #random_unique        <- as.numeric(unlist(sample(unique_id, replace = T))) %>% data.frame() #random draw based on id not on total sample
        ####
        random_unique1       <-random_unique %>%
          as.data.frame()%>%
          #dplyr::rename(cluster_id = random_unique1) %>%
          dplyr::rename(cluster_id = ".") %>%
          dplyr::arrange(cluster_id) %>%
          dplyr::mutate(
            cluster_id = as.numeric(cluster_id),
            
          ) %>%
          # dplyr::group_by(b_cluster_id) %>% #old way
          dplyr::group_by(cluster_id) %>%
          dplyr::mutate(
            time = as.numeric(rep(1:dplyr::n()))
          ) %>%
          dplyr::summarise(
            weight = max(time) #weight is to summarize how many times of being sampled
          ) %>%
          dplyr::ungroup()
        
        IVs2_x <- as.data.frame(xs) %>%
          dplyr::rename(time = time,
                        cluster_id = cluster_id)
        
        
        smp2 <- IVs2_x %>% #update the id infor
          dplyr::select(time, cluster_id) %>%
          dplyr::right_join(random_unique1, by="cluster_id") %>%  #depends on the random_unique1
          dplyr::group_split(cluster_id) %>%
          purrr::map(function(data){
            rep(data$time, times=max(data$weight)) %>% #replicate by weight
              as.data.frame() %>%
              dplyr::rename(
                time = "."
              ) %>%
              dplyr::mutate(cluster_id = rep(unique(data$cluster_id))) %>%
              dplyr::group_by(time) %>%
              dplyr::mutate(
                seq = rep(1:dplyr::n()),
                b_cluster_id = as.integer(paste0(cluster_id, seq))
              ) %>%
              dplyr::select(
                time, b_cluster_id, cluster_id
              )}) %>% do.call(rbind.data.frame, .)%>%
          group_by(b_cluster_id)%>%
          mutate(cluster_id_rpql=cur_group_id()) # update the cluster_id with a sequential number due to the rpql function
        
        
        
        
        if (family == "gaussian"){
          #Assuming outcome's normal distribution
          outcome2 <- MM %>%
            dplyr::select(cluster_id, time, ycol) %>%
            purrr::map_dfc(., as.numeric) %>% #for continuous outcome, can do this
            dplyr::mutate_at(., vars(ycol), scale) %>% #DV needs to be scale
            dplyr::arrange(cluster_id) %>%
            dplyr::rename(time = time, cluster_id = cluster_id)
        } else{ # if the outcome variable is discrete, DV do not scale
          outcome2 <- data %>%
            dplyr::select(cluster_id, time, ycol) %>%
            dplyr::mutate_at(., vars(c(cluster_id, time)), as.numeric) %>%
            arrange(cluster_id) %>%
            dplyr::rename(time = time, cluster_id = cluster_id) #DV no scale
          outcome2_lmer=outcome2
          
        }
        
        
        IVs2       <- smp2 %>% dplyr::left_join(IVs2_x, by=c("cluster_id", "time"))
        DV2        <- smp2 %>% dplyr::left_join(outcome2, by=c("cluster_id", "time"))
        
        x1   <- IVs2 %>% dplyr::arrange(cluster_id, time) %>% as.matrix()
        y1   <- DV2 %>% dplyr::arrange(cluster_id, time) %>%
          dplyr::select(b_cluster_id, time,cluster_id_rpql, ycol) %>%
          as.matrix()
        
        x2_data <- IVs2 %>%
          dplyr::ungroup() %>%
          dplyr::left_join(DV2, by=c("cluster_id", "time", "b_cluster_id","cluster_id_rpql")) %>%
          dplyr::arrange(cluster_id, b_cluster_id,cluster_id_rpql,time) %>%
          dplyr::select(-time)
        
        x2_data_selective_x_only <- x2_data %>%
          dplyr::select(all_of(var_model_selection_x))%>%
          purrr::map_dfc(., as.numeric) %>% #for no categorical var
          dplyr::select(starts_with("x"))
        
        x2_data_selective_z_only <- x2_data %>%
          dplyr::select(all_of(var_model_selection_x))%>%
          dplyr::mutate_all(., as.numeric)%>%
          dplyr::select(starts_with("z"))
        
        x2_lmer_use <- x2_data %>%
          dplyr::ungroup() %>%
          dplyr::select(all_of(var_model_selection_x), ycol, b_cluster_id,cluster_id_rpql) %>%
          purrr::map_dfc(., as.numeric) #assume all continuous predictors
        
        # fix
        name_fix  <- gtools::mixedsort(grep("x",check_name$name,value=T)) #replace rownames function
        loc_x1    <- which(check_name$name=="x1") # used for simulation
        name_fix  <- name_fix[name_fix!="x1"]
        fix_ind   <- grep("x", check_name$name)
        fix_ind   <- fix_ind[fix_ind!=loc_x1]
        
        # random
        name_random  <- gtools::mixedsort(grep("z",check_name$name,value=T))
        #loc_z1       <- which(check_name$name=="z1") #only for simulation, need CHANGE for package
        #name_random  <- name_random[name_random!="z1"]
        random_ind   <- grep("z",check_name$name)
        #random_ind   <- random_ind[random_ind!=loc_z1]
        
        # proportional of fixed and random which is non-zero
        # number of non-zero parameters
        fix_len <- length(fix_ind)
        ran_len <- length(random_ind)
        total_len <- fix_len+ran_len
        
        
        size_fixed  <- q1_fix #q1 needs to set up a default value for package
        size_random <- q1_random
        
        
        w1          <- gtools::mixedsort(sample(name_fix,size=size_fixed)) #update by name rather than order
        w1=c("x1",w1)
        # the proportion of fix and random with different number of selection times
        for (i in 1:rep){
          
          w2          <- gtools::mixedsort(sample(name_random,size=size_random))
          
          # update
          # number of time vairables selected
          
          fre=cbind(as.data.frame(c(w1,w2)),rep(1,q1_fix+q1_random+1))
          
          names(fre)=c("name","fre")
          
          fre <- name_res %>%
            dplyr::left_join(fre, by="name") %>%
            dplyr::mutate(
              unselect = rep(0),
              fre = dplyr::coalesce(as.numeric(fre), unselect)
            ) %>%
            dplyr::select(-unselect)
          
          
          
          
          
          x2_data2_x          <- cbind(x2_data_selective_x_only[, w1])
          x2_data2_z         <- cbind(x2_data_selective_z_only[, w2])
          
          
          
          
          x2_data3 <- names(cbind(x2_data_selective_x_only,(x2_data_selective_z_only)))%>%
            
            as.data.frame() %>%
            dplyr::rename(name = ".") %>%
            tibble::rownames_to_column(., "number")
          
          
          x2_data4 <- names(cbind(x2_data2_x,x2_data2_z)) %>%
            as.data.frame() %>%
            dplyr::rename(name = ".") %>%
            left_join(x2_data3,by="name")
          
          
          b_fix_effect <- x2_data4 %>% dplyr::filter(name %in% fixed_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
          b_random_effect <- x2_data4 %>% dplyr::filter(name %in% random_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
          
          x <- as.matrix(x2_lmer_use[, b_fix_effect])
          y <- x2_lmer_use[, ycol] %>% dplyr::pull(ycol)
          z <- as.matrix(x2_lmer_use[, b_random_effect])
          
          
          
          list1 <- list(id = list(cluster_id_rpql = x2_lmer_use$cluster_id_rpql), X=x, Z=list(cluster_id_rpql=z), Y=y)
          
          if(pen.type=="adl"){# use lmer to calculate the initial value
            if (family %in% "gaussian"){#for lmer DV continue
              fit1  <- lme4::lmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use)
            } #cate needs to be factor
            if (family %in% "binomial"){#for glmer DV binary
              fit1  <- lme4::glmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use,family = "binomial")} #cate needs to be factor}
            
            xx    <- x %>% as.data.frame()
            fit_sat1 <- build.start.fit(fit1, gamma = 2, cov.groups = NULL) #no categorical predictor --> CHANGE for package
            fit_final1 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                           lambda = lam, pen.type = "adl", hybrid.est = TRUE,
                                           pen.weights = fit_sat1$pen.weights, start = fit_sat1,restarts = 50),
                                   error=function(e){return(NA)})
           } else if(pen.type=="lasso"){ # do not use lmer as initial value
            fit_final1 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                           lambda = lam, pen.type = "lasso", hybrid.est = TRUE,
                                           restarts = 50),
                                   error=function(e){return(NA)})
            
          }
          
          if(anyNA(fit_final1)){fc=fc+1 }
          if(fc>bt*0.05){break} # if 5% bootscraping cannot converage, break the loop and give NA
          if(!anyNA(fit_final1)){
            result      <- summary(fit_final1$best.fit[[ci_criteria]])
            fixcoef     <- as.data.frame(result$fixef) %>%  tibble::rownames_to_column(., "name") %>% dplyr::rename(value = `result$fixef`)
            randcov     <- as.matrix(result$ran.cov$cluster_id_rpql)
            randcoef    <- diag(randcov)%>% as.data.frame() %>% cbind(b_random_effect,.)  %>% dplyr::rename(value = ".", name=b_random_effect)
            togethercoef<- rbind(fixcoef, randcoef)
            
            fre_list=left_join(fre_list,fre,by="name")
            
            name_res=as.data.frame(c(fixed_effect,random_effect))
            names(name_res)="name"
            x2_data5 <- name_res %>%
              dplyr::left_join(togethercoef, by="name") %>%
              dplyr::mutate(
                unselect = rep(0),
                value = dplyr::coalesce(as.numeric(value), unselect)
              ) %>%
              dplyr::select(name, value) %>%
              #dplyr::arrange(name) %>%
              t() %>%
              as.data.frame() %>%
              header.true()
            
            
          } #end of non-null fit_final1
          else{
            fc=fc+1
            x2_data5 <- x2_data3 %>%
              dplyr::mutate(
                value = rep(NA)
              ) %>%
              dplyr::select(name, value) %>%
              dplyr::arrange(name) %>%
              t() %>%
              as.data.frame() %>%
              header.true()
          }#end of null fit_final1
          
          
          
          M[,(s-1)*bt+(j-1)*rep+i] <- as.numeric(t(x2_data5))
          
        }
        
      }#end of bootstrapping
    }#end of imputation
    
    rownames(M)=names(x2_data5)
    fre_res=fre_list%>%
      rowwise()%>%
      mutate(total=sum(c_across(where(is.numeric)),na.rm=T))%>%
      ungroup()%>%
      dplyr::select(name,total)
    
    
    
    round1_res <- list(est=as.data.frame(M), fre=fre_res,family=family)
    
    return(round1_res)
  } # end of  warning num_w<=q1_random
}#end of function1

rpql_selection2 <- function(q2_fix=5,q2_random=3, #the number of randomly selected predictors
                            im=10,
                            bt=10,
                            data=data,
                            fixed_effect,
                            random_effect,
                            ycol,
                            id_vector,
                            no_scale,
                            cluster_id,
                            time,
                            family = "gaussian",
                            #cov.groups=NULL,
                            #categorical_levels = c(0),
                            pen.type="adl",
                            result=res1,
                            lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70)),
                            verbose = TRUE,
                            random_draw_sample = 0.5,
                            ci_criteria,
                            use_lmer=T,# whether the weights of ALasso from lmer
                            use_fre_fun1=T # use the frequency list to calculate the important index from function1
){
  #start function2
  
  if(is.null(q2_fix)){
    q2_fix=floor(length(fixed_effect)/2)
  }
  
  if(is.null(q2_random)){
    q2_random=length(t(unique(data%>%dplyr::select(time))))-1
  }
  num_w= length(t(unique(data%>%dplyr::select(time))))
  if(num_w<=q2_random){
    stop("number of random effect used in each bootscrapping step need to be less than the within group size")
  } else{
    
    result_matrix<- vector(mode = "list", length = im*bt)
    
    var_list_x   <-c(fixed_effect, random_effect, id_vector)
    var_list              <- c(var_list_x, ycol)
    var_model_selection_x <- setdiff(var_list_x, id_vector)
    fl                    <- length(fixed_effect)
    rl                    <- length(random_effect)
    rlcov                 <- (rl*(rl+1))/2-rl
    rep=floor(q2_fix/q2_random) # proportion between fix and random
    
    M    <- result$est%>%
      purrr::map_dfc(., as.character)%>%
      purrr::map_dfc(., as.numeric)#from function1
    m    <-dim(M)[2]
    
    #identify the cate predictors
    p    <-dim(data[,var_list_x])[2]
    n    <-dim(data[,var_list_x])[1]
    
    #x_name<-result$x_name #from function1
    x_name<-c(fixed_effect, random_effect) #this name is for prob
    
    
    x_name_rand <- outer(random_effect, random_effect, FUN = "paste0") #add on Feb 2022
    x_name_full <- c(fixed_effect, diag(x_name_rand), x_name_rand[lower.tri(x_name_rand, diag = F)]) #add on Feb 2022, the order of first 13 should be matached to N
    # set.seed(seed)
    
    v            <-dim(data[,var_model_selection_x])[2]
    name_res=as.data.frame(c(fixed_effect,random_effect))
    names(name_res)="name"
    fre_list=name_res # list of number of variables selected
    
    check_name <- x_name %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = gsub("\\d{1}$","\\" ,name),
      )
    
    check_name_full <-  x_name_full %>% as.data.frame() %>%
      dplyr::rename(name = ".") %>%
      tibble::rownames_to_column(., "number") %>%
      dplyr::mutate(
        V = paste0("V", number),
        no_number = sub("(.).*", "\\1",  name),
      )
    
    check_namename=as.character(check_name$name)
    
    PO_list <- list() #for probability result table
    CC_list <- list() #for coefficient result table
    
    M=M%>%t()%>%as.data.frame()
    M<-data.table::setnames(M, old=colnames(M), new=check_namename) %>% as.matrix()
    
    #IMPORTANT!!!
    converge_bt  <- as.numeric(dim(na.omit(M))[1]) #any missing indicated non-converge from function1
    
    #update calculate the importance index by using actual frequency
    if(use_fre_fun1){
      
      Impms        <-as.data.frame(colSums(abs(M), na.rm = T))%>%
        mutate(name=check_namename)
      colnames(Impms)=c("sum","name")
      Impms=
        Impms%>%
        left_join(result$fre,by="name")%>%
        dplyr::select(name,sum,total)%>%
        rowwise()%>%
        mutate(prob=sum/total)%>%
        mutate(prob=if_else(total==0,0,prob))%>%
        ungroup()%>%
        dplyr::select(name,prob)
      
      
      
      nms=Impms$name
      Impms=as.matrix(Impms[,2])
      rownames(Impms)=nms
      
    } else{ # using total number of bootstrapping
      
      Impms        <-as.matrix(colSums(abs(M), na.rm = T)/converge_bt )
    }
    
    
    
    
    
    
    
    # N    <-array(0,dim=c(bt*im,v,1))
    # NN   <-matrix(0,bt*im,v+rlcov) #need to include the covariance, # add on Feb 2022
    #NN <- data.table::setnames(NN, old=names(NN), new=check_name_full$name)
    
    
    
    for (s in 1:im){
      if (verbose) message(sprintf("Starting imputation % 3i of %i", s, im))
      MM           <- data%>%filter(.imp==s)%>%dplyr::select(var_list)
      
      y            <- as.matrix(MM[,ycol])
      scale_x_name <- setdiff(var_list_x, no_scale)
      xs <- MM %>%
        dplyr::select(all_of(var_list_x) )%>%
        mutate_at(., vars(scale_x_name), as.numeric) %>%
        mutate_at(., vars(scale_x_name), scale)
      j=1
      fc2=0
      while(j<=bt){
        if (verbose) message(sprintf("Starting bootstrap % 3i of %i", j, bt))
        unique_id2 <- as.data.frame(unique(data[, cluster_id]))
        names(unique_id2)="group_id"
        
        random_unique <- as.numeric(sample(unique_id2$group_id, round(dim(unique_id2)[1])*random_draw_sample,  replace = T))
        random_unique1<- random_unique %>%
          as.data.frame()%>%
          dplyr::rename(cluster_id = ".") %>%
          dplyr::arrange(cluster_id) %>%
          dplyr::mutate(
            cluster_id = as.numeric(cluster_id),
          ) %>%
          dplyr::group_by(cluster_id) %>%
          dplyr::mutate(
            time = as.numeric(rep(1:dplyr::n()))
          ) %>%
          dplyr::summarise(
            weight = max(time) #weight is to summarize how many times of being sampled
          ) %>%
          dplyr::ungroup()
        
        
        IVs2_x <- as.data.frame(xs) %>%
          dplyr::rename(time = time,
                        cluster_id = cluster_id)
        
        
        
        
        smp2 <- IVs2_x %>% #update the id infor
          dplyr::select(time, cluster_id) %>%
          dplyr::right_join(random_unique1, by="cluster_id") %>%  #depends on the random_unique1
          dplyr::group_split(cluster_id) %>%
          purrr::map(function(data){
            rep(data$time, times=max(data$weight)) %>% #replicate by weight
              as.data.frame() %>%
              dplyr::rename(
                time = "."
              ) %>%
              dplyr::mutate(cluster_id = rep(unique(data$cluster_id))) %>%
              dplyr::group_by(time) %>%
              dplyr::mutate(
                seq = rep(1:dplyr::n()),
                b_cluster_id = as.integer(paste0(cluster_id, seq))
              ) %>%
              dplyr::select(
                time, b_cluster_id, cluster_id
              )}) %>% do.call(rbind.data.frame, .)%>%
          group_by(b_cluster_id)%>%
          mutate(cluster_id_rpql=cur_group_id())
        
        if (family == "gaussian"){
          outcome2 <- MM %>%
            dplyr::select(cluster_id, time, ycol) %>%
            purrr::map_dfc(., as.numeric) %>% #for continuous outcome, can do this
            dplyr::mutate_at(., vars(ycol), scale) %>% #DV needs to be scale
            dplyr::arrange(cluster_id) %>%
            dplyr::rename(time = time, cluster_id = cluster_id)
          
        }else{
          outcome2 <- MM %>%
            dplyr::select(cluster_id, time, ycol) %>%
            dplyr::mutate_at(., vars(c(cluster_id, time)), as.numeric) %>%
            arrange(cluster_id) %>%
            dplyr::rename(time = time, cluster_id = cluster_id)
        }
        
        
        IVs2 <- smp2 %>% dplyr::left_join(IVs2_x, by=c("cluster_id", "time"))
        DV2  <- smp2 %>% dplyr::left_join(outcome2, by=c("cluster_id", "time"))
        
        #if(dim(IVs2)[1] != dim(DV2)[1]){warning("The number of rows of x is not as same as the number of the outcome")}
        
        x1   <- IVs2 %>% dplyr::arrange(cluster_id, time) %>% as.matrix()
        y1   <- DV2 %>% dplyr::arrange(cluster_id, time) %>%
          dplyr::select(b_cluster_id, time, ycol) %>%
          as.matrix()
        
        x2_data <- IVs2 %>%
          dplyr::ungroup() %>%
          dplyr::left_join(DV2, by=c("cluster_id", "time", "b_cluster_id","cluster_id_rpql")) %>%
          dplyr::arrange(cluster_id, b_cluster_id,cluster_id_rpql, time) %>%
          dplyr::select(-time)
        
        #the following condition is for all continuous predictors random draw
        x2_data_selective_x_only <- x2_data %>%
          dplyr::select(all_of(var_model_selection_x))%>%
          dplyr::mutate_all(., as.numeric)%>%
          # change 9/5 choose w1 and z1 also
          dplyr::select(starts_with("x")) #, -x1) #CHANGE for package, specify the intercept
        
        x2_data_selective_z_only <- x2_data %>%
          dplyr::select(all_of(var_model_selection_x))%>%
          dplyr::mutate_all(., as.numeric)%>%
          # change 9/5 choose w1 and z1 also
          dplyr::select(starts_with("z"))#,-z1 ) #CHANGE for package, specify the intercept
        
        x2_lmer_use <- x2_data %>%
          dplyr::ungroup() %>%
          dplyr::select(all_of(var_model_selection_x), ycol, b_cluster_id,cluster_id_rpql) %>%
          dplyr::mutate_all(., as.numeric) #assume all continuous predictors
        
        ## method 2 half sample from fixed effect half from random effect
        # fix
        name_fix    <- gtools::mixedsort(grep("x",check_name$name,value=T)) #update 06/30, replace rownames function
        # loc_x1=which(rownames(Impms)=="x1")
        # name_fix=name_fix[name_fix!="x1"]
        fix_ind=grep("x",rownames(Impms))
        # fix_ind=fix_ind[fix_ind!=loc_x1]
        impms_fix=Impms[fix_ind]
        # random
        name_random <- gtools::mixedsort(grep("z",rownames(Impms),value=T))
        # loc_z1=which(rownames(Impms)=="z1")
        # name_random=name_random[name_random!="z1"]
        random_ind=grep("z",rownames(Impms))
        # random_ind=random_ind[random_ind!=loc_z1]
        impms_random=Impms[random_ind]
        
        
        size_fixed  <- q2_fix #q2 also needs to set up a default value for package
        size_random <- q2_random
        
        w1          <- sample(name_fix,size=size_fixed,prob=impms_fix)#update by name rather than order
        w1=gtools::mixedsort(w1)
        # w1=c("x1",w1)
        # for (i in 1:rep){
        w2          <- sample(name_random,size=size_random,prob=impms_random)
        w2=gtools::mixedsort(w2)
        
        fre=cbind(as.data.frame(c(w1,w2)),rep(1,q2_fix+q2_random))
        names(fre)=c("name","fre")
        
        fre <- name_res %>%
          dplyr::left_join(fre, by="name") %>%
          dplyr::mutate(
            unselect = rep(0),
            fre = dplyr::coalesce(as.numeric(fre), unselect)
          ) %>%
          dplyr::select(-unselect)
        
        
        x2_data2_x          <- cbind(x2_data_selective_x_only[, w1])
        x2_data2_z         <- cbind(x2_data_selective_z_only[, w2])
        
        x2_data3 <- names(cbind(x2_data_selective_x_only,(x2_data_selective_z_only)))%>%
          
          
          as.data.frame() %>%
          dplyr::rename(name = ".") %>%
          tibble::rownames_to_column(., "number")
        
        
        x2_data4 <- names(cbind(x2_data2_x,x2_data2_z)) %>%
          as.data.frame() %>%
          dplyr::rename(name = ".") %>%
          left_join(x2_data3,by="name")
        
        
        
        
        # formula2        <- paste0(ycol, "~", paste0(colnames(x2_data2), collapse = "+"),"-1", "+(1|b_cluster_id)")
        b_fix_effect <- x2_data4 %>% dplyr::filter(name %in% fixed_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
        b_random_effect <- x2_data4 %>% dplyr::filter(name %in% random_effect) %>% dplyr::select(name) %>% dplyr::pull(name)
        
        x <- as.matrix(x2_lmer_use[, b_fix_effect])
        y <- x2_lmer_use[, ycol] %>% dplyr::pull(ycol)
        z <- as.matrix(x2_lmer_use[, b_random_effect])
        
        
        list2 <- list(id = list(cluster_id_rpql = x2_lmer_use$cluster_id_rpql), X=x, Z=list(cluster_id_rpql=z), Y=y)
        
        if(pen.type=="adl" & use_lmer==T){ # use lmer to calculate initial value
          if (family == "gaussian"){#for lmer DV continue
            fit2  <- lme4::lmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use)
          } # # gaussian
          if (family %in% "binomial"){#for glmer DV binary
            fit2  <- lme4::glmer(y ~ x-1 + (z-1 | cluster_id_rpql), data=x2_lmer_use,family = "binomial")} #binomial
          xx    <- x %>% as.data.frame()
          
          fit_sat2 <- build.start.fit(fit2, gamma = 2, cov.groups = NULL) #no categorical predictor --> CHANGE for package
          fit_final2 <- tryCatch(rpqlseq(y = list2$Y, X = list2$X, Z = list2$Z, id = list2$id,
                                         lambda = lam, pen.type = "adl", hybrid.est = TRUE,
                                         pen.weights = fit_sat2$pen.weights, start = fit_sat2,restarts = 50),
                                 error=function(e){return(NA)})
        } else if(pen.type=="adl" & use_lmer==F){ # do not use lmer
          weight_f=cbind(names=rownames(Impms),as.data.frame(Impms))
          colnames(weight_f)[2]="prob"
          weight_f=weight_f%>%
            filter(names %in% colnames(x))%>%
            dplyr::select(prob)%>%
            t()%>%as.data.frame()%>%
            dplyr::select(colnames(x))%>%
            as.vector()
          weight_f=as.numeric(weight_f)^-1
          # weight_f=as.vector(t(weight_f))
          weight_r=cbind(names=rownames(Impms),as.data.frame(Impms))
          colnames(weight_r)[2]="prob"
          weight_r=weight_r%>%
            # rename(prob=V1)%>%
            filter(names %in% colnames(z))%>%
            dplyr::select(prob)%>%
            t()%>%as.data.frame()%>%
            dplyr::select(colnames(z))
          weight_r=weight_r^-1
          weight_r=as.vector(t(weight_r))
          
          wt=list(fixed=weight_f,random=list(cluster_id_rpql=weight_r))
          
          
          fit_final2 <- tryCatch(rpqlseq(y = list2$Y, X = list2$X, Z = list2$Z, id = list2$id,
                                         lambda = lam, pen.type = "adl", hybrid.est = TRUE,pen.weights=wt,
                                         restarts = 50),
                                 error=function(e){return(NA)})
          
        } else if(pen.type=="lasso"){
          fit_final2 <- tryCatch(rpqlseq(y = list1$Y, X = list1$X, Z = list1$Z, id = list1$id,
                                         lambda = lam, pen.type = "lasso", hybrid.est = TRUE,
                                         restarts = 50),
                                 error=function(e){return(NA)})
        }
        
        
        
        if(all(is.na(fit_final2))){j=j} else{
          j=j+1
        }
        # if(anyNA(fit_final2)){fc2=fc2+1 }
        # if(fc2>bt*0.05){break} # if 5% bootscraping cannot converage, break the loop and give NA
        
        if(!anyNA(fit_final2)){
          
          final_result <- fit_final2$best.fit[[ci_criteria]] #save all results
          result2      <- summary(fit_final2$best.fit[[ci_criteria]]) #save summary result
          fixcoef2     <- as.data.frame(result2$fixef) %>%  tibble::rownames_to_column(., "name") %>% dplyr::rename(value = `result2$fixef`)
          randcov2     <- as.matrix(result2$ran.cov$cluster_id_rpql)
          randcoef2    <- diag(randcov2)%>%
            as.data.frame() %>%
            cbind(name=as.character(b_random_effect),.)%>%
            dplyr::rename(value = ".")
          randcoef3    <-  randcov2[lower.tri(randcov2, diag = T)]
          togethercoef <- rbind(fixcoef2, randcoef2)
          fre_list=left_join(fre_list,fre,by="name")
          
          name_res=as.data.frame(c(fixed_effect,random_effect))
          names(name_res)="name"
          x2_data5_cc <- name_res %>% #this coefficient table does not include the covariance values
            dplyr::left_join(togethercoef, by="name") %>%
            dplyr::mutate(
              unselect = rep(NA),
              value = dplyr::coalesce(as.numeric(value), unselect)
            ) %>%
            dplyr::select(name, value) %>%
            # dplyr::arrange(name) %>%
            t() %>%
            as.data.frame() %>%
            header.true()
          
          
          x2_data5_po <- name_res %>%
            dplyr::left_join(togethercoef, by="name") %>%
            #purrr::map_dfc(., as.numeric) %>%
            dplyr::mutate(
              unselect = rep(0),
              value1 = dplyr::coalesce(as.numeric(value), unselect),
              value = dplyr::if_else(value1 != 0, 1, 0)
            ) %>%
            dplyr::select(name, value) %>%
            # dplyr::arrange(name) %>%
            t() %>%
            as.data.frame()%>%
            header.true()
          
        }else{
          final_result <- list() #null list
          # fc2=fc2+1
          x2_data5 <- x2_data3 %>%
            dplyr::mutate(
              value = rep(NA_real_)
            ) %>%
            dplyr::select(name, value) %>%
            dplyr::arrange(name) %>%
            t() %>%
            as.data.frame() %>%
            header.true()
          x2_data5_cc <- x2_data5
          x2_data5_po <- x2_data5
          
          
        }#end of managing result2
        
        #rm(random_unique, random_unique1, w1,w2) #delete the random draw id
        
        PO_list[[j+(s-1)*bt]] <- x2_data5_po
        CC_list[[j+(s-1)*bt]] <- x2_data5_cc
        result_matrix[[(s-1)*bt+j]] <- final_result
        
        # }
      } #end of boostrapping
      
    }#end of imputation
    PO <- do.call(rbind, PO_list) %>% as.data.frame() %>% purrr::map_dfc(., as.character)%>%purrr::map_dfc(., as.numeric)
    CC <- do.call(rbind, CC_list) %>% as.data.frame() %>% purrr::map_dfc(., as.character)%>%purrr::map_dfc(., as.numeric)
    converge_bt  <- as.numeric(dim(na.omit(PO))[1]) #any missing indicated non-converge from function1
    
    
    
    
    
    fre_res=fre_list%>%
      rowwise()%>%
      mutate(total=sum(c_across(where(is.numeric)),na.rm=T))%>%
      ungroup()%>%
      dplyr::select(name,total)
    
    
    
    coefficient_matrix   <- as.matrix(abs(colSums(CC, na.rm = T)/converge_bt) )
    probability_matrix   <- as.matrix(colSums(PO, na.rm = T)/converge_bt)
    
    
    
    #make a little bit improvement:
    ds_prob <- as.data.frame(probability_matrix) %>%
      data.table::setDT(keep.rownames = T) %>%
      dplyr::rename(prob = `V1`, name = rn) %>%
      dplyr::arrange(desc(prob))
    
    ds_coef <- as.data.frame(coefficient_matrix) %>%
      data.table::setDT(keep.rownames = T) %>%
      dplyr::rename(estimate = `V1`, name = rn) %>%
      dplyr::arrange(desc(abs(estimate)))
    
    
    
    
    
    
    round2_res<-
      list(probability = list(matrix = PO, summary = ds_prob),
           coefficient = list(matrix = CC, summary = ds_coef),
           raw_rpql_result  =  result_matrix,
           frequency=fre_res,
           fre_matrix=fre_list)
    
    
    
    
    
    return(round2_res)
    
  }
  
}#end of function2


