library(mvtnorm)
library(magrittr) #Pipes
library(VIM)
library(lme4)
library(lmerTest)
library(officer)
library(lares)
library(dplyr)
library(MASS)
library(missMethods) #generating missing values for simulation
library(Matrix)
library(rpql)


generate_glmm <-function(beta, time_var=T, n_time=NULL, n_group, sigma_value, mean_value, between_sub_var, time_matrix, x_z_same = F, xz_cor = 0.5,
                         group_size_low=30,group_size_upper=30,multi_cor=0.2,family="gaussian",mix=F){
  q <- length(diag(between_sub_var))
  p <- length(beta)

  beta   <- as.matrix(beta, nrow = length(beta), ncol=1)
  beta_t <- t(beta)

  # update equal/unequal group_size
  if(group_size_low==group_size_upper){# equal group_size
    group_size=rep(group_size_low,n_group)
  } else { #unequal group_size
    group_size=sample(c(group_size_low:group_size_upper),n_group,replace = T)}

  if(time_var==F){#cross-sectional
    check_xz_dim <-  abs(p-q) #only when both x and z to be the same in dim,  can let X and Z to be same
    if(x_z_same == T & check_xz_dim==0){
      sig=diag(p-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig[6:7,6:7]=multi_cor

      if(mix){
        sig[6:7,6:7]=multi_cor*(-1)
        sig[4,2]=sig[2,4]=xz_cor*(-1)
      }

      diag(sig)=1


      H        <- abs(outer(1:(p-1), 1:(p-1), "-"))
      I        <- abs(outer(1:(q-1), 1:(q-1), "-"))

      x <- rmvnorm(sum(group_size),rep(0,p-1),sigma=0^H)


      x_m_new <- x
      c1 <- var(x_m_new)
      # cholesky decomposition to get independence
      chol1 <- solve(chol(c1))
      newx <-  x_m_new %*% chol1
      chol2 <- chol(sig)
      x_m <- newx %*% chol2

      x=cbind(1,x_m)
      z <- x

    }else{ #X and Z are diff
      H <- abs(outer(1:(p-1), 1:(p-1), "-"))
      I <- abs(outer(1:(q-1), 1:(q-1), "-"))
      # x <- cbind(1,rmvnorm(n_group*n_time,rep(0,p-1),sigma=xz_cor^H))
      # z <- cbind(1,rmvnorm(n_group*n_time,rep(0,q-1),sigma=xz_cor^I))

      sig=diag(p-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig[1,3]=sig[2,4]=sig[3,1]=sig[4,2]=multi_cor
      sig[6:7,6:7]=multi_cor
      diag(sig)=1

      # if(mix){
      #   sig[3,1]=sig[4,2]=sig[1,3]=sig[2,4]=xz_cor*(-1)
      # }

      sig_r=diag(q-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig_r[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig_r[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig[1,3]=sig[2,4]=sig[3,1]=sig[4,2]=multi_cor
      sig_r[6:7,6:7]=multi_cor
      diag(sig_r)=1

      x <- rmvnorm(sum(group_size),rep(0,p-1),sigma=0^H)


      x_m_new <- x
      c1 <- var(x_m_new)
      # cholesky decomposition to get independence
      chol1 <- solve(chol(c1))
      newx <-  x_m_new %*% chol1
      chol2 <- chol(sig)
      x_m <- newx %*% chol2

      x=cbind(1,x_m)

      z <- rmvnorm(sum(group_size),rep(0,q-1),sigma=0^H)


      z_m_new <- z
      # find the current correlation matrix
      c1z <- var(z_m_new)
      # cholesky decomposition to get independence
      chol1z <- solve(chol(c1z))
      newz <-  z_m_new %*% chol1
      chol2z <- chol(sig_r)
      z_m <- newz %*% chol2

      z=cbind(1,z_m)


    }


  }
  else{#when time_var is true --> longitudinal set the number of wave in n_time variable n_time=group_size_low=group_size_upper

    #when time_var is true --> longitudinal set the number of wave in n_time variable n_time=group_size_low=group_size_upper
    if(x_z_same==F){
      n_time   <- length(time_matrix)
      #x_z_same <- F #force X and Z are diff
      H        <- abs(outer(1:(p-2), 1:(p-2), "-"))
      I        <- abs(outer(1:(q-2), 1:(q-2), "-"))

      # update 8/18 change the correlation matrix for x
      sig=diag(p-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig[6:7,6:7]=multi_cor
      diag(sig)=1


      # change the correlation matrix for z
      sig_r=diag(q-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig_r[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig_r[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig_r[6:7,6:7]=multi_cor
      diag(sig_r)=1

      if(p < 3){x_m <- rep(time_matrix, n_group)
      }else{
        # update 8/18 change the correlation between time variable and other variable
        x_m <- rmvnorm(n_group*n_time,rep(0,p-2),sigma=0^H)
        x1=rep(time_matrix, n_group)
        x_m_new <- cbind(scale(x1),x_m)
        # find the current correlation matrix
        c1 <- var(x_m_new)
        # cholesky decomposition to get independence
        chol1 <- solve(chol(c1))
        newx <-  x_m_new %*% chol1
        chol2 <- chol(sig)
        x_m <- newx %*% chol2 * sd(x1) + mean(x1)

      } #assuming the first two columns are constant
      x <- cbind(1, x_m)

      if(q < 3){z_m <- rep(time_matrix, n_group)
      }else{ # update 8/18 change the correlation between time variable and other variable
        z_m <- rmvnorm(n_group*n_time,rep(0,q-2),sigma=0^I)
        z1=rep(time_matrix, n_group)
        z_m_new <- cbind(scale(z1),z_m)
        # find the current correlation matrix
        c1z <- var(z_m_new)
        # cholesky decomposition to get independence
        chol1z <- solve(chol(c1z))
        newz <-  z_m_new %*% chol1z
        chol2z <- chol(sig_r)
        z_m <- newz %*% chol2z * sd(z1) + mean(z1)
      }

      z <- cbind(1, z_m)

    } else{

      n_time   <- length(time_matrix) #
      #x_z_same <- F #force X and Z are diff
      H        <- abs(outer(1:(p-2), 1:(p-2), "-"))
      # I        <- abs(outer(1:(q-2), 1:(q-2), "-"))

      # update 8/18 change the correlation matrix for x
      sig=diag(p-1)*(1-xz_cor)+xz_cor
      # correlation x2-x4-->0.9
      sig[1:2,1:2]=multi_cor
      # correlation x5,x6-->0.9
      sig[3:4,3:4]=multi_cor
      # x7 x8 -->0.9
      sig[6:7,6:7]=multi_cor
      diag(sig)=1


      if(p < 3){x_m <- rep(time_matrix, n_group)
      }else{
        # update 8/18 change the correlation between time variable and other variable
        x_m <- rmvnorm(n_group*n_time,rep(0,p-2),sigma=0^H)
        x1=rep(time_matrix, n_group)
        x_m_new <- cbind(scale(x1),x_m)
        # find the current correlation matrix
        c1 <- var(x_m_new)
        # cholesky decomposition to get independence
        chol1 <- solve(chol(c1))
        newx <-  x_m_new %*% chol1
        chol2 <- chol(sig)
        x_m <- newx %*% chol2 * sd(x1) + mean(x1)

      } #assuming the first two columns are constant
      x <- cbind(1, x_m)
      z=x}
  }#end of time_var
  x_name <- paste0("x", c(1:p)) %>% c()
  z_name <- paste0("z", c(1:q)) %>% c()
  ds_x <- as.data.frame(x) %>% data.table::setnames(., old=names(.), new=x_name)
  ds_z <- as.data.frame(z) %>% data.table::setnames(., old=names(.), new=z_name)


  b0 <- MASS::mvrnorm(
    n  = n_group,
    mu = rep(0, q),
    Sigma = between_sub_var
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      group_id = rep(1:n_group)
    ) %>%
    tidyr::pivot_longer(cols = -group_id) %>%
    dplyr::select(value) %>%
    as.matrix()


  Z <- matrix(0, nrow=sum(group_size), ncol=q*n_group)
  re=0
  ri=0
  for (i in 1:n_group) {

    ri <- re+1
    re <- ri+group_size[i]-1

    ci <- (i-1)*q+1
    ce <- (i-1)*q+q
    #
    Z[ri : re, ci: ce] <- z[ri : re, ]
  }


  E_matrix <- list()
  for (i in 1:n_group) {
    E_matrix[[i]] <- rnorm(group_size[i], mean = 0, sd = sigma_value) %>%
      as.data.frame() %>%
      dplyr::rename(epsilon = ".") %>%
      dplyr::mutate(
        within_sub_id = rep(1:group_size[i]),
        group_id = rep(i)
      )
  }

  E_matrix <- do.call(rbind, E_matrix)
  E <- as.matrix(E_matrix$epsilon)

  if(family=="gaussian"){ # y is continue
    Y <- x%*%beta+Z%*%b0+E
    y <- as.data.frame(Y) %>% dplyr::rename(y= value)}

  if(family=="binomial"){ # y is binary
    # y create by the gendat.glmm by the rpql library
    Y <- rpql::gendat.glmm(id = list(cluster = E_matrix$group_id), X = ds_x, beta = beta,
                           Z = list(cluster = ds_z), D = list(cluster = between_sub_var), phi = 1, family = binomial())

    y=Y$y }

  # a=lme4::VarCorr(fit)
  # diag(a[[1]])





  ds <- cbind(ds_x, ds_z, E_matrix, y) %>% as.data.frame()

  return(ds)
}#end of function


calculate_cov <- function(tao, cor_coeff){
  q <- length(tao)
  n_off_dia <- q*(q-1)/2
  n_nonzero_tao <- length(which(tao!=0))
  correlation <- matrix(0, nrow = q, ncol = q)
  cor_matrix <- diag(n_nonzero_tao)
  off_cor <- cor_coeff
  cor_matrix[lower.tri(cor_matrix)] <- off_cor[which(off_cor!=0)]
  cor_matrix[upper.tri(cor_matrix)] <- t(cor_matrix)[upper.tri(t(cor_matrix))]
  correlation[1:n_nonzero_tao, 1:n_nonzero_tao]<-cor_matrix

  tao_t <- t(tao)
  var_cov0 <- tao%*%tao_t
  var_cov <- correlation*var_cov0
}


# ---- generate-complete-data --------------------------------------------------

beta=c(0.1, 0.5, 0.4, 0.3, 0.2,0,0, 0,0,0,rep(0, 20))


tao_use         <- c(3,2,1, rep(0, 27))

cor_coeff_use   <- c(0.8, 0.2, 0.5, rep(0, 27))

between_sub_var <- calculate_cov(tao_use,  cor_coeff_use)
sigma_value     <- 1
mean_value      <- 0
n_group         <- 50 #group size
n_time          <-30 #within group: cluster-size
group_size_low=30
group_size_upper=30
family="gaussian"
time_var        <- F
xz_cor          <- 0.05 #


generate_list <- replicate(n=500, generate_glmm(beta, time_var=F, n_time, n_group, sigma_value,
                                                mean_value, between_sub_var, time_matrix,
                                                x_z_same = T, xz_cor=0.05,
                                                group_size_low,multi_cor =-0.05,
                                                group_size_upper,family,mix=T),
                           simplify = F)




