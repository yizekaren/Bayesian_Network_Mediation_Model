
BNMM <- function(y, x, d, z, num_visit, mcmc_samples = 1000, num_subnetwork){
  # y: vector of outcome (n*1)
  # x: vector of treatment (n*1)
  # d: design matrix for covariate (n*(number of covariate + 1))
  # z: array of connectivity matrices (n*num_visit_maximum*num_region*num_region)
  # num_visit: number of visit for each subject (n*1)
  # mcmc_samples: number of MCMC iterations
  # num_subnetwork: number of subnetworks
  
  n <- length(y)
  pd <- dim(d)[2] 
  num_visit_maximum <- dim(z)[2]
  num_region <- dim(z)[3]
  
  # initialization
  num_mediator <- ((num_subnetwork+1)*num_subnetwork)/2
  gamma_indicator_y <- matrix(0, nrow = mcmc_samples, ncol = num_mediator)
  gamma_indicator_theta <- matrix(0, nrow = mcmc_samples, ncol = num_mediator)
  p_gamma_y <- 0.5
  p_gamma_theta <- 0.5
  
  gamma_indicator_y[1,] <- sample(c(0, 1), size = num_mediator, replace = T, prob = p_gamma_y)
  gamma_indicator_theta[1,] <- sample(c(0, 1), size = num_mediator, replace = T, prob = p_gamma_theta)
  
  gamma_temp_y <- gamma_indicator_y[1,]
  gamma_temp_theta <- gamma_indicator_theta[1,]
  
  gamma_indicator_temp_decision_y<- rep(0, times = num_mediator)
  gamma_indicator_temp_decision_theta<- rep(0, times = num_mediator)
  
  beta_x <- matrix(0, nrow = mcmc_samples, ncol = pd)
  beta_x[1, ] <- rnorm(pd)
  sigma2alpha <- 10^2
  
  beta_m <- matrix(0,nrow = mcmc_samples, ncol = num_mediator)
  beta_m[1, which(gamma_temp_y==1)] <- rnorm(sum(gamma_temp_y==1))
  
  beta_z <- rep(0, times = mcmc_samples)
  beta_z[1] <- rnorm(1)
  sigma20 <- 10^2 # sigma_{zy}^2
  
  alpha_x <- array(0, dim = c(mcmc_samples, num_mediator, pd)) #w_matrix
  sigma2m <- 10^2 # sigma_{xm}^2
  
  alpha_z <- matrix(0, nrow = mcmc_samples, ncol = num_mediator)
  alpha_z[1, which(gamma_temp_theta==1)] <- rnorm(sum(gamma_temp_theta==1))
  
  sigma2a <- rep(0, times = mcmc_samples) 
  sigma2a[1] <- 1^2
  a3 <- 0.1
  b3 <- 0.1
  
  sigma2b <- rep(0, times = mcmc_samples)
  sigma2b[1] <- 1^2
  a4 <- 0.1
  b4 <- 0.1
  
  sigma2w <- rep(0, times = mcmc_samples) # sigma_1^2
  sigma2w[1] <- 1^2
  a2 <- 2
  b2 <- 2
  
  sigma2 <- rep(0, times = mcmc_samples) # sigma_2^2
  sigma2[1] <- 1^2
  a1 <- 1
  b1 <- 1
  
  sigma2k <- rep(0, times = mcmc_samples)
  sigma2k[1] <- 1^2
  ak <- 1
  bk <- 1
  
  sigma2qr_matrix <- array(0, dim = c(mcmc_samples, num_subnetwork, num_subnetwork)) # sigma_qr
  sigma2qr_matrix[1,,] <- matrix(1, nrow = num_subnetwork, ncol = num_subnetwork)
  a0 <- 1
  b0 <- 1
  
  paic <- matrix(0, nrow = mcmc_samples, ncol = num_subnetwork)
  paic_sp_pos <- rep(0, times = num_subnetwork)
  dir_vec <-  rep(3, times = num_subnetwork)
  dir_vec_pos<- rep(0, times = num_subnetwork)
  paic[1,]<- rdirichlet(1, dir_vec)
  
  c_latent_ini<- matrix(0, nrow = num_region, ncol = num_subnetwork)
  for (i in 1:num_region){
    c_latent_ini[i,] <- t(rmultinom(1, size = 1, prob = paic[1,]))
  }
  
  c_latent <- array(0, dim = c(mcmc_samples, num_region, num_subnetwork))
  c_latent[1,,] <- c_latent_ini
  chosen_region<- seq(from = 1, to = num_region, by = 5)
  for (i in 1:length(chosen_region)){
    c_latent[1, chosen_region[i], ]<- t(rmultinom(1, size = 1, prob = paic[1,]))
  }
  
  mik_qr_ini <- array(0, dim=c(n, num_visit_maximum, num_mediator))
  for (i in 1:n){
    for (j in 1:num_visit[i]){
      for (s_subnetwork in 1:num_subnetwork){
        for (t_subnetwork in s_subnetwork:num_subnetwork){
          s_region<- which(c_latent[1,,s_subnetwork]==1)
          t_region<- which(c_latent[1,,t_subnetwork]==1)
          min_l_p<- min(s_subnetwork,t_subnetwork)
          max_l_p<- max(s_subnetwork,t_subnetwork)
          k_cal<- (min_l_p-1)*num_subnetwork+ max_l_p- (min_l_p*(min_l_p-1)/2)
          mik_qr_ini[i,j,k_cal]<-  mean(z[i,j,s_region,t_region], na.rm = T)
        }
      }
    }
  }
  
  mik_qr<- array(0, dim = c(mcmc_samples, n, num_visit_maximum, num_mediator))
  mik_qr[1,,,] <- mik_qr_ini
  
  mi_qr_ini <- array(0, dim = c(n, num_mediator))
  for (j in 1:n){
    for (k in 1:num_mediator){
      mi_qr_ini[j,k]<- mean(mik_qr_ini[j,(1:num_visit[j]),k]) 
    }
  }
  
  mi_qr <- array(0, dim = c(mcmc_samples, n, num_mediator))
  mi_qr[1, , ] <- mi_qr_ini
  
  
  final_class <- rep(0, times=num_region)
  
  
  #### MCMC running ####
  for (i in 2:mcmc_samples) {
    print(paste0("Iteration: ", i))
    
    #gamma_for_y
    tp_num<- 0
    fp_num<- 0
    
    var_y<- diag(sigma2w[(i-1)],n)
    
    for (k in 1:num_mediator){
      #loglik for 1
      gamma_temp_y[k]<- 1
      mean_y_1<- d%*%beta_x[(i-1),]+mi_qr[(i-1),,]%*%(beta_m[(i-1),]*gamma_temp_y)+beta_z[(i-1)]*x
      val_1<- dmvnorm(t(y),mean_y_1,var_y,log = T) + log(p_gamma_y)
      
      gamma_temp_y[k]<- 0
      mean_y_0<- d%*%beta_x[(i-1),]+mi_qr[(i-1),,]%*%(beta_m[(i-1),]*gamma_temp_y)+beta_z[(i-1)]*x
      val_0<- dmvnorm(t(y),mean_y_0,var_y,log = T) + log(1-p_gamma_y)
      
      max_num<- max(val_1,val_0)
      val_1<- val_1- max_num
      val_0<- val_0- max_num
      
      val_1<- exp(val_1)/sum(exp(val_1)+exp(val_0))
      gamma_indicator_y[i,k]<- rbinom(1,1,val_1)
      gamma_temp_y[k]<-  gamma_indicator_y[i,k]
      
      #compute result for gamma[k]
      if (mean(gamma_indicator_y[1:i,k])>0.5){
        gamma_indicator_temp_decision_y[k]<- 1
      } else {
        gamma_indicator_temp_decision_y[k]<- 0
      }
      
    }
    
    # gamma_for_theta (mi_qr)
    var_thetak<- diag(sigma2[(i-1)],n)
    
    for (k in 1:num_mediator){
      #loglik for 1
      gamma_temp_theta[k]<- 1
      mean_thetak_1<- d%*%alpha_x[(i-1),k,]+ (alpha_z[(i-1),k]*gamma_temp_theta[k])*x
      val_1<- dmvnorm(t(mi_qr[(i-1),,k]),mean_thetak_1,var_thetak,log=T)+log(p_gamma_theta)
      
      #loglik for 0
      gamma_temp_theta[k]<- 0
      mean_thetak_0<- d%*%alpha_x[(i-1),k,]+ (alpha_z[(i-1),k]*gamma_temp_theta[k])*x
      val_0<- dmvnorm(t(mi_qr[(i-1),,k]),mean_thetak_0,var_thetak,log=T)+log(1-p_gamma_theta)
      
      max_num<- max(val_1,val_0)
      val_1<- val_1- max_num
      val_0<- val_0- max_num
      
      val_1<- exp(val_1)/sum(exp(val_1)+exp(val_0))
      gamma_indicator_theta[i,k]<- rbinom(1,1,val_1)
      gamma_temp_theta[k]<-  gamma_indicator_theta[i,k]
      
      #compute result for gamma[k]
      if (mean(gamma_indicator_theta[1:i,k])>0.5){
        gamma_indicator_temp_decision_theta[k]<- 1
      } else {
        gamma_indicator_temp_decision_theta[k]<- 0
      }
      
      if ((gamma_indicator_temp_decision_y[k]*gamma_indicator_temp_decision_theta[k])==1){
        if (gamma_indicator_true[k]==1){
          tp_num<- tp_num+1
        } else {
          fp_num<- fp_num+1
        }
      }
    }
    
    num_correct_mediator<- sum(gamma_indicator_temp_decision_y*gamma_indicator_temp_decision_theta==gamma_indicator_true)
    
    #beta_m
    gamma_chosen_y<- which(gamma_temp_y==1)
    gamma_not_chosen_y<- which(gamma_temp_y==0)
    
    if (length(gamma_not_chosen_y)>0){
      beta_m[i,gamma_not_chosen_y]<- mvrnorm(n=1,mu=rep(0,times=length(gamma_not_chosen_y)),Sigma=diag(sigma2a[(i-1)],length(gamma_not_chosen_y)))
    }
    
    if (length(gamma_chosen_y)>0){
      t_vec<- y- d%*%beta_x[(i-1),]-beta_z[(i-1)]*x 
      sigma_pos_inv<- (sigma2a[(i-1)]*t(mi_qr[(i-1),,gamma_chosen_y])%*%mi_qr[(i-1),,gamma_chosen_y]+diag(sigma2w[(i-1)],length(gamma_chosen_y)))/(sigma2w[(i-1)]*sigma2a[(i-1)])
      sigma_pos<- solve(sigma_pos_inv)
      mu_pos<- (sigma_pos%*%t(mi_qr[(i-1),,gamma_chosen_y])%*%t_vec)/sigma2w[(i-1)]
      
      beta_m[i,gamma_chosen_y]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
    }
    
    
    #alpha_z
    gamma_chosen_theta<- which(gamma_temp_theta==1)
    gamma_not_chosen_theta<- which(gamma_temp_theta==0)
    
    if (length(gamma_not_chosen_theta)>0){
      alpha_z[i,gamma_not_chosen_theta]<- mvrnorm(n=1,
                                                  mu=rep(0,times=length(gamma_not_chosen_theta)),
                                                  Sigma=diag(sigma2b[(i-1)],length(gamma_not_chosen_theta)))
    }
    
    if (length(gamma_chosen_theta)>0){
      sigma_pos_inv<- diag((sigma2b[(i-1)]*sum(x^2)+sigma2[(i-1)])/(sigma2[(i-1)]*sigma2b[(i-1)]),length(gamma_chosen_theta))
      sigma_pos<- solve(sigma_pos_inv)
      
      sum_beta_theta<- 0
      for (j in 1:n){
        t_vec<- mi_qr[(i-1),j,gamma_chosen_theta]- alpha_x[(i-1),gamma_chosen_theta,]%*%d[j,]
        sum_beta_theta<- sum_beta_theta+ x[j]*t_vec
      }
      
      mu_pos<- (sigma_pos%*%sum_beta_theta)/sigma2[(i-1)]
      
      alpha_z[i,gamma_chosen_theta]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
    }
    
    
    #beta_z
    sigma_pos_inv<- (sum(x^2)/sigma2w[(i-1)]+ (1/sigma20))
    sigma_pos<- 1/sigma_pos_inv
    
    t_vec<- y-d%*%beta_x[(i-1),]-mi_qr[(i-1),,]%*%(beta_m[i,]*gamma_indicator_y[i,])
    sum_alpha0<- sum(x*t_vec)
    
    mu_pos_alpha0<- (sigma_pos*sum_alpha0)/sigma2w[(i-1)]
    
    beta_z[i]<- rnorm(1,mean=mu_pos_alpha0,sd=sqrt(sigma_pos))
    
    #beta_x
    sigma_pos_inv<- (sigma2alpha*t(d)%*%d+diag(sigma2w[(i-1)],pd))/(sigma2w[(i-1)]*sigma2alpha)
    sigma_pos<- solve(sigma_pos_inv)
    
    t_vec<- y- mi_qr[(i-1),,]%*%(beta_m[i,]*gamma_indicator_y[i,])- x*beta_z[i]
    
    mu_pos<- (sigma_pos%*%t(d)%*%t_vec)/sigma2w[(i-1)]
    
    beta_x[i,]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
    
    #sigma2a
    sigma2a[i]<- 1/rgamma(1,(a3+(num_mediator/2)),(b3+(sum(beta_m[i,]^2)/2)))
    
    #sigma2b
    sigma2b[i]<- 1/rgamma(1,(a4+(num_mediator/2)),(b4+(sum(alpha_z[i,]^2)/2)))
    
    #sigma2w
    t_vec<- y-d%*%beta_x[i,]-mi_qr[(i-1),,gamma_chosen_y]%*%as.matrix(beta_m[i,gamma_chosen_y])-beta_z[i]*x
    sigma2w[i]<- 1/rgamma(1,(a2+ (n/2)),(b2+ (sum(t_vec^2)/2)))
    
    
    #alpha_x
    sigma_pos_inv<- (sigma2m*t(d)%*%d+diag(sigma2[(i-1)],pd))/(sigma2[(i-1)]*sigma2m)
    sigma_pos<- solve(sigma_pos_inv)
    for (k in 1:num_mediator){
      t_vec<- mi_qr[(i-1),,k]- ((alpha_z[i,k]*gamma_indicator_theta[i,k])*x)
      mu_pos<- (sigma_pos%*%t(d)%*%t_vec)/sigma2[(i-1)]
      alpha_x[i,k,]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
    }
    
    
    #mi_qr: update for each individual theta_i
    for (k in 1:n){
      sigma_pos_inv<- (((beta_m[i,]*gamma_indicator_y[i,])%*%t(beta_m[i,]*gamma_indicator_y[i,]))/sigma2w[i])+ 
        diag((num_visit[k]/sigma2k[(i-1)])+(1/sigma2[(i-1)]),num_mediator)
      sigma_pos<- solve(sigma_pos_inv)
      
      t_num<- y[k]-t(d[k,])%*%beta_x[i,]-beta_z[i]*x[k]
      l_vec<- alpha_x[i,,]%*%d[k,]+x[k]*(alpha_z[i,]*gamma_indicator_theta[i,])
      
      sum_m<- apply(mik_qr[(i-1),k,,],2,sum)
      
      mu_pos<- sigma_pos%*%((as.numeric(t_num)*(beta_m[i,]*gamma_indicator_y[i,])/sigma2w[i])+(sum_m/sigma2k[(i-1)])+(l_vec/sigma2[(i-1)]))
      
      mi_qr[i,k,]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
    }
    
    
    #sigma2
    sum_sigma2<- 0
    for (k in 1: num_mediator){
      mu_vec<- d%*%alpha_x[i,k,]+ (alpha_z[i,k]*gamma_indicator_theta[i,k])*x
      sum_sigma2<- sum_sigma2+ sum((mi_qr[i,,k]- mu_vec)^2)
    }
    
    sigma2[i]<- 1/rgamma(1,(a1+(n*num_mediator/2)),(b1+(sum_sigma2/2)))
    
    #paic: paic<- matrix(0,nrow=mcmc_samples,ncol=num_subnetwork)
    
    for (s in 1:num_subnetwork){
      dir_vec_pos[s]<- sum(c_latent[(i-1),,s])+ dir_vec[s]
    }
    
    paic[i,]<- rdirichlet(1, dir_vec_pos)
    
    
    #c_latent<- array(0,dim=c(mcmc_samples,num_region,num_subnetwork)): 
    # cannot do in parallel, must update in time
    correct_c_latent_end<- 0
    c_latent[i,,]<- c_latent[(i-1),,]
    
    for (s in 1:num_region){
      
      for (p in 1: num_subnetwork){
        
        cloglik<- 0
        for (l in 1:num_region){
          
          l_subnetwork<- which(c_latent[i,l,]==1)
          min_l_p<- min(p,l_subnetwork)
          max_l_p<- max(p,l_subnetwork)
          k_cal<- (min_l_p-1)*num_subnetwork+ max_l_p- (min_l_p*(min_l_p-1)/2)
          
          matrix_sum<- dnorm(z[,,s,l],mik_qr[(i-1),,,k_cal],
                             sqrt(sigma2qr_matrix[(i-1),min_l_p,max_l_p]),log = T)
          
          cloglik<- sum(matrix_sum)+cloglik
        }
        
        paic_sp_pos[p]<- cloglik+ c_latent[i,s,p]*log(paic[i,p])
      }
      
      #exp(log(paic_sp_pos)- max(log(paic_sp_pos)))/sum(exp(log(paic_sp_pos)- max(log(paic_sp_pos))))
      paic_sp_pos<- paic_sp_pos- max(paic_sp_pos)
      paic_sp_pos<- exp(paic_sp_pos)/sum(exp(paic_sp_pos))
      c_latent[i,s,]<- t(rmultinom(1, size = 1, prob = paic_sp_pos))
      
    }
    
    for (j in 1:num_region){
      final_class[j]<- which.max(apply(c_latent[(1:i),j,],2,sum))
    }
    
    #mijk: both matrix and vector get updated: can do in parallel, matrix for region, vector for theta
    for (p in 1:num_subnetwork){
      for (q in p:num_subnetwork){
        if (p!=q){
          for (j in 1:n){
            for (s in 1:num_visit[j]){
              sum_mpq<- 0
              num_mpq<- 0
              for (e in which(c_latent[i,,p]==1)){
                for (f in which(c_latent[i,,q]==1)){
                  
                  sum_mpq<- sum_mpq+ z[j,s,e,f]
                  num_mpq<- num_mpq+1
                  
                }
              }
              
              sigma_pos_inv<- (num_mpq/sigma2qr_matrix[(i-1),p,q])+(1/sigma2k[(i-1)])
              sigma_pos<- solve(sigma_pos_inv)
              
              k_cal<- (p-1)*num_subnetwork+ q- (p*(p-1)/2)
              mu_pos<- sigma_pos*((sum_mpq/sigma2qr_matrix[(i-1),p,q])+(mi_qr[i,j,k_cal]/sigma2k[(i-1)]))
              mik_qr[i,j,s,k_cal]<- rnorm(1,mu_pos,sqrt(sigma_pos))
            }
            
          }
        } else {
          for (j in 1:n){
            for (s in 1:num_visit[j]){
              sum_mpq<- 0
              num_mpq<- 0
              for (e in which(c_latent[i,,p]==1)){
                for (f in which(c_latent[i,,q]==1)){
                  if (f>=e){
                    sum_mpq<- sum_mpq+ z[j,s,e,f]
                    num_mpq<- num_mpq+1
                  }
                  
                  
                }
              }
              
              sigma_pos_inv<- (num_mpq/sigma2qr_matrix[(i-1),p,q])+(1/sigma2k[(i-1)])
              sigma_pos<- solve(sigma_pos_inv)
              
              k_cal<- (p-1)*num_subnetwork+ q- (p*(p-1)/2)
              mu_pos<- sigma_pos*((sum_mpq/sigma2qr_matrix[(i-1),p,q])+(mi_qr[i,j,k_cal]/sigma2k[(i-1)]))
              mik_qr[i,j,s,k_cal]<- rnorm(1,mu_pos,sqrt(sigma_pos))
            }
            
          }
          
        }
        
        
      }
      
    }
    
    
    #sigma2k
    
    sum_sigma2k<- 0
    for (j in 1:n){
      for (s in 1:num_visit[j]){
        sum_sigma2k<- sum_sigma2k+ sum((mik_qr[i,j,s,]-mi_qr[i,j,])^2)
      }
    }
    
    sigma2k[i]<- 1/rgamma(1,(ak+(sum(num_visit)/2)),(bk+(sum_sigma2k/2)))
    
    #sigma2(p,q): can do in parallel/ no need to use vector form
    for (p in 1:num_subnetwork){
      for (q in p:num_subnetwork){
        k_cal<- (p-1)*num_subnetwork+ q- (p*(p-1)/2)
        sum_sigm2pq<- 0
        num_sigm2pq<- 0
        if (p!=q){
          for (j in 1:n){
            for (s in 1:num_visit[j]){
              for (e in which(c_latent[i,,p]==1)){
                for (f in which(c_latent[i,,q]==1)){
                  sum_sigm2pq<- sum_sigm2pq+ (z[j,s,e,f]-mik_qr[i,j,s,k_cal])^2
                  num_sigm2pq<- num_sigm2pq+1
                }
              }
            }
          }
          
          
          #sigma2qr_matrix<- array(0,dim=c(mcmc_samples,num_subnetwork,num_subnetwork))
          sigma2qr_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
        } else {
          for (j in 1:n){
            for (s in 1:num_visit[j]){
              for (e in which(c_latent[i,,p]==1)){
                for (f in which(c_latent[i,,q]==1)){
                  if (f>=e){
                    sum_sigm2pq<- sum_sigm2pq+ (z[j,s,e,f]-mik_qr[i,j,s,k_cal])^2
                    num_sigm2pq<- num_sigm2pq+1
                  }
                }
              }
            }
          }
          
          #sigma2qr_matrix<- array(0,dim=c(mcmc_samples,num_subnetwork,num_subnetwork))
          sigma2qr_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
        }
        
      }
      
    }
    
  }
}


