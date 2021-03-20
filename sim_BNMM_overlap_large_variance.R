

############################################### Data Generation ########################################################

path <- "/../../"

###seed to control true parameters the same 
args<- c(50002)
set.seed(50002)

#library
library(MASS)
library(mvtnorm)
library(gdata)
library(gtools)

###global settings
n<- 50
num_visit_maximum<- 6
num_visit<- rep(0,n)
for (i in 1:n){
  num_visit[i]<- num_visit_maximum
}

num_region<- 100
num_subnetwork<- 10 #assume we know this
num_mediator<- ((num_subnetwork+1)*num_subnetwork)/2

#gamma_indicator
p_true<- 0.15
gamma_indicator_true_theta<- rep(0,times=num_mediator)
gamma_indicator_true_theta[seq(from=5,to=40,by=5)]<- 1
num_true_mediator_theta<- sum(gamma_indicator_true_theta==1)
index_gamma_indicator_theta_true<- which(gamma_indicator_true_theta==1)

gamma_indicator_true_y<- rep(0,times=num_mediator)
gamma_indicator_true_y[c(seq(from=10,to=30,by=5),40,50,55)]<- 1
num_true_mediator_y<- sum(gamma_indicator_true_y==1)
index_gamma_indicator_y_true<- which(gamma_indicator_true_y==1)

gamma_indicator_true<- rep(0,times=num_mediator)#  true
gamma_indicator_true[c(seq(from=10,to=30,by=5),40)]<- 1
num_true_mediator<- sum(gamma_indicator_true==1)
index_gamma_indicator_true<- which(gamma_indicator_true==1)

#true parameters
alpha0_true<- 1.5
alpha_theta_true<- rep(2,times=num_true_mediator_y)
beta_theta_true<- seq(from=1.5, to=2.5, length.out=num_true_mediator_theta)

#direct effect
sum(alpha_theta_true[1:6]*beta_theta_true[c(2:6,8)])

#latent class label
c_latent_true<- matrix(0,nrow=num_region,num_subnetwork)
dir_true<- rep(3,times=num_subnetwork)
paic_true<- rdirichlet(1, dir_true)
for (i in 1:num_region){
  c_latent_true[i,]<- t(rmultinom(1, size = 1, prob = paic_true))
}


###input structure: observed
#vector y: outcome, n*1
#vector x: treatment, n*1
#array z (after fisher transformation from the original correlated table): dim represents subject, visit, a picture respectively 

### seed to simulate x,y,z
set.seed(as.numeric(args[1]))

###design vectors/matrices
#x: treatment
#x<- rbinom(n,1,0.7)
x<- rnorm(n)


#theta_true
theta_matrix_true<- matrix(0,nrow=n,ncol=num_mediator)
sigma2_true<- 0.5^2

for (k in 1:num_mediator){
  theta_matrix_true[,k]<- rnorm(n,0.3,sqrt(sigma2_true))
}

mean_theta<- matrix(0,nrow=n,ncol=num_true_mediator_theta)
for (k in 1:num_true_mediator_theta){
  mean_theta[,k]<- rep(0.3,times=n)+ beta_theta_true[k]*x
  theta_matrix_true[,index_gamma_indicator_theta_true[k]]<- rmvnorm(1,mean_theta[,k],diag(sigma2_true,n))
  
}

#y
y<- rep(0,times=n)
sigma2w_true<- 1^2
for (i in 1:n){
  mean_y<-  t(theta_matrix_true[i,index_gamma_indicator_y_true])%*%alpha_theta_true+ x[i]*alpha0_true
  y[i]<- rnorm(1,mean_y,sqrt(sigma2w_true))
}

#mij
m_connectivity_strength_vector_true<- array(0,dim=c(n,num_visit_maximum,num_mediator))
sigma2k_true<- 0.5^2
for (i in 1:n){
  for (j in 1:num_visit[i]){
    for (k in 1:num_mediator){
      m_connectivity_strength_vector_true[i,j,k]<- rnorm(n=1,mean=theta_matrix_true[i,k],sqrt(sigma2k_true))
    }
  }
}

#generate z
sigma2pq_true<- 0.5^2
z<- array(0,dim=c(n,num_visit_maximum,num_region,num_region))
for (i in 1:n){
  for (j in 1:num_visit[i]){
    for (s in 1:num_region){
      for (t in s:num_region){
        s_subnetwork<- which(c_latent_true[s,]==1)
        t_subnetwork<- which(c_latent_true[t,]==1)
        min_l_p<- min(s_subnetwork,t_subnetwork)
        max_l_p<- max(s_subnetwork,t_subnetwork)
        k_cal<- (min_l_p-1)*num_subnetwork+ max_l_p- (min_l_p*(min_l_p-1)/2)
        z[i,j,s,t]<- rnorm(1,m_connectivity_strength_vector_true[i,j,k_cal],sqrt(sigma2pq_true))
        z[i,j,t,s]<- z[i,j,s,t]
      }
    }
  }
}

############################################### Model Fitting ########################################################

#initialized values 
mcmc_samples_1<- 1000
gamma_indicator_y<- matrix(0,nrow=mcmc_samples_1,ncol=num_mediator)
gamma_indicator_theta<- matrix(0,nrow=mcmc_samples_1,ncol=num_mediator)
gamma_indicator_y[1,]<- gamma_indicator_true
gamma_indicator_theta[1,]<- gamma_indicator_true
p_gamma_y<- 0.5
p_gamma_theta<- 0.5

gamma_temp_y<- gamma_indicator_y[1,]
gamma_temp_theta<- gamma_indicator_theta[1,]

gamma_indicator_temp_decision_y<- rep(0,times=num_mediator)
gamma_indicator_temp_decision_theta<- rep(0,times=num_mediator)

alpha_theta<- matrix(0,nrow=mcmc_samples_1,ncol=num_mediator)
alpha_theta[1,index_gamma_indicator_y_true]<- alpha_theta_true

beta_theta<- matrix(0,nrow=mcmc_samples_1,ncol=num_mediator)
beta_theta[1,index_gamma_indicator_theta_true]<- beta_theta_true

#for real-time monitor
alpha_theta_mean<- rep(0,times=num_mediator)
beta_theta_mean<- rep(0,times=num_mediator)

alpha0<- rep(0,times=mcmc_samples_1)
alpha0[1]<- alpha0_true
sigma20<- 10^2

sigma2a<- rep(0,times=mcmc_samples_1)

sigma2a[1]<- 1^2
a3<- 0.1
b3<- 0.1

sigma2b<- rep(0,times=mcmc_samples_1)

sigma2b[1]<- 1^2
a4<- 0.1
b4<- 0.1


sigma2w<- rep(0,times=mcmc_samples_1)
sigma2w[1]<- 1^2
a2<- 2
b2<- 2

theta_matrix<- array(0,dim=c(mcmc_samples_1,n,num_mediator))
theta_matrix[1,,]<- theta_matrix_true

m_connectivity_strength_vector<- array(0,dim=c(mcmc_samples_1,n,num_visit_maximum,num_mediator))
m_connectivity_strength_vector[1,,,]<- m_connectivity_strength_vector_true

c_latent<- array(0,dim=c(mcmc_samples_1,num_region,num_subnetwork))
for (i in 1:mcmc_samples_1){
  c_latent[i,,]<- c_latent_true
}

sigma_connectivity_strength_matrix<- array(0,dim=c(mcmc_samples_1,num_subnetwork,num_subnetwork))
sigma_connectivity_strength_matrix[1,,]<- matrix(sigma2pq_true,nrow=num_subnetwork,ncol=num_subnetwork)
a0<- 1
b0<- 1

sigma2<- rep(0,times=mcmc_samples_1)
sigma2[1]<- sigma2_true
a1<- 1
b1<- 1


sigma2k<- rep(0,times=mcmc_samples_1)
sigma2k[1]<- sigma2k_true
ak<- 1
bk<- 1

d<- matrix(1,nrow=n,ncol=1)
pd<- dim(d)[2]

w_matrix<- array(0,dim=c(mcmc_samples_1,num_mediator,pd))
sigma2m<- 10^2

for (i in 2:mcmc_samples_1) {
  
  #gamma_for_y
  
  tp_num<- 0
  fp_num<- 0
  
  var_y<- diag(sigma2w[(i-1)],n)
  
  for (k in 1:num_mediator){
    #loglik for 1
    gamma_temp_y[k]<- 1
    mean_y_1<- theta_matrix[(i-1),,]%*%(alpha_theta[(i-1),]*gamma_temp_y)+alpha0[(i-1)]*x
    val_1<- dmvnorm(t(y),mean_y_1,var_y,log = T) + log(p_gamma_y)
    
    gamma_temp_y[k]<- 0
    mean_y_0<- theta_matrix[(i-1),,]%*%(alpha_theta[(i-1),]*gamma_temp_y)+alpha0[(i-1)]*x
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
  
  # gamma_for_theta
  
  var_thetak<- diag(sigma2[(i-1)],n)
  
  for (k in 1:num_mediator){
    #loglik for 1
    gamma_temp_theta[k]<- 1
    mean_thetak_1<- d%*%w_matrix[(i-1),k,]+ (beta_theta[(i-1),k]*gamma_temp_theta[k])*x
    val_1<- dmvnorm(t(theta_matrix[(i-1),,k]),mean_thetak_1,var_thetak,log=T)+log(p_gamma_theta)
    
    #loglik for 0
    gamma_temp_theta[k]<- 0
    mean_thetak_0<- d%*%w_matrix[(i-1),k,]+ (beta_theta[(i-1),k]*gamma_temp_theta[k])*x
    val_0<- dmvnorm(t(theta_matrix[(i-1),,k]),mean_thetak_0,var_thetak,log=T)+log(1-p_gamma_theta)
    
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
  
  #alpha_theta
  
  gamma_chosen_y<- which(gamma_temp_y==1)
  gamma_not_chosen_y<- which(gamma_temp_y==0)
  
  if (length(gamma_not_chosen_y)>0){
    alpha_theta[i,gamma_not_chosen_y]<- mvrnorm(n=1,mu=rep(0,times=length(gamma_not_chosen_y)),Sigma=diag(sigma2a[(i-1)],length(gamma_not_chosen_y)))
  }
  
  if (length(gamma_chosen_y)>0){
    t_vec<- y- alpha0[(i-1)]*x
    sigma_pos_inv<- (sigma2a[(i-1)]*t(theta_matrix[(i-1),,gamma_chosen_y])%*%theta_matrix[(i-1),,gamma_chosen_y]+diag(sigma2w[(i-1)],length(gamma_chosen_y)))/(sigma2w[(i-1)]*sigma2a[(i-1)])
    sigma_pos<- solve(sigma_pos_inv)
    mu_pos<- (sigma_pos%*%t(theta_matrix[(i-1),,gamma_chosen_y])%*%t_vec)/sigma2w[(i-1)]
    
    alpha_theta[i,gamma_chosen_y]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
  }
  
  
  #beta_theta
  
  gamma_chosen_theta<- which(gamma_temp_theta==1)
  gamma_not_chosen_theta<- which(gamma_temp_theta==0)
  
  if (length(gamma_not_chosen_theta)>0){
    beta_theta[i,gamma_not_chosen_theta]<- mvrnorm(n=1,mu=rep(0,times=length(gamma_not_chosen_theta)),Sigma=diag(sigma2b[(i-1)],length(gamma_not_chosen_theta)))
  }
  
  if (length(gamma_chosen_theta>0)){
    sigma_pos_inv<- diag((sigma2b[(i-1)]*sum(x^2)+sigma2[(i-1)])/(sigma2[(i-1)]*sigma2b[(i-1)]),length(gamma_chosen_theta))
    sigma_pos<- solve(sigma_pos_inv)
    
    sum_beta_theta<- 0
    for (j in 1:n){
      t_vec<- theta_matrix[(i-1),j,gamma_chosen_theta]- w_matrix[(i-1),gamma_chosen_theta,]*d[j,]
      sum_beta_theta<- sum_beta_theta+ x[j]*t_vec
    }
    
    mu_pos<- (sigma_pos%*%sum_beta_theta)/sigma2[(i-1)]
    
    beta_theta[i,gamma_chosen_theta]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
  }
  
  
  #alpha0
  sigma_pos_inv<- (sum(x^2)/sigma2w[(i-1)]+ (1/sigma20))
  sigma_pos<- 1/sigma_pos_inv
  
  t_vec<- y- theta_matrix[(i-1),,]%*%(alpha_theta[i,]*gamma_indicator_y[i,])
  sum_alpha0<- sum(x*t_vec)
  
  mu_pos_alpha0<- (sigma_pos*sum_alpha0)/sigma2w[(i-1)]
  
  alpha0[i]<- rnorm(1,mean=mu_pos_alpha0,sd=sqrt(sigma_pos))
  
  #sigma2a
  
  sigma2a[i]<- 1/rgamma(1,(a3+ (num_mediator/2)),(b3+ (sum(alpha_theta[i,]^2)/2)))
  
  #sigma2b
  
  sigma2b[i]<- 1/rgamma(1,(a4+ (num_mediator/2)),(b4+ (sum(beta_theta[i,]^2)/2)))
  
  #sigma2w
  
  t_vec<- y-theta_matrix[(i-1),,gamma_chosen_y]%*%as.matrix(alpha_theta[i,gamma_chosen_y])-alpha0[i]*x
  sigma2w[i]<- 1/rgamma(1,(a2+ (n/2)),(b2+ (sum(t_vec^2)/2)))
  
  
  #w_matrix
  
  for (k in 1:num_mediator){
    t_vec<- theta_matrix[(i-1),,k]- ((beta_theta[i,k]*gamma_indicator_theta[i,k])*x)
    sigma_pos_inv<- (sigma2m*t(d)%*%d+diag(sigma2[(i-1)],pd))/(sigma2[(i-1)]*sigma2m)
    sigma_pos<- solve(sigma_pos_inv)
    mu_pos<- (sigma_pos%*%t(d)%*%t_vec)/sigma2[(i-1)]
    
    w_matrix[i,k,]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
  }
  
  
  #theta_matrix: update for each individual theta_i
  for (k in 1:n){
    sigma_pos_inv<- (((alpha_theta[i,]*gamma_indicator_y[i,])%*%t(alpha_theta[i,]*gamma_indicator_y[i,]))/sigma2w[i])+ diag((num_visit[k]/sigma2k[(i-1)])+(1/sigma2[(i-1)]),num_mediator)
    sigma_pos<- solve(sigma_pos_inv)
    
    t_num<- y[k]-alpha0[i]*x[k]
    l_vec<- w_matrix[i,,]*d[k,]+x[k]*(beta_theta[i,]*gamma_indicator_theta[i,])
    
    sum_m<- apply(m_connectivity_strength_vector[(i-1),k,,],2,sum)
    
    mu_pos<- sigma_pos%*%((as.numeric(t_num)*(alpha_theta[i,]*gamma_indicator_y[i,])/sigma2w[i])+(sum_m/sigma2k[(i-1)])+(l_vec/sigma2[(i-1)]))
    
    theta_matrix[i,k,]<- mvrnorm(n=1,mu=mu_pos,Sigma=sigma_pos)
  }
  
  
  #sigma2
  
  sum_sigma2<- 0
  for (k in 1: num_mediator){
    mu_vec<- d%*%w_matrix[i,k,]+ (beta_theta[i,k]*gamma_indicator_theta[i,k])*x
    sum_sigma2<- sum_sigma2+ sum((theta_matrix[i,,k]- mu_vec)^2)
  }
  
  sigma2[i]<- 1/rgamma(1,(a1+(n*num_mediator/2)),(b1+(sum_sigma2/2)))
  
 
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
            
            sigma_pos_inv<- (num_mpq/sigma_connectivity_strength_matrix[(i-1),p,q])+(1/sigma2k[(i-1)])
            sigma_pos<- solve(sigma_pos_inv)
            
            k_cal<- (p-1)*num_subnetwork+ q- (p*(p-1)/2)
            mu_pos<- sigma_pos*((sum_mpq/sigma_connectivity_strength_matrix[(i-1),p,q])+(theta_matrix[i,j,k_cal]/sigma2k[(i-1)]))
            m_connectivity_strength_vector[i,j,s,k_cal]<- rnorm(1,mu_pos,sqrt(sigma_pos))
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
            
            sigma_pos_inv<- (num_mpq/sigma_connectivity_strength_matrix[(i-1),p,q])+(1/sigma2k[(i-1)])
            sigma_pos<- solve(sigma_pos_inv)
            
            k_cal<- (p-1)*num_subnetwork+ q- (p*(p-1)/2)
            mu_pos<- sigma_pos*((sum_mpq/sigma_connectivity_strength_matrix[(i-1),p,q])+(theta_matrix[i,j,k_cal]/sigma2k[(i-1)]))
            m_connectivity_strength_vector[i,j,s,k_cal]<- rnorm(1,mu_pos,sqrt(sigma_pos))
          }
          
        }
        
      }
      
    }
    
  }
  
  
  #sigma2k
  
  sum_sigma2k<- 0
  for (j in 1:n){
    for (s in 1:num_visit[j]){
      sum_sigma2k<- sum_sigma2k+ sum((m_connectivity_strength_vector[i,j,s,]-theta_matrix[i,j,])^2)
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
                sum_sigm2pq<- sum_sigm2pq+ (z[j,s,e,f]-m_connectivity_strength_vector[i,j,s,k_cal])^2
                num_sigm2pq<- num_sigm2pq+1
              }
            }
          }
        }
        
        
        #sigma_connectivity_strength_matrix<- array(0,dim=c(mcmc_samples_1,num_subnetwork,num_subnetwork))
        sigma_connectivity_strength_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
      } else {
        for (j in 1:n){
          for (s in 1:num_visit[j]){
            for (e in which(c_latent[i,,p]==1)){
              for (f in which(c_latent[i,,q]==1)){
                if (f>=e){
                  sum_sigm2pq<- sum_sigm2pq+ (z[j,s,e,f]-m_connectivity_strength_vector[i,j,s,k_cal])^2
                  num_sigm2pq<- num_sigm2pq+1
                }
              }
            }
          }
        }
        
        
        #sigma_connectivity_strength_matrix<- array(0,dim=c(mcmc_samples_1,num_subnetwork,num_subnetwork))
        sigma_connectivity_strength_matrix[i,p,q]<- 1/rgamma(1,(a0+(num_sigm2pq/2)),(b0+(sum_sigm2pq/2)))
      }
      
    }
    
  }
  
  print(i)
  print(c("alpha0_true:",alpha0_true))
  print(c("alpha0_mu:",mu_pos_alpha0))
  print(c("alpha0_mean:",mean(alpha0[1:i])))
  print(c("tot_num_of_true_significant_mediator:",num_true_mediator))
  print(c("tot_num_of_significant_mediator:",sum(gamma_indicator_temp_decision_y*gamma_indicator_temp_decision_theta)))
  print(c("tot_prop_of_correct_mediator%:",(num_correct_mediator/num_mediator)*100))
  print(c("true_positive:",tp_num))
  print(c("false_positive:",fp_num))
  
}


#save the results as .rdata
save.image(file=paste0(path, 
                       "mediation_simulation_n_full_version_overlap_assuming_c_true_large_variance_",args[1],".rdata"))


