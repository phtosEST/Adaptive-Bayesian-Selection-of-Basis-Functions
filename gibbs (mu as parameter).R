library(fda)
library(invgamma)
library(dqrng)
library(mvnfast)


# Simulated Data
Xt = seq(0,1,length=100)
basisBspline_Simulated_Data = create.bspline.basis(range(Xt),norder=4,nbasis=10)
B_simulated_data = getbasismatrix(Xt, basisBspline_Simulated_Data, nderiv=0)
coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0)
set.seed(1234)
argvals = lapply(1:50,function(s) Xt)
raw_data =lapply(1:50,function(s) as.numeric(B_simulated_data%*%coef)+rnorm(length(argvals[[s]]),0,0.1))
y = raw_data




for(Kind in c(5:15,20,25,30)){
  Nint=10000
  K = Kind 
  m = 5
  n = unlist(lapply(argvals,length))
  
  # B-SPLINE MATRIX 
  # each element of the list is a matrix n_i by K
  B = lapply(argvals,function(s) getbasismatrix(Xt,create.bspline.basis(range(Xt),norder=4,nbasis=K), nderiv=0))
  str(B)
  
  
  # Hyperparameters
  lambda1 = 0
  lambda2 = 0
  delta1 = 0
  delta2 = 0
  
  
  # Parameters - initial value of the chain
  tau2 = 1#5
  sigma2 = 1#5
  beta_names = as.vector(outer(paste0("beta_",1:K,"_"),1:m,paste0))
  beta=rep(-1,m*K)#rep(1,m*K)
  theta_names=as.vector(outer(paste0("theta_",1:K,"_"),1:m,paste0))
  theta=rep(1/5,m*K)#rep(4/5,m*K)
  mu_names=as.vector(outer(paste0("mu_",1:K,"_"),1:m,paste0))
  mu=rep(1/5,m*K)#rep(4/5,m*K)
  
  
  # Latent variable - initial value of the chain
  Z_names = as.vector(outer(paste0("Z_",1:K,"_"),1:m,paste0))
  Z = vector()
  set.seed(1234)
  sapply(1:(m*K),function(ki){
    Z[ki] <<- sample(c(0,1),1)
  })
  # run the next line only for chain2 
  # Z=abs(Z-1)
  
  vec_par = c(beta,theta,mu,sigma2,tau2)
  mat = matrix(NA,ncol = length(vec_par), nrow = Nint)
  # matrix with chains in columns 
  mat[1,] = vec_par
  colnames(mat) = c(beta_names,theta_names,mu_names,"sigma2","tau2")
  matZ = matrix(NA,ncol = length(Z), nrow = Nint)
  matZ[1,] = Z
  colnames(matZ) = Z_names
  
  
  set.seed(Sys.time())
  start=Sys.time()
  for(iteration in 2:Nint){
    print(iteration)
    
    beta_old = mat[iteration-1,grep("^beta_[0-9]+_[0-9]+$",colnames(mat))]
    
    s1 = sum(n[1:m])/2 +((m*K)/2)+ delta1
    save_vals = sum(sapply(1:m,function(i){
      beta_i_old = mat[iteration-1,grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mat))]
      Z_i_old = matZ[iteration-1,grep(paste0("^Z_[0-9]+_",i,"$"),colnames(matZ))]
      value = sum(sapply(1:n[i], function(j){(y[[i]][j]-sum(Z_i_old*beta_i_old*B[[i]][j,]))^2}))
      return(value)
    }))
    s2 = (save_vals+(sum(beta_old^2)/mat[iteration-1,"tau2"])+2*delta2)/2
    sigma2_new = invgamma::rinvgamma(1,shape=s1,rate=s2)
    mat[iteration,"sigma2"] = sigma2_new
    
    
    
    
    t1 = ((m*K)/2)+lambda1
    t2 = ((sum(beta_old^2)/sigma2_new)+2*lambda2)/2
    tau2_new = invgamma::rinvgamma(1,shape=t1,rate=t2)
    mat[iteration,"tau2"] = tau2_new
    
    
    sapply(1:m, function(i){
      
      beta_i_old = mat[iteration-1,grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mat))]
      Z_i_mixed = matZ[iteration,grep(paste0("^Z_[0-9]+_",i,"$"),colnames(matZ))]
      resume2i = matZ[iteration-1,grep(paste0("^Z_[0-9]+_",i,"$"),colnames(matZ))]
      ind_zi = is.na(Z_i_mixed)
      Z_i_mixed[ind_zi] = resume2i[ind_zi]
      
      sapply(1:K, function(k){
        
        Z_i_mixed_0 = Z_i_mixed
        Z_i_mixed_1 = Z_i_mixed
        
        thetaki_old = mat[iteration-1,grep(paste0("^theta_",k,"_",i,"$"),colnames(mat))]
        
        p_mu=thetaki_old
        muki_new=1
        while(muki_new==0|muki_new>=0.6){
          muki_new = ifelse(p_mu==0.5,dqrunif(1)
                            ,log(((2*p_mu-1)*dqrunif(1)-p_mu+1)/(1-p_mu))/log(p_mu/(1-p_mu)))
        }
        mat[iteration,grep(paste0("^mu_",k,"_",i,"$"),colnames(mat))] <<- muki_new
        
        Z_i_mixed_0[k]=0
        Z_i_mixed_1[k]=1
        save_Zvals = sum(sapply(1:n[i],function(j){
          (y[[i]][j]-sum(beta_i_old*Z_i_mixed_1*B[[i]][j,]))^2 - (y[[i]][j]-sum(beta_i_old*Z_i_mixed_0*B[[i]][j,]))^2
        }))
        prob_posteriori_Zki = ifelse(thetaki_old==1,1,ifelse(thetaki_old==0,0,thetaki_old/((1-thetaki_old)*exp(save_Zvals/(2*sigma2_new))+thetaki_old)))
        zki_new = ifelse(dqrunif(1)<prob_posteriori_Zki,1,0)
        matZ[iteration,grep(paste0("^Z_",k,"_",i,"$"),colnames(matZ))] <<- zki_new 
        Z_i_mixed[grep(paste0("^Z_",k,"_",i,"$"),names(Z_i_mixed))] <<- zki_new
        thetaki_new = rbeta(1,shape1=muki_new+zki_new,shape2=(1-muki_new)-zki_new+1)
        mat[iteration,grep(paste0("^theta_",k,"_",i,"$"),colnames(mat))] <<- thetaki_new 
      })
    })
    
    
    
    sapply(1:m, function(i){
      Z_i_new = matZ[iteration,grep(paste0("^Z_[0-9]+_",i,"$"),colnames(matZ))]
      G_values = lapply(1:n[i], function(j){
        G_ij = Z_i_new*B[[i]][j,]
        return(list(matrix_G=G_ij%*%t(G_ij),vector_G=G_ij*y[[i]][j]))
      })
      inv_D_i = solve(diag((1/tau2_new),nrow=K)  + Reduce('+',lapply(G_values, '[[', 1)))
      meanMN = inv_D_i%*%Reduce('+',lapply(G_values, '[[', 2))
      varMN = sigma2_new*inv_D_i
      
      mat[iteration,grep(paste0("^beta_[0-9]+_",i,"$"),colnames(mat))] <<- rmvn(1, meanMN,varMN, ncores = 8)
      
    })
    
    
  }
  Sys.time()-start
  
  chains =  list(mat,matZ)
  # save(chains,file = paste0(".../B-spline/Sigma_1por10/mat_BsplineAuto_m",m,"_K",K,"(chain1).Rdata"))
}


########################################################################################
########################################################################################
