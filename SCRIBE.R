####################################
### Metropolis hasting algorithm ###
####################################

MH_bg <- function(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg, Lambda, Pi, Nu, Sigma, step, burn, K){
  bg_sum <- matrix(0, ncol=N_batch, nrow=Gene_N)
  bgsq_sum <- matrix(0, ncol=N_batch, nrow=Gene_N)
  ita <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Zexpbg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bg_propose <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Phi <- pnorm(Pi)
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg[,batches]), FUN='*')
  }
  log_int_z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bg[Y_zero])-1))
  log_int_z[!Y_zero] <- Y[!Y_zero]*log(Lambda_bg[!Y_zero])-Lambda_bg[!Y_zero]
  log_dense <- matrix(0, ncol=N_batch, nrow=Gene_N)
  log_dense_propose <- matrix(0, ncol=N_batch, nrow=Gene_N)
  for(batches in 1:N_batch){
    log_dense[, batches] <- -(bg[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
  }
  for(i in 1:burn){
    bg_propose <- bg + matrix(rnorm(N_batch*Gene_N, mean=0, sd=step), ncol=N_batch)
    for(batches in 1:N_batch){
      Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_propose[,batches]), FUN='*')
    }
    log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bg[Y_zero])-1))
    log_int_z[!Y_zero] <- Y[!Y_zero]*log(Lambda_bg[!Y_zero])-Lambda_bg[!Y_zero]
    for(batches in 1:N_batch){
      log_dense_propose[, batches] <- -(bg_propose[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
    }
    log_ratio <- matrix(sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0)), ncol=N_batch)
    ratio <- exp(log_ratio)
    hidden <- matrix(runif(N_batch*Gene_N), ncol=N_batch)
    accept <- ratio>hidden
    bg[accept] <- bg_propose[accept]
    log_dense[accept] <- log_dense_propose[accept]
  }
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg[,batches]), FUN='*')
  }
  for(i in 1:K){
    bg_propose <- bg + matrix(rnorm(N_batch*Gene_N, mean=0, sd=step), ncol=N_batch)
    for(batches in 1:N_batch){
      Lambda_bg_propose[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_propose[,batches]), FUN='*')
    }
    log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bg_propose[Y_zero])-1))
    log_int_z[!Y_zero] <- Y[!Y_zero]*log(Lambda_bg_propose[!Y_zero])-Lambda_bg_propose[!Y_zero]
    for(batches in 1:N_batch){
      log_dense_propose[, batches] <- -(bg_propose[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
    }
    log_ratio <- matrix(sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0)), ncol=N_batch)
    ratio <- exp(log_ratio)
    hidden <- matrix(runif(N_batch*Gene_N), ncol=N_batch)
    accept <- ratio>hidden
    bg[accept] <- bg_propose[accept]
    for(batches in 1:N_batch){
      Lambda_bg[accept[,batches], bat_ind == batch_name[batches]] <- Lambda_bg_propose[accept[,batches],bat_ind == batch_name[batches]]
    }
    log_dense[accept] <- log_dense_propose[accept]
    bg_sum <- bg_sum + bg
    bgsq_sum <- bgsq_sum + bg^2
    expbg <- exp(bg)
    explambda <- exp(-Lambda_bg)
    z_temp <- Phi*explambda/(1+Phi*(explambda-1))
    z_temp[!Y_zero] <- 1
    Z <- Z + z_temp
    ita <- ita + Pi + z_temp*exp(-Pi^2/2)/(sqrt(2*pi)*Phi) - (1-z_temp)*exp(-Pi^2/2)/(sqrt(2*pi)*(1-Phi)) 
    for(batches in 1:N_batch){
      Zexpbg[,bat_ind == batch_name[batches]] <- Zexpbg[,bat_ind == batch_name[batches]] + sweep(z_temp[, bat_ind == batch_name[batches]], 1, expbg[,batches], FUN='*')
    }
  }
  return(list(bg_sum=bg_sum, bgsq_sum=bgsq_sum, ita=ita, Zexpbg=Zexpbg, Z=Z))
}

#####################
### main function ###
#####################

fitSCRIBE <- function(Y, bat_ind, bio_ind, burn=50, burn_start=500, K=500, max_iter=100){

  Y <- as.matrix(Y)
  Y_zero <- Y <= 0
  Cell_N <- ncol(Y)
  Gene_N <- nrow(Y)
  gene_len <- rowMeans(Y)/mean(Y)
  group_name <- unique(bio_ind)
  N_group <- length(group_name)
  batch_name = unique(bat_ind)
  N_batch <- length(batch_name)
  
  tot_read <- rep(0, Cell_N)
  for(batches in 1:N_batch){
    tot_read[bat_ind==batch_name[batches]] <- colMeans(Y[,bat_ind==batch_name[batches]])/mean(Y[,bat_ind==batch_name[batches]])
  }
  
  #start point
  Gamma <- rep(0, N_batch)
  Alpha <- rep(0, N_batch)
  Beta <- rep(0, N_batch)
  Nu <- rep(0, N_batch)
  Sigma <- rep(0.5, N_batch)
  Mu <- matrix(0,nrow=Gene_N,ncol=N_group)
  for(groups in 1:N_group){
    Mu[,groups] <- rowMeans(sweep(Y[,bio_ind==group_name[groups]], 2, 1/tot_read[bio_ind==group_name[groups]], FUN='*'))/(gene_len)
  }
  bg <- matrix(0, ncol=N_batch, nrow=Gene_N)
  for(batches in 1:N_batch){
    bg[,batches] <- rnorm(Gene_N, mean=Nu[batches], sd=Sigma[batches])
  }
  
  Group_Matrix <- matrix(0, ncol=Cell_N, nrow=N_group)
  for(groups in 1:N_group){
    Group_Matrix[groups,]<-bio_ind==group_name[groups]
  }
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  Pi <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi[, bat_ind == batch_name[batches]], 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]]), FUN="+"), 1, Beta[batches]*log(gene_len), FUN="+")
  }
  
  # E step
  step = 0.1
  MH <- MH_bg(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg, Lambda, Pi, Nu, Sigma, step, burn_start, K)
  bg_av <- MH$bg_sum/K
  bgsq_av <- MH$bgsq_sum/K
  ita <- MH$ita/K
  Zexpbg <- MH$Zexpbg
  
  #M step
  Gamma1 <- Gamma
  Alpha1 <- Alpha
  Beta1 <- Beta
  Nu1 <- Nu
  Sigma1 <- Sigma
  Mu1 <- Mu
  for(groups in 1:N_group){
    Mu[,groups] <- K*rowSums(Y[,bio_ind==group_name[groups]])/(gene_len*rowSums(sweep(Zexpbg[,bio_ind==group_name[groups]], 2, tot_read[bio_ind==group_name[groups]], FUN="*")))
  }
  
  for(batches in 1:N_batch){
    batch_ita <- ita[,bat_ind == batch_name[batches]]
    response <- as.vector(batch_ita)
    logrc <- rep(log(tot_read[bat_ind == batch_name[batches]]), each=Gene_N)
    loglg <- rep(log(gene_len), sum(bat_ind == batch_name[batches]))
    LM_Results <- lm(response~logrc+loglg)
    Gamma[batches] <- LM_Results$coefficients[[1]]
    Alpha[batches] <- LM_Results$coefficients[[2]]
    Beta[batches] <- LM_Results$coefficients[[3]]
  }
  Nu <- colMeans(bg_av)
  Sigma <- sqrt(colMeans(bgsq_av)-Nu^2)
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  Pi <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi[, bat_ind == batch_name[batches]], 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]]), FUN="+"), 1, Beta[batches]*log(gene_len), FUN="+")
  }
  stepthreshold <- step/2
  stepsq = ifelse(as.vector(bgsq_av-bg_av^2)>stepthreshold^2, as.vector(bgsq_av-bg_av^2), stepthreshold^2)
  step <- sqrt(stepsq)
  j <- 1
  while(j<=max_iter & max(abs((Alpha1-Alpha)/Alpha), abs((Beta1-Beta)/Beta), abs((Gamma1-Gamma)/Gamma), abs((Nu1-Nu)/Nu))>0.001){
    # E step
    MH <- MH_bg(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg_av, Lambda, Pi, Nu, Sigma, step, burn, K)
    bg_av <- MH$bg_sum/K
    bgsq_av <- MH$bgsq_sum/K
    ita <- MH$ita/K
    Zexpbg <- MH$Zexpbg
    
    #M step
    Gamma1 <- Gamma
    Alpha1 <- Alpha
    Beta1 <- Beta
    Nu1 <- Nu
    Sigma1 <- Sigma
    Mu1 <- Mu
    
    for(groups in 1:N_group){
      Mu[,groups] <- K*rowSums(Y[,bio_ind==group_name[groups]])/(gene_len*rowSums(sweep(Zexpbg[,bio_ind==group_name[groups]], 2, tot_read[bio_ind==group_name[groups]], FUN="*")))
    }
    
    for(batches in 1:N_batch){
      batch_ita <- ita[,bat_ind == batch_name[batches]]
      response <- as.vector(batch_ita)
      logrc <- rep(log(tot_read[bat_ind == batch_name[batches]]), each=Gene_N)
      loglg <- rep(log(gene_len), sum(bat_ind == batch_name[batches]))
      LM_Results <- lm(response~logrc+loglg)
      Gamma[batches] <- LM_Results$coefficients[[1]]
      Alpha[batches] <- LM_Results$coefficients[[2]]
      Beta[batches] <- LM_Results$coefficients[[3]]
    }
    
    Nu <- colMeans(bg_av)
    Sigma <- sqrt(colMeans(bgsq_av)-Nu^2)
    Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
    Pi <- matrix(0, ncol=Cell_N, nrow=Gene_N)
    for(batches in 1:N_batch){
      Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi[, bat_ind == batch_name[batches]], 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]]), FUN="+"), 1, Beta[batches]*log(gene_len), FUN="+")
    }
    stepthreshold <- stepthreshold/2
    stepsq = ifelse(as.vector(bgsq_av-bg_av^2)>stepthreshold^2, as.vector(bgsq_av-bg_av^2), stepthreshold^2)
    step <- sqrt(stepsq)
    cat("iteration", j, '\n')
    cat("Alpha=", Alpha, '\n')
    cat("Beta=", Beta, '\n')
    cat("Gamma=", Gamma, '\n')
    cat("Nu=", Nu, '\n')
    cat("Sigma=", Sigma, '\n')
    j <- j + 1
  }
  Lambda_bg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_av[,batches]), FUN='*')
  }
  meanNu <- mean(Nu)
  bg_av = bg_av - meanNu
  Nu = Nu - meanNu
  Mu = Mu * exp(meanNu)
  Z <- MH$Z/K
  Lambda <- sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix
  # dropout imputation
  cat("start imputing ... ...")
  Y_impute <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Y_impute[Lambda_bg==0] <- 0
  Y_impute[Lambda_bg!=0] <- Y[Lambda_bg!=0] * sqrt(Lambda[Lambda_bg!=0]/Lambda_bg[Lambda_bg!=0]) - sqrt(Lambda[Lambda_bg!=0]*Lambda_bg[Lambda_bg!=0]) + Lambda[Lambda_bg!=0]
  Y_impute[Y_impute<0] <- 0
  Y_impute[Y_zero] <- Z[Y_zero] * Y_impute[Y_zero] + (1 - Z[Y_zero]) * rpois(sum(Y_zero), Lambda[Y_zero])
  return(list(Mu=Mu, Alpha=Alpha, Beta=Beta, Gamma=Gamma, Nu=Nu, Sigma=Sigma, 
              Y_impute = Y_impute, bg=bg_av, Groups=group_name, Batches=batch_name))
}