# asymptotics n_E to infty
logpost_offset = function(beta, y, X, prior_sd, offset){
  # browser()
  eta = X%*%beta + offset
  -(sum(y*eta) - sum(log(1+exp(eta))) - sum(beta^2)/2/prior_sd^2 - (log(prior_sd) + log(2*pi)/2)*length(beta))
}
logpost_gr_offset = function(beta, y, X, prior_sd, offset){
  eta = c(X%*%beta) + offset
  -(t(X)%*%y - colSums(1/(1+exp(-eta))*X)) + beta/prior_sd^2
}

marg_lik_laplace_map_offset = function(y, X, prior_sd, offset){
  # browser()
  optim_res = optim(rep(0,ncol(X)), fn = logpost_offset, gr = logpost_gr_offset, y = y, X = X, prior_sd = 10, offset = offset,
                    method = "BFGS" )
  beta = optim_res$par
  
  eta = c(X%*%beta) + offset
  p = 1/(1+exp(-eta))
  
  S = t(X)%*%diag(p*(1-p))%*%X + 1/prior_sd^2*diag(length(beta))
  -0.5*determinant(S, logarithm = T)$modulus + ncol(S)/2*log(2*pi) - optim_res$value
}

marg_lik_laplace_prob_offset = function(y, X, prior_sd, offset){
  # browser()
  optim_res = optim(rep(0,ncol(X)), fn = logpost_offset, gr = logpost_gr_offset, y = y, X = X, 
                    prior_sd = 10, offset = offset, method = "BFGS" )
  beta = optim_res$par
  
  eta = c(X%*%beta) + offset
  p = 1/(1+exp(-eta))
  
  Sigma = solve(t(X)%*%diag(p*(1-p))%*%X + 1/prior_sd^2*diag(length(beta)))
  
  # get index
  all_vars = colnames(X)
  idx = c(which(all_vars == "casectrl"), 
          which(all_vars == "kps90_100:casectrl"), 
          which(all_vars == "mgmt:casectrl"), 
          which(all_vars == "kps90_100:mgmt:casectrl"))
  
  beta_use = beta[idx]
  Sigma_use = Sigma[idx,idx]
  contrast_mat = model.matrix( ~ x2*x3-1, expand.grid(x2 = factor(0:1), x3 = factor(0:1)))
  contrast_mat[,1] = 1
  beta_new = contrast_mat%*%beta_use
  Sigma_new = contrast_mat%*%Sigma_use%*%t(contrast_mat)
  
  1 - mvtnorm::pmvnorm(lower = rep(-Inf, 4), upper = rep(0, 4), mean = c(beta_new), sigma = Sigma_new)
}

permute_driver_offset = function(data_in, prior_sd, beta_e, self = F, prob = F, n_perm = 1e3){
  # browser()
  # eor and kps90_100:mgmt are the two terms need refitting
  # beta_e is for intercept, Age, Sex, eor, kps90_100, mgmt, kps90_100:mgmt
  X_raw = model.matrix( outcome ~ Age + Sex + eor + kps90_100*mgmt*casectrl, data = data_in)
  offset = c(X_raw[,c("(Intercept)", "Age", "Sex", "eorSTR", "kps90_100", "mgmt", "kps90_100:mgmt")]%*%beta_e)
  
  if(self & prob == FALSE){
    X = X_raw[,c("casectrl", "kps90_100:casectrl", "mgmt:casectrl", "kps90_100:mgmt:casectrl")]
  }else if(self & prob){
    X = X_raw[,c("(Intercept)", "casectrl", "kps90_100:casectrl", "mgmt:casectrl", "kps90_100:mgmt:casectrl")]
  }else{
    X = X_raw[,c("(Intercept)","eorSTR", "kps90_100:mgmt", "casectrl", "kps90_100:casectrl", 
                 "mgmt:casectrl", "kps90_100:mgmt:casectrl")]
  }
  y = data_in$outcome
  
  if(prob){
    marg_lik_obs = marg_lik_laplace_prob_offset(y = y, X = X, prior_sd = prior_sd, offset = offset)
  }else{
    marg_lik_obs = marg_lik_laplace_map_offset(y = y, X = X, prior_sd = prior_sd, offset = offset)
  }
  
  in_casectrl = data_in$casectrl
  marg_lik_perm = sapply(1:n_perm, function(rep){
    data_in$casectrl = gtools::permute(in_casectrl)
    if(self & prob == FALSE){
      X = model.matrix( outcome ~ casectrl + 
                          kps90_100:casectrl + mgmt:casectrl +
                          kps90_100:mgmt:casectrl, data = data_in)[,-1]
      colnames(X) = c("casectrl", "kps90_100:casectrl", "mgmt:casectrl",
                          "kps90_100:mgmt:casectrl")
    }else if(self & prob){
      X = model.matrix( outcome ~ casectrl + 
                          kps90_100:casectrl + mgmt:casectrl +
                          kps90_100:mgmt:casectrl, data = data_in)
      colnames(X)[-1] = c("casectrl", "kps90_100:casectrl", "mgmt:casectrl",
                      "kps90_100:mgmt:casectrl")
    }
    else{
      X = model.matrix( outcome ~ eor + kps90_100:mgmt + casectrl + 
                          kps90_100:casectrl + mgmt:casectrl +
                          kps90_100:mgmt:casectrl, data = data_in )#[,-1]
    }
    if(prob){
      marg_lik_obs = marg_lik_laplace_prob_offset(y = y, X = X, prior_sd = prior_sd, offset = offset)
    }else{
      marg_lik_obs = marg_lik_laplace_map_offset(y = y, X = X, prior_sd = prior_sd, offset = offset)
    }
  })
  mean(marg_lik_perm > marg_lik_obs)
}
