or_calc = function(p1, p2){
  log(p1/(1-p1)) - log((p2/(1-p2)))
}

# imputation, AVAGLIO, other
impute_data_set = function(dataset){
  fit = coxph(Surv(dataset$os, dataset$os_status) ~ Sex + Age + kps90_100 + eor + mgmt, data = dataset)
  data_imp = dataset %>% filter( os_status == 0 & os < 365 )
  
  surv_pred = survfit(fit, newdata = data_imp)
  idx_time = which(surv_pred$time == 365)
  
  outcome_imp = 1 - sapply(1:nrow(data_imp), function(idx){
    obs_time = data_imp$os[idx]
    p = surv_pred$surv[,idx][idx_time]/surv_pred$surv[,idx][which(surv_pred$time==obs_time)]
    sample(0:1, 1, prob = c(p, 1-p))
  })
  
  out1 = dataset %>% filter(!(os_status == 0 & os < 365))
  out1$outcome = 1*(out1$os>=365)
  out1 = out1 %>% select(record, Age, Sex, kps90_100, eor, mgmt, outcome)
  out2 = data_imp %>% select(record, Age, Sex, kps90_100, eor, mgmt)
  out2$outcome = outcome_imp
  
  rbind(out1, out2)
}

# permutation module
compute_marg_lik = function(y, eta_null, X_perm, beta_perm){
  # calculate marginal likelihood
  # browser()
  eta = X_perm%*%t(beta_perm) + eta_null
  logp = -log(1+exp(-eta))
  
  mean(colSums(logp - eta*(1-y)))
}

marg_lik_laplace = function(glm_res, prior_sd){
  # browser()
  Sigma = vcov(glm_res)
  0.5*determinant(Sigma, logarithm = T)$modulus + ncol(Sigma)/2*log(2*pi) + 
    logLik(glm_res) + sum(dnorm(x = coef(glm_res), mean = 0, sd = prior_sd, log = T))
}

logpost = function(beta, y, X, prior_sd){
  eta = X%*%beta
  -(sum(y*eta) - sum(log(1+exp(eta))) - sum(beta^2)/2/prior_sd^2 - (log(prior_sd) + log(2*pi)/2)*length(beta))
}
logpost_gr = function(beta, y, X, prior_sd){
  eta = c(X%*%beta)
  -(t(X)%*%y - colSums(1/(1+exp(-eta))*X)) + beta/prior_sd^2
}

# ED-PT m
marg_lik_laplace_map = function(y, X, prior_sd){
  # browser()
  optim_res = optim(rep(0,ncol(X)), fn = logpost, gr = logpost_gr, y = y, X = X, prior_sd = 10, method = "BFGS" )
  beta = optim_res$par
  
  eta = c(X%*%beta)
  p = 1/(1+exp(-eta))
  
  S = t(X)%*%(X*p*(1-p)) + 1/prior_sd^2*diag(length(beta))
  -0.5*determinant(S, logarithm = T)$modulus + ncol(S)/2*log(2*pi) - optim_res$value
}

# ED-PT mtilde1
marg_lik_laplace_prob = function(y, X, prior_sd, threshold = 0){
  # browser()
  optim_res = optim(rep(0,ncol(X)), fn = logpost, gr = logpost_gr, y = y, X = X, prior_sd = 10, method = "BFGS" )
  beta = optim_res$par
  
  eta = c(X%*%beta)
  p = 1/(1+exp(-eta))
  
  Sigma = chol2inv(chol(t(X)%*%(X*p*(1-p)) + 1/prior_sd^2*diag(length(beta))))
  
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
  
  1 - mvtnorm::pmvnorm(lower = rep(-Inf, 4), upper = rep(threshold, 4), mean = c(beta_new), sigma = Sigma_new)
}

# ED-PT mtilde2
marg_lik_laplace_regret = function(y, X, p_all, prior_sd, n_idx, n_mc = 1e3){
  # browser()
  optim_res = optim(rep(0,ncol(X)), fn = logpost, gr = logpost_gr, y = y, X = X, prior_sd = 10, method = "BFGS" )
  beta = optim_res$par
  p = c(1/(1+exp(-X%*%beta)))
  Sigma = chol2inv(chol(t(X)%*%(X*p*(1-p)) + 1/prior_sd^2*diag(length(beta))))
  beta_sample = MASS::mvrnorm(n = n_mc, mu = beta, Sigma = Sigma)
  
  # get index
  all_vars = colnames(X)
  idx = c(which(all_vars == "casectrl"), 
          which(all_vars == "kps90_100:casectrl"), 
          which(all_vars == "mgmt:casectrl"), 
          which(all_vars == "kps90_100:mgmt:casectrl"))
  
  # compute p for treated and control
  X_use = X[1:n_idx,]
  eta1 = X_use%*%t(beta_sample)
  X_use[,idx] = 0
  eta2 = X_use%*%t(beta_sample)
  
  p1 = 1/(1+exp(-eta1))
  p2 = 1/(1+exp(-eta2))
  
  pdiff = p1 - p2
  pdiff[pdiff<0] = 0
  
  mean(pdiff)
}

permute_driver = function(data_merge, prior_sd, self = F, prob = c("marginal", "prob", "regret"), 
                          threshold = 0, p_all = NA, n_perm = 1e3){
  # browser()
  if(length(unique(data_merge$source)) == 1 | self & prob != "prob" ){
    X = model.matrix( outcome ~ (Age + Sex + eor) + kps90_100*mgmt*casectrl, data = data_merge)
  }else if(self & prob == "prob"){
    X = model.matrix( outcome ~ source + (Age + Sex + eor) + kps90_100*mgmt*casectrl, data = data_merge)
  }else{
    X = model.matrix( outcome ~ Age + Sex + eor + eor:source + kps90_100 + mgmt + kps90_100:mgmt:source + 
                        kps90_100*mgmt*casectrl, data = data_merge)
  }
  y = data_merge$outcome
  if(prob == "prob"){
    marg_lik_obs = marg_lik_laplace_prob(y = y, X = X, prior_sd = prior_sd, threshold = threshold)
  }else if(prob == "marginal"){
    marg_lik_obs = marg_lik_laplace_map(y = y, X = X, prior_sd = prior_sd)
  }else{
    marg_lik_obs = marg_lik_laplace_regret(y = y, X = X, prior_sd = prior_sd, p_all = p_all, n_idx = sum(data_merge$source==1) )
  }
  
  in_idx = which(data_merge$source == 1)
  in_casectrl = data_merge$casectrl[in_idx]
  marg_lik_perm = sapply(1:n_perm, function(rep){
    data_merge$casectrl[in_idx] = gtools::permute(in_casectrl)
    if(length(unique(data_merge$source)) == 1 | self == T & prob != "prob"){
      X = model.matrix( outcome ~ (Age + Sex + eor) + kps90_100*mgmt*casectrl, data = data_merge)
    }else if( self & prob == "prob"){
      X = model.matrix( outcome ~ source + (Age + Sex + eor) + kps90_100*mgmt*casectrl, data = data_merge)
    }else{
      X = model.matrix( outcome ~ Age + Sex + eor + eor:source + kps90_100 + mgmt + kps90_100:mgmt:source + 
                          kps90_100*mgmt*casectrl, data = data_merge)
    }
    if(prob == "prob"){
      marg_lik_obs = marg_lik_laplace_prob(y = y, X = X, prior_sd = prior_sd, threshold = threshold)
    }else if(prob == "marginal"){
      marg_lik_obs = marg_lik_laplace_map(y = y, X = X, prior_sd = prior_sd)
    }else{
      marg_lik_obs = marg_lik_laplace_regret(y = y, X = X, prior_sd = prior_sd, p_all = p_all,
                                             n_idx = sum(data_merge$source==1))
    }
  })
  mean(marg_lik_perm > marg_lik_obs)
}

# data generation
data_generate_module = function(data_in_use, data_out_use, n_in, n_out, 
                                ratio, alt = F, lor = rep(0,4), out = c('self', "ext")){
  n_case = round( n_in * ratio )
  n_ctrl = round( n_in * (1-ratio))
  
  data_in = rbind( data_in_use[sample(1:nrow(data_in_use), n_case, replace = T),],
                   data_in_use[sample(1:nrow(data_in_use), n_ctrl, replace = T),] )
  data_in$casectrl = rep(c(1,0), c(n_case, n_ctrl))
  data_in$source = 1
  
  p_orig = (data_in %>% group_by(subtype) %>% summarise(p = mean(outcome)))$p
  odds_orig = p_orig/(1-p_orig)
  p_new = sapply(seq_along(lor), function(i){
    if(lor[i]>0){
      (exp(lor[i])-1)*odds_orig[i]/(exp(lor[i])*odds_orig[i] + 1)
    }else{
      (1-exp(lor[i]))/(1 + odds_orig[i]*exp(lor[i]))
    }
  })
  
  # tweaking subgroup-specific treatment effect
  # group-specific
  if(alt){
    for(group in 1:4){
      if(lor[group] > 0){
        idx =  which(data_in$outcome == 0 & data_in$casectrl == 1 & data_in$subtype == group)
        new_Y = rbinom(length(idx), size = 1, p = p_new[group])
        data_in$outcome[idx] = new_Y
      }else{
        idx = which(data_in$outcome == 1 & data_in$casectrl == 1 & data_in$subtype == group)
        new_Y = rbinom(length(idx), size = 1, p = p_new[group])
        data_in$outcome[idx] = 1 - new_Y
      }
    }
  }
  
  if(out != "self"){
    data_out = data_out_use[sample(1:nrow(data_out_use), n_out, replace = T),]
  }else{
    data_out = data_in_use[sample(1:nrow(data_in_use), n_out, replace = T),]
  }
  data_out$casectrl = 0
  data_out$source = 0
  
  
  list(data_in = data_in, data_out = data_out)
}

permute_test = function(data_all, prior_sd, ID = F, self = F, prob = c("marginal", "prob", "regret"), 
                        threshold = 0, n_perm = 100){
  # browser()
  data_in = data_all$data_in
  data_out = data_all$data_out
  data_merge = rbind(data_in, data_out)
  
  coding = data_in$kps90_100 + 2*data_in$mgmt
  p_all = table(coding)/length(coding)
  
  if(ID){
    p_val = permute_driver(data_in, prior_sd, self = self, n_perm = n_perm, prob = prob, threshold = threshold, 
                           p_all = p_all)
  }else{
    p_val = permute_driver(data_merge, prior_sd, self = self, n_perm = n_perm, prob = prob, threshold = threshold,
                           p_all = p_all)
  }
  # p_val = permute_driver(data_in, prior_sd, n_perm = n_perm)
  
  p_val
}

other_test = function(data_all){
  data_in = data_all$data_in
  data_out = data_all$data_out
  data_merge = rbind(data_in, data_out)
  
  # # no ED LR test
  model_null = glm( outcome ~ Age + Sex + eor + kps90_100*mgmt, family = "binomial", data = data_in)
  model_full = glm( outcome ~ Age + Sex + eor + kps90_100*mgmt*casectrl, family = "binomial", data = data_in)
  p_in = anova(model_null, model_full, test = "Chisq")$`Pr(>Chi)`[2]
  
  # with ED LR test
  model_null_m = glm( outcome ~ Age + Sex + eor + kps90_100*mgmt, family = "binomial", data = data_merge)
  model_full_m = glm( outcome ~ Age + Sex + eor + kps90_100*mgmt*casectrl, family = "binomial", data = data_merge)
  p_m = anova(model_null_m, model_full_m, test = "Chisq")$`Pr(>Chi)`[2]
  
  # ID Z test
  Z_ID = prop.test(x = c(sum(data_in$outcome*data_in$casectrl), sum(data_in$outcome*(1-data_in$casectrl))), 
            n = c(sum(data_in$casectrl), sum(1-data_in$casectrl)))
  p_in_z = Z_ID$p.value
  
  # ID + ED Z test
  Z_m = prop.test(x = c(sum(data_merge$outcome*data_merge$casectrl), sum(data_merge$outcome*(1-data_merge$casectrl))), 
                   n = c(sum(data_merge$casectrl), sum(1-data_merge$casectrl)))
  p_m_z = Z_m$p.value
  
  c(p_in, p_m, p_in_z, p_m_z)
}

matching_test = function(data_all){
  data_merge = rbind(data_all$data_in, data_all$data_out)
  t1 = MatchIt::matchit(casectrl ~ Age + Sex + kps90_100 + eor + mgmt, data = data_merge,
                        method = "nearest", distance = "glm")
  t1_d = MatchIt::match.data(t1)
  fit = glm(outcome ~ casectrl*(Age + Sex + kps90_100 + eor + mgmt), family = "binomial", 
            data = t1_d, weights = weights)
  res = marginaleffects::avg_comparisons(fit, variables = "casectrl", vcov = ~subclass, 
                                         newdata = subset(t1_d, casectrl == 1),
                                         wts = "weights")
  res$p.value
}

ipw_test = function(data_all, type = c("weight", "regression")){
  data_merge = rbind(data_all$data_in, data_all$data_out)
  # ipw step
  weight_model = ipw::ipwpoint(exposure = source, family = "binomial", link = "logit",
                               denominator =  ~ eor + Age + Sex + mgmt + kps90_100, data = data_merge)
  w_all = weight_model$ipw.weights
  w_all[data_merge$source==0] = (1-1/w_all[data_merge$source==0])*w_all[data_merge$source==0]*
    sum(1 - data_merge$source)/sum(data_merge$source)
  w_all[data_merge$source==0] = 1
  
  if(type == "weight"){
    model_res = lm(outcome ~ casectrl, data = data_merge, weights = w_all)
    sd_res = sqrt(sandwich::vcovHC(model_res)[2,2])
  }else{
    model_res = lm(outcome ~ casectrl + eor + Age + Sex, data = data_merge, weights = w_all)
    sd_res = sqrt(sandwich::vcovHC(model_res)[2,2])
  }
  t = abs(coef(model_res)[2]/sd_res)
  2*(1-pnorm(t))
}