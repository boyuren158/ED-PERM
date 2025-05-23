# functions for S1
data_sim = function( n_all, r_all, re_all, mu0, gamma, beta_1, gamma_1, mu0e ){
  # browser()
  K_all = seq_along(n_all)
  group_all = rep(K_all, n_all*(1 + r_all))
  a_all = rep(c(1,0), c(sum(n_all), sum(n_all*r_all)))
  a_all = dqrng::dqsample(a_all)
  D_in = data.frame(a = a_all, group = as.factor(group_all))
  D_mt = model.matrix( ~ a*group, data = D_in )
  coef_all = c(mu0, gamma, beta_1, gamma_1)
  y = rnorm(n = nrow(D_in), mean = D_mt%*%coef_all, sd = 1)
  y_nocov = y - model.matrix( ~ group, data = D_in)[,-1,drop = F]%*%beta_1 - mu0
  D_in$y = y
  
  beta1_all = c(0, beta_1)
  ye_all = lapply( K_all, function(k){
    rnorm(n = n_all[k]*re_all[k], mean = mu0 + mu0e + beta1_all[k], sd = 1)
  })
  list(D_in = D_in, ye_all = ye_all,
       y_nocov = y_nocov)
}

# generic m
m_calc = function( y_in, a_in, X_in, y_E, X_E, sigma2 ){
  Xt = cbind(1, X_in, a_in, X_in*a_in)
  XEt = cbind(1, X_E)
  Vinv = t(Xt)%*%Xt + diag(rep(1/sigma2, ncol(Xt)))
  Vinv[1:(1+ncol(X_in)), 1:(1+ncol(X_in))] = Vinv[1:(1+ncol(X_in)), 1:(1+ncol(X_in))] +
    t(XEt)%*%XEt
  mu = t(Xt)%*%y_in
  mu[1:(1+ncol(X_in))] = t(XEt)%*%y_E + mu[1:(1+ncol(X_in))]
  c(1/2*t(mu)%*%chol2inv(chol(Vinv))%*%mu - 1/2*determinant(Vinv, logarithm = T)$modulus)
}

# ID only statistic
m_calc_in = function(a, x, y, sigma2){
  X = cbind(1, a, x, a*x)
  Vinv = t(X)%*%X + 1/sigma2*diag(ncol(X))
  Z = t(X)%*%y
  c(1/2*t(Z)%*%solve(Vinv)%*%Z) - 1/2*determinant(Vinv, logarithm = T)$modulus
}

permute_test_sg = function(data, sigma2, n_perm){
  # browser()
  D_in = data$D_in
  X_in = model.matrix( ~ x, data = data.frame(x = as.factor(D_in$group)) )[,-1,drop = F]
  X_E = model.matrix( ~ x, data = data.frame(x = as.factor(rep(1:length(data$ye_all), 
                                                               sapply(data$ye_all, length)))) )[,-1,drop = F]
  m_obs = m_calc(y_in = D_in$y, a_in = D_in$a, X_in = X_in,
                 y_E = unlist(data$ye_all), X_E = X_E, sigma2 = sigma2)
  m_obs_ID = m_calc_in(D_in$a, X_in, D_in$y, sigma2)
  
  m_perm = sapply(1:n_perm, function(perm){
    a_new = dqrng::dqsample(D_in$a, size = length(D_in$a))
    m1 = m_calc( y_in = D_in$y, a_in = a_new, X_in = X_in,
                 y_E = unlist(data$ye_all), X_E = X_E, sigma2 = sigma2 )
    m2 = m_calc_in(a_new, X_in, D_in$y, sigma2)
    c(m1, m2)
  })
  
  p_perm = mean(m_perm[1,] > m_obs)
  p_perm_id = mean(m_perm[2,] > m_obs_ID)
  
  #alternative test
  #ID only Wald-test
  m1_id = lm(y ~ a*as.factor(group), data = data$D_in)
  S_id = solve(summary(m1_id)$cov.unscaled[c(2,4), c(2,4)])
  chi2_id = coef(m1_id)[c(2,4)]%*%(S_id%*%coef(m1_id)[c(2,4)])
  p_f_id = 1 - pchisq(chi2_id, df = 2)
  
  #ID + ED Wald-test
  data_merge = rbind( data$D_in, do.call(rbind, lapply(seq_along(data$ye_all), function(x){
    data.frame(a = 0, group = x, y = data$ye_all[[x]])
  })))
  data_merge$group = as.factor(data_merge$group)
  
  m1_ied = lm(y ~ a*group, data = data_merge)
  S_ied = solve(summary(m1_ied)$cov.unscaled[c(2,4), c(2,4)])
  chi2_ied = coef(m1_ied)[c(2,4)]%*%(S_ied%*%coef(m1_ied)[c(2,4)])
  p_f_ied = 1 - pchisq(chi2_ied, df = 2)
  
  #orcale
  in_tmp = data$D_in
  in_tmp$y = data$y_nocov
  data_o = in_tmp[in_tmp$a == 1,]
  data_o$group = as.factor(data_o$group)
  m1_o = lm(y ~ group, data = data_o)
  M_o = model.matrix(m1_o)
  S_o = t(M_o)%*%M_o
  chi2_o = coef(m1_o)%*%(S_o%*%coef(m1_o))
  p_f_o = 1 - pchisq(chi2_o, df = 2)
  
  c(p_perm, p_perm_id, p_f_id, p_f_ied, p_f_o)
}

# functions for S2
m_calc_ex2 = function( y_in, a_in, X_in, y_E, X_E, sigma2 ){
  
  Xt = cbind(1, X_in, a_in, X_in[,1]*a_in)
  XEt = cbind(1, X_E)
  Vinv = t(Xt)%*%Xt + diag(rep(1/sigma2, ncol(Xt)))
  Vinv[1:(1+ncol(X_in)), 1:(1+ncol(X_in))] = Vinv[1:(1+ncol(X_in)), 1:(1+ncol(X_in))] +
    t(XEt)%*%XEt
  mu = t(Xt)%*%y_in
  mu[1:(1+ncol(X_in))] = t(XEt)%*%y_E + mu[1:(1+ncol(X_in))]
  c(1/2*t(mu)%*%chol2inv(chol(Vinv))%*%mu - 1/2*determinant(Vinv, logarithm = T)$modulus)
}

m_calc_in_ex2 = function(a, x, y, sigma2){
  X = cbind(1, a, a*x[,1], x)
  Vinv = t(X)%*%X + 1/sigma2*diag(ncol(X))
  Z = t(X)%*%y
  c(1/2*t(Z)%*%solve(Vinv)%*%Z) - 1/2*determinant(Vinv, logarithm = T)$modulus
}

data_sim_continuous_subg = function( n, r, re, rho, p, mu0, gamma, beta_1, gamma_1, mu0e ){
  # beta_1 is for the X
  # gamma_1 is for X*a
  # gamma is for a
  X = cbind( rep(0:1, c(n*(1+r)*rho[1], n*(1+r)*rho[2])), matrix(rnorm(p*n*(1+r)), ncol = p))
  X_E = cbind( rep(0:1, c(n*re*rho[1], n*re*rho[2])), matrix(rnorm(p*n*re), ncol = p))
  
  a = c(rep(0:1, c(n*rho[1], n*r*rho[1])),
        rep(0:1, c(n*rho[2], n*r*rho[2])))
  a_E = rep(0, n*re)
  
  y = cbind(1, X, a, X*a)%*%c(mu0, beta_1, gamma, gamma_1) + rnorm(n = n*(1+r))
  y_E = cbind(1, X_E)%*%c(mu0, beta_1) + mu0e + rnorm(n = n*re)
  
  return(list(D_in = list(a = a, y = y, X = X),
              D_E = list(y_E = y_E, X_E = X_E),
              y_o = y - X%*%beta_1 - mu0))
}

permute_test_continuous_subg = function(data, sigma2, n_perm){
  # browser()
  D_in = data$D_in
  D_E = data$D_E
  
  m_obs = m_calc_ex2(y_in = D_in$y, a_in = D_in$a, X_in = D_in$X,
                     y_E = D_E$y_E, X_E = D_E$X_E, sigma2 = sigma2)
  m_obs_ID = m_calc_in_ex2(D_in$a, D_in$X, D_in$y, sigma2)
  
  m_perm = sapply(1:n_perm, function(perm){
    a_new = dqrng::dqsample(D_in$a, size = length(D_in$a))
    m1 = m_calc_ex2(y_in = D_in$y, a_in = a_new, X_in = D_in$X,
                    y_E = D_E$y_E, X_E = D_E$X_E, sigma2 = sigma2)
    m2 = m_calc_in_ex2(a_new, D_in$X, D_in$y, sigma2)
    c(m1, m2)
  })
  
  p_perm = mean(m_perm[1,] > m_obs)
  p_perm_id = mean(m_perm[2,] > m_obs_ID)
  
  #alternative test
  #ID only Wald test
  D_in_df = data.frame(X = D_in$X, y = D_in$y, a = D_in$a)
  m1_id = lm(y ~  (.-a) + a + a:X.1, data = D_in_df)
  idx = c(length(coef(m1_id))-1, length(coef(m1_id)))
  S_id = solve(summary(m1_id)$cov.unscaled[idx, idx])
  chi2_id = coef(m1_id)[idx]%*%(S_id%*%coef(m1_id)[idx])
  p_f_id = 1 - pchisq(chi2_id, df = 2)
  
  #ID z-test no covariates
  m2_id = lm(y ~ a*X.1, data = D_in_df)
  m2_id_reduced = lm(y ~ X.1, data = D_in_df)
  p_f_id_z = anova(m2_id_reduced, m2_id)$`Pr(>F)`[2]
  
  
  c(p_perm, p_perm_id, p_f_id, p_f_id_z)
}

# functions for negative effect examples
sim_asym = function(n, r, re, p_all, a_all, b_all, beta){
  # n is the size of experimental group in total
  K = length(p_all)
  # first all controls
  mean_ctrl = rep(beta, n*p_all*r)
  Y_ctrl = rnorm(length(mean_ctrl), mean = mean_ctrl, sd = 1)
  
  # then all exp
  mean_exp = rep(beta + a_all, n*p_all)
  Y_exp = rnorm(length(mean_exp), mean = mean_exp, sd = 1)
  
  # then ED
  mean_ctrl_ED = rep(beta + b_all, n*p_all*re)
  Y_ctrl_ED = rnorm(length(mean_ctrl_ED), mean = mean_ctrl_ED, sd = 1)
  
  list(ID = data.frame(y = c(Y_ctrl, Y_exp),
                       a = rep(0:1, c(length(Y_ctrl), length(Y_exp))),
                       group = c(rep(1:K, n*p_all*r), rep(1:K, n*p_all))) %>% arrange(group, a),
       ED = data.frame(y = Y_ctrl_ED,
                       group = c(rep(1:K,n*p_all*re))))
}

m_statistic = function(y, a, ybar, yebar, r, re){
  (sum(y*a)/sum(a) - ((1+r)*ybar + re*yebar)/(1+r+re))^2
}

# mtilde 1
# probability in Theta-tilde
p_statistic = function(y, a, ybar, yebar, r, re, neg = F){
  n1 = sum(a)
  n = length(a)
  ne = n/(1+r)*re
  mean = ((n + ne)*sum(y*a)/n1 - (n*ybar + ne*yebar))/(n - n1 + ne)
  var = (n + ne)/n1/(n - n1 + ne)
  z = mean/sqrt(var)
  if(neg){
    pnorm(-z, log.p = T)
  }else{
    pnorm(z, log.p = T)
  }
}

perm_asym_sim = function(data, sigma2 = 10, n_perm = 1e3){
  # browser()
  X_in = model.matrix( ~ as.factor(group), data$ID )[,-1,drop = F]
  X_E = model.matrix( ~ as.factor(group), data$ED )[,-1,drop = F]
  m_obs = m_calc(y_in = data$ID$y, a_in = data$ID$a, X_in = X_in,
                   y_E = data$ED$y, X_E = X_E, sigma2 = sigma2)
  
  m_perm = sapply(1:n_perm, function(perm){
    a_new = dqrng::dqsample(data$ID$a, size = length(data$ID$a))
    m_calc(y_in = data$ID$y, a_in = a_new, X_in = X_in,
            y_E = data$ED$y, X_E = X_E, sigma2 = sigma2)
  })
  
  mean(m_perm > m_obs)
}

perm_asym_sim_pos = function(data, sigma2 = 10, n_perm = 1e3){
  all_group = unique(data$ID$group)
  r = 1/(sum(data$ID$a)/sum(1-data$ID$a))
  re = nrow(data$ED)/sum(data$ID$a)
  # each column is a group
  bar_all = sapply(all_group, function(i){
    c(mean(data$ID$y[data$ID$group==i]),
      mean(data$ED$y[data$ED$group==i]))
  })
  all_stat = sapply(all_group, function(i){
    n1 = sum(data$ID$a[data$ID$group == i])
    n = sum(data$ID$group == i)
    ne = sum(data$ED$group == i)
    m = m_statistic(data$ID$y[data$ID$group==i], data$ID$a[data$ID$group==i],
                    bar_all[1,i], bar_all[2,i], r, re) - log(n1) - log(n - n1 + ne)
    p = p_statistic(data$ID$y[data$ID$group==i], data$ID$a[data$ID$group==i],
                    bar_all[1,i], bar_all[2,i], r, re, neg = T)
    c(m, p)
  })
  
  m_obs = log1p(-expm1(sum(all_stat[2,])) - 1) # + sum(all_stat[1,])
  # sum(all_stat[1,]) + sum(all_stat[2,])
  
  m_perm = sapply(1:n_perm, function(iter){
    a_new = dqrng::dqsample(data$ID$a)
    all_stat = sapply(all_group, function(i){
      n = sum(data$ID$group == i)
      ne = sum(data$ED$group == i)
      n1 = sum(a_new[data$ID$group == i])
      m = m_statistic(data$ID$y[data$ID$group==i], a_new[data$ID$group==i],
                      bar_all[1,i], bar_all[2,i], r, re) - log(n1) - log(n - n1 + ne)
      p = p_statistic(data$ID$y[data$ID$group==i], a_new[data$ID$group==i],
                      bar_all[1,i], bar_all[2,i], r, re, neg = T)
      c(m, p)
    })
    
    m_obs =  log1p(-expm1(sum(all_stat[2,])) - 1) # + sum(all_stat[1,])
    # sum(all_stat[1,] + all_stat[2,])
  })
  
  mean(m_perm > m_obs)
}

# mtilde 2
dec_m = function(a, y, ne, yebar){
  n1 = sum(a)
  n = length(a)
  y1bar = sum(a*y)/n1
  pos_mean = ((n + ne)*y1bar - (sum(y) + ne*yebar))/(n - n1 + ne)
  var = (n + ne)/n1/(n - n1 + ne)
  
  truncnorm::etruncnorm(a = 0, mean = pos_mean, sd = sqrt(var))*(1 - pnorm(0, mean = pos_mean, sd = sqrt(var)))
}

dec_perm = function(data, n_perm){
  group_all = unique(data$ID$group)
  ne_all = table(data$ED$group)
  yebar_all = sapply(group_all, function(g) mean(data$ED$y[data$ED$group == g]))
  
  m_obs = sum(sapply(group_all, function(g){
    dec_m(data$ID$a[data$ID$group==g], data$ID$y[data$ID$group==g],
          ne_all[g], yebar_all[g])
  }))
  m_perm = sapply(1:n_perm, function(iter){
    a_new = dqrng::dqsample(data$ID$a)
    sum(sapply(group_all, function(g){
      dec_m(a_new[data$ID$group==g], data$ID$y[data$ID$group==g],
            ne_all[g], yebar_all[g])
    }))
  })
  mean(m_perm >= m_obs)
}