library(dplyr)
library(survival)
library(gtools)
set.seed(1)
source("permutation_src.R")

# we don't do fisher's p-val combination for now
data_raw = read.csv("GBM-Data-Tools/Full Excel Export of the GBM ECA v1 REDCAP 10-20-22.csv")
data_raw = data_raw %>% filter( !is.na(Sex) & !is.na(Age) & !is.na(kps90_100) & !is.na(eor) & !is.na(mgmt))
data_raw_lst = split(data_raw, data_raw$database)

data_in = data_raw_lst$AVAGLIO
data_out = data_raw_lst$DFCI
data_out = data_out %>% filter(StandardRT_TMZ == 1)

data_in_use = impute_data_set(data_in) %>% filter(eor != "BX")
data_out_use = impute_data_set(data_out) %>% filter(eor != "BX")

data_in_use$subtype = as.factor(4 - data_in_use$kps90_100 - 2*data_in_use$mgmt)
data_out_use$subtype = as.factor(4 - data_out_use$kps90_100 - 2*data_out_use$mgmt)

match_custom = function(data_in, data_out){
  # match on subtype, sex, eor and nearest neighbor on age
  data_1 = data_in[data_in$casectrl==1,]
  data_ctrl = rbind(data_in[data_in$casectrl==0,], data_out)
    
  Y_0 = sapply(1:nrow(data_1), function(idx){
    x = data_1[idx,]
    data_2 = data_ctrl[data_ctrl$subtype == x[[8]] & data_ctrl$Sex == x[[3]] & data_ctrl$eor == x[[5]],]
    if(nrow(data_2)==0){
      data_2 = data_ctrl[data_ctrl$subtype == x[[8]] & data_ctrl$Sex == x[[3]],]
    }
    idx = order(abs(data_2$Age - x[[2]]))[1]
    data_2$outcome[idx]
  })

  test_res <- tryCatch( mcnemar.test(as.factor(data_1$outcome), as.factor(Y_0))$p.value,
                        error = function(e) prop.test(x = c(sum(data_1$outcome), sum(Y_0)), 
                                                      n = c(nrow(data_1), length(Y_0)))$p.value )
  test_res
}

ipw_fun = function(data_in, data_out, type = c("weight", "regression")){
  # browser()
  data_merge = rbind(data_in, data_out)
  # ipw step
  weight_model = ipw::ipwpoint(exposure = source, family = "binomial", link = "logit",
           denominator =  ~ eor + Age + Sex + mgmt + kps90_100, data = data_merge)
  w_all = weight_model$ipw.weights
  w_all[data_merge$source==0] = (1-1/w_all[data_merge$source==0])*w_all[data_merge$source==0]*
    sum(1 - data_merge$source)/sum(data_merge$source)
  w_all[data_merge$source==0] = 1
  
  # data_merge$weights = w_all
  # w_all[data_merge$casectrl==0] = w_all[data_merge$casectrl==0]/2
  
  if(type == "weight"){
    # model_res = survey::svyglm(outcome ~ casectrl, 
    #                    design = svydesign( ~ 1, weights = ~ weights, data = data_merge))
    # summary(model_res)$coefficients[2,4]
    model_res = lm(outcome ~ casectrl, data = data_merge, weights = w_all)
    sd_res = sqrt(sandwich::vcovHC(model_res)[2,2])
  }else{
    # model_res = survey::svyglm(outcome ~ casectrl + eor + Age + Sex, 
    #                    design = svydesign( ~ 1, weights = ~ weights, data = data_merge))
    # summary(model_res)$coefficients[2,4]
    model_res = glm(outcome ~ casectrl + eor + Age + Sex, data = data_merge, weights = w_all)
    sd_res = sqrt(sandwich::vcovHC(model_res)[2,2])
  }
  t = abs(coef(model_res)[2]/sd_res)
  # c(coef(model_res)[2], 2*(1-pnorm(t)))
  2*(1-pnorm(t))
}

fisher_combine = function(p_all){
  res = sum(-2*log(p_all))
  1 - pchisq(res, df = 2*length(p_all))
}

library(snowfall)
sfInit(parallel = T, cpus = 8)
sfLibrary(dplyr)
sfLibrary(survey)
sfLibrary(MatchIt)
sfLibrary(marginaleffects)
sfSource("permutation_src.R")
sfExport("ipw_fun")
p_all_dfci = lapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "ext")
    data_in_lst = split(data_all$data_in, data_all$data_in$subtype)
    data_out_lst = split(data_all$data_out, data_all$data_out$subtype)
    
    data_merge = rbind(data_all$data_in, data_all$data_out)
    t1 = MatchIt::matchit(casectrl ~ Age + Sex + kps90_100 + eor + mgmt, data = data_merge,
                     method = "nearest", distance = "glm")
    t1_d = match.data(t1)
    fit = lm(outcome ~ casectrl*(Age + Sex + kps90_100 + eor + mgmt), data = t1_d, weights = weights)
    res = avg_comparisons(fit, variables = "casectrl", vcov = ~subclass, newdata = subset(t1_d, casectrl == 1),
                    wts = "weights")
    res$p.value
  })
})

p_all_dfci_neg = lapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "ext",
                                    lor = c(-2,0,0,0), alt = T)
    data_in_lst = split(data_all$data_in, data_all$data_in$subtype)
    data_out_lst = split(data_all$data_out, data_all$data_out$subtype)
    
    data_merge = rbind(data_all$data_in, data_all$data_out)
    t1 = MatchIt::matchit(casectrl ~ Age + Sex + kps90_100 + eor + mgmt, data = data_merge,
                          method = "nearest", distance = "glm")
    t1_d = match.data(t1)
    fit = lm(outcome ~ casectrl*(Age + Sex + kps90_100 + eor + mgmt), data = t1_d, weights = weights)
    res = avg_comparisons(fit, variables = "casectrl", vcov = ~subclass, newdata = subset(t1_d, casectrl == 1),
                          wts = "weights")
    res$p.value
  })
})

p_all_ava_neg = lapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "self",
                                    lor = c(-2,0,0,0), alt = T)
    data_in_lst = split(data_all$data_in, data_all$data_in$subtype)
    data_out_lst = split(data_all$data_out, data_all$data_out$subtype)
    
    data_merge = rbind(data_all$data_in, data_all$data_out)
    t1 = MatchIt::matchit(casectrl ~ Age + Sex + kps90_100 + eor + mgmt, data = data_merge,
                          method = "nearest", distance = "glm")
    t1_d = match.data(t1)
    fit = lm(outcome ~ casectrl*(Age + Sex + kps90_100 + eor + mgmt), data = t1_d, weights = weights)
    res = avg_comparisons(fit, variables = "casectrl", vcov = ~subclass, newdata = subset(t1_d, casectrl == 1),
                          wts = "weights")
    res$p.value
  })
})

sapply(p_all_dfci_neg, function(x) mean(x<0.05))
sapply(p_all_ava_neg, function(x) mean(x<0.05))

p_all_ava = lapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "self")
    data_in_lst = split(data_all$data_in, data_all$data_in$subtype)
    data_out_lst = split(data_all$data_out, data_all$data_out$subtype)
    
    data_merge = rbind(data_all$data_in, data_all$data_out)
    t1 = MatchIt::matchit(casectrl ~ Age + Sex + kps90_100 + eor + mgmt, data = data_merge,
                          method = "nearest", distance = "glm")
    t1_d = match.data(t1)
    fit = lm(outcome ~ casectrl*(Age + Sex + kps90_100 + eor + mgmt), data = t1_d, weights = weights)
    res = avg_comparisons(fit, variables = "casectrl", vcov = ~subclass, newdata = subset(t1_d, casectrl == 1),
                          wts = "weights")
    res$p.value
  })
})

p_all_dfci_ipw = sapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "ext")
    ipw_fun(data_all$data_in, data_all$data_out, type = "weight")
  })
})
colMeans(p_all_dfci_ipw<0.05)

p_all_dfci_ipw_neg = sapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "ext",
                                    lor = c(-2,0,0,0), alt = T)
    ipw_fun(data_all$data_in, data_all$data_out, type = "weight")
  })
})
colMeans(p_all_dfci_ipw_neg<0.05)

p_all_ava_ipw = sapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "self")
    ipw_fun(data_all$data_in, data_all$data_out, type = "weight")
  })
})
colMeans(p_all_ava_ipw<0.05)

p_all_ava_ipw_neg = sapply(seq(50,250,50), function(n){
  print(n)
  sfSapply(1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n, ratio = 2/3, out = "self",
                                    lor = c(-2,0,0,0), alt = T)
    ipw_fun(data_all$data_in, data_all$data_out, type = "weight")
  })
})
colMeans(p_all_ava_ipw_neg<0.05)

pwr_neg_causal = rbind(sapply(p_all_dfci_neg, function(x) mean(x<0.05)),
                       sapply(p_all_ava_neg, function(x) mean(x<0.05)),
                       colMeans(p_all_dfci_ipw_neg<0.05),
                       colMeans(p_all_ava_ipw_neg<0.05))
pwr_neg_causal_df = data.frame(power = c(pwr_neg_causal),
                               method = c("Test-Matching", "Test-Matching", "Test-IPW", "Test-IPW"),
                               source = c("ID + ED"),
                               ne = rep(seq(50,250,50), each = 4),
                               scenario = "6. Negative (-2,0,0,0)",
                               ED = rep(c("DFCI", "AVAGLIO")))
