binary_marg_multi = function(n_pos, n_all){
  browser()
  if(is.null(nrow(n_pos))){
    n_pos = matrix(n_pos, nrow = 1)
    n_all = matrix(n_all, nrow = 1)
  }
  n_c = n_all[,2, drop = F] + n_all[,3, drop = F]
  s_c = n_pos[,2, drop = F] + n_pos[,3, drop = F]
  
  sum(-lchoose(n_all[,1], n_pos[,1]) - lchoose(n_c, s_c))
}

binary_sim_multi = function(n_tot, r_all, p_ctrl_all, p_delta_all, n_e, p_e){
  # n_tot: vector, # of subjects in each group
  # r_all: vector, randomization ratio in each group, control/case
  # p_ctrl_all: vector, ctrl response rate in each group
  # p_delta_all: vector, trt - ctrl response rate in each group
  # n_e: vector, # of ED in each group
  # p_e: vector, response rate in ED in each group
  # note: we only consider external controls
  
  df = data.frame(group = rep(seq_along(n_tot), n_tot),
                  trt = unlist(lapply(seq_along(r_all), function(idx){
                    r = r_all[idx]
                    p = 1/(r+1)
                    rbinom(n_tot[idx], 1, prob = p)
                  })))
  p_all = rep(p_ctrl_all, n_tot)
  p_delta = rep(p_delta_all, n_tot)
  p_all = p_all + p_delta*df$trt
  df$response = rbinom(sum(n_tot), size = 1, prob = p_all)
  
  list(df = df, e_all = n_e, e_pos = rbinom(n = length(p_e), size = n_e, prob = p_e))
}

permutation_binary_multi = function(df, e_pos, e_all, m = binary_marg_multi, n_perm = 1e3){
  # df is the individual level ID, colnames group, trt, response
  # e_pos, # of responders in ED, per group
  # e_all, # of subjects in ED, per group
  n_pos = cbind(table(df$group, factor(df$trt, levels = c(1,0)), df$response)[,,2], e_pos)
  n_all = cbind(table(df$group, factor(df$trt, levels = c(1,0))), e_all)
  m_obs = m(n_pos, n_all)
  m_perm = sapply(1:n_perm, function(rep){
    trt = factor(dqrng::dqsample(df$trt, size = nrow(df)), level = c(1,0))
    n_pos = cbind(table(df$group, trt, df$response)[,,2], e_pos)
    n_all = cbind(table(df$group, trt), e_all)
    m(n_pos, n_all)
  })
  mean(m_obs<=m_perm)
}

# type I error
library(snowfall)
sfInit(parallel = T, cpus = 8)
sfExport("binary_sim_multi", "permutation_binary_multi", "binary_marg_multi")
sfLibrary(dqrng)
system.time(type1e <- lapply(seq(-0.1,0.1,0.02), function(bias){
  print(bias)
  res = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0,0), n_e = c(375, 375), p_e = c(0.5+bias, 0.5+bias))
    p_perm = permutation_binary_multi(sim_res$df, sim_res$e_pos, sim_res$e_all, n_perm = 1e3)
    p_perm_ID = permutation_binary_multi(sim_res$df, rep(0, length(sim_res$e_all)), 
                                         rep(0, length(sim_res$e_pos)), n_perm = 1e3)
    m1 = glm( response ~ trt * group, family = "binomial", data = sim_res$df)
    m2 = glm( response ~ group, family = "binomial", data = sim_res$df)
    p_in_lrt = anova(m2, m1, test = "LRT")$`Pr(>Chi)`[2]

    library(dplyr)
    df_collapse = sim_res$df %>% group_by(trt, group) %>% 
      summarise(n_resp = sum(response), n = length(response))
    df_collapse$n_resp[df_collapse$trt==0] = df_collapse$n_resp[df_collapse$trt==0] + sim_res$e_pos
    df_collapse$n[df_collapse$trt==0] = df_collapse$n[df_collapse$trt==0] + sim_res$e_all
    df_collapse$y = df_collapse$n_resp/df_collapse$n
    
    m1_c = glm( y  ~ trt*group, family = "binomial", data = df_collapse, weights = n )
    m2_c = glm( y  ~ group, family = "binomial", data = df_collapse, weights = n )
    p_c_lrt = anova(m2_c, m1_c, test = "LRT")$`Pr(>Chi)`[2]
    c(p_perm, p_perm_ID, p_in_lrt, p_c_lrt)
  })
}))

oracle_multi = function(p_ctrl, data){
  # browser()
  dt = data.table::data.table(data)
  group_data = dt[,.(s = sum(response), n = length(response), p = mean(response)), by = group]
  test_res = prop.test(group_data$s, group_data$n, p = p_ctrl, correct = F)
  test_res$p.value
}

wald_multi = function(s_e, n_e, data){
  # browser()
  dt = data.table::data.table(data)
  group_data = dt[,.(s = sum(response), n = length(response), p = mean(response)), .(group, trt)]
  x_all = sapply(unique(group_data$group), function(i){
    tmp = group_data[group == i,]
    prop.test(x = c(tmp[trt == 1,]$s, tmp[trt == 0,]$s + s_e[i]),
              n = c(tmp[trt == 1,]$n, tmp[trt == 0,]$n + n_e[i]),
              correct = F)$statistic
  })
  pchisq(sum(x_all),df = length(s_e), lower.tail = F)
}

sfExport("oracle_multi")
sfExport("wald_multi")
system.time(type1e_oracle <- sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0,0), n_e = c(375, 375), p_e = c(0.5, 0.5))
    oracle_multi(c(0.5,0.5), sim_res$df[sim_res$df$trt == 1,])
  }))

mean(type1e_oracle<0.05)

system.time(type1e_wald <- lapply(seq(-0.1,0.1,0.02), function(bias){
  print(bias)
  res = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0,0), n_e = c(375, 375), p_e = c(0.5+bias, 0.5+bias))
    p_in_wald = wald_multi(c(0,0), c(0,0), sim_res$df)
    p_c_wald = wald_multi(sim_res$e_pos, sim_res$e_all, sim_res$df)
    c(p_in_wald, p_c_wald)
  })
}))

system.time(power <- lapply(seq(-0.1,0.1,0.02), function(bias){
  print(bias)
  res = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0.25,0.1), n_e = c(375, 375), p_e = c(0.5+bias, 0.5+bias))
    p_perm = permutation_binary_multi(sim_res$df, sim_res$e_pos, sim_res$e_all, n_perm = 1e3)
    p_perm_ID = permutation_binary_multi(sim_res$df, rep(0, length(sim_res$e_all)), 
                                         rep(0, length(sim_res$e_pos)), n_perm = 1e3)
    m1 = glm( response ~ trt * group, family = "binomial", data = sim_res$df)
    m2 = glm( response ~ group, family = "binomial", data = sim_res$df)
    p_in_lrt = anova(m2, m1, test = "LRT")$`Pr(>Chi)`[2]
    
    library(dplyr)
    df_collapse = sim_res$df %>% group_by(trt, group) %>% 
      summarise(n_resp = sum(response), n = length(response))
    df_collapse$n_resp[df_collapse$trt==0] = df_collapse$n_resp[df_collapse$trt==0] + sim_res$e_pos
    df_collapse$n[df_collapse$trt==0] = df_collapse$n[df_collapse$trt==0] + sim_res$e_all
    df_collapse$y = df_collapse$n_resp/df_collapse$n
    
    m1_c = glm( y  ~ trt*group, family = "binomial", data = df_collapse, weights = n )
    m2_c = glm( y  ~ group, family = "binomial", data = df_collapse, weights = n )
    p_c_lrt = anova(m2_c, m1_c, test = "LRT")$`Pr(>Chi)`[2]
    c(p_perm, p_perm_ID, p_in_lrt, p_c_lrt)
  })
}))

system.time(power <- lapply(seq(-0.1,0.1,0.02), function(bias){
  print(bias)
  res = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0.25,0.1), n_e = c(375, 375), p_e = c(0.5+bias, 0.5+bias))
    p_perm = permutation_binary_multi(sim_res$df, sim_res$e_pos, sim_res$e_all, n_perm = 1e3)
    p_perm_ID = permutation_binary_multi(sim_res$df, rep(0, length(sim_res$e_all)), 
                                         rep(0, length(sim_res$e_pos)), n_perm = 1e3)
    c(p_perm, p_perm_ID)
  })
}))

power_wald <- lapply(seq(-0.1,0.1,0.02), function(bias){
  print(bias)
  res = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0.25,0.1), n_e = c(375, 375), p_e = c(0.5+bias, 0.5+bias))
    p_in_wald = wald_multi(c(0,0), c(0,0), sim_res$df)
    p_c_wald = wald_multi(sim_res$e_pos, sim_res$e_all, sim_res$df)
    c(p_in_wald, p_c_wald)
  })
})

power_oracle = sfSapply(1:1e4, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0.25,0.1), n_e = c(375, 375), p_e = c(0.5, 0.5))
    oracle_multi(c(0.5,0.5), sim_res$df)
  })

res_05 = sapply(type1e, function(x) rowMeans(x<0.05))
res_05_wald = sapply(type1e_wald, function(x) rowMeans(x<0.05))
res_05[2,] = mean(res_05[2,])
res_05[3,] = mean(res_05_wald[1,])
res_05 = rbind(res_05,mean(type1e_oracle<0.05))

res_01 = sapply(type1e, function(x) rowMeans(x<0.01))
res_01_wald = sapply(type1e_wald, function(x) rowMeans(x<0.01))
res_01[2,] = mean(res_01[2,])
res_01[3,] = mean(res_01_wald[1,])
res_01 = rbind(res_01,mean(type1e_oracle<0.01))

power_05 = sapply(power, function(x) rowMeans(x<0.05))
power_05_wald = sapply(power_wald, function(x) rowMeans(x<0.05))
power_05[2,] = mean(power_05[2,])
power_05[3,] = mean(power_05_wald[1,])
power_05[4,] = power_05_wald[2,]
power_05 = rbind(power_05, mean(power_oracle < 0.05))

plt_type1e = data.frame(Method = c(rep(c("Permutation", "Wald"), each = 2), "Oracle"),
                        Data = factor(c(rep(c("ID", "ID + ED"),2), "ID"), levels = c("ID + ED", "ID")),
                        bias = rep(seq(-0.1,0.1,0.02), each = 5),
                        type1e = c(c(res_05), c(res_01)),
                        alpha = (rep(c(0.05, 0.01), each = length(res_05))))
plt_power = data.frame(Method = c(rep(c("Permutation", "Wald"), each = 2), "Oracle"),
                       Data = factor(c("ID + ED", "ID", "ID", "ID + ED", "ID"), levels = c("ID + ED", "ID")),
                       bias = rep(seq(-0.1,0.1,0.02), each = 5),
                       power = c(power_05))

library(ggplot2)
library(dplyr)
p1 = ggplot(plt_type1e %>% filter(alpha == 0.01), aes(x = bias, y = type1e, linetype = Data, color = Method)) + 
  geom_smooth(se = F, linewidth = 1) + geom_hline(aes(yintercept = alpha), linetype = 4) +
  facet_grid(.~alpha) + theme_bw() + ylab("Type I Error Rate") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.position = "none") + ylim(c(0,0.2))
p2 = ggplot(plt_type1e %>% filter(alpha == 0.05), aes(x = bias, y = type1e, linetype = Data, color = Method)) + 
  geom_smooth(se = F, linewidth = 1) + geom_hline(aes(yintercept = alpha), linetype = 4) +
  facet_grid(.~alpha) + theme_bw() + ylab("Type I Error Rate") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.position = "none") + ylim(c(0,0.4))

p3 = ggplot(plt_power, aes(x = bias, y = power, linetype = Data, color = Method)) + geom_smooth(se = F, linewidth = 1) + 
  theme_bw() + ylab("Power") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12) ) + ylim(c(0.4,1))
p3

pdf("figures/two_binary_example_multi.pdf", width = 14, height = 4)
gridExtra::grid.arrange(p1, p2, p3, widths = c(4,4,6), heights = c(4) )
dev.off()

saveRDS(list(type1e, type1e_wald, type1e_oracle, power, power_wald, power_oracle), "Robj/binary_multi.rds")


# modified S1 type of simulations
binary_marg_multi_m1 = function(n_pos, n_all){
  # browser()
  if(is.null(nrow(n_pos))){
    n_pos = matrix(n_pos, nrow = 1)
    n_all = matrix(n_all, nrow = 1)
  }
  n_c = n_all[,2, drop = F] + n_all[,3, drop = F]
  s_c = n_pos[,2, drop = F] + n_pos[,3, drop = F]
  
  mu_exp = (n_pos[,1]+1)/(n_all[,1]+2)
  mu_ctl = (s_c + 1)/(n_c + 2)
  
  var_exp = (n_pos[,1]+1)*(n_all[,1] - n_pos[,1] + 1)/(n_all[,1]+2)^2/(n_all[,1]+3)
  var_ctl = (s_c + 1)*(n_c - s_c + 1)/(n_c + 2)^2/(n_c+3)
  
  mu_diff = mu_exp - mu_ctl
  var_diff = var_exp + var_ctl
  
  -expm1(sum(pnorm(0, mean = mu_diff, sd = sqrt(var_diff), log.p = T)))
}

binary_marg_multi_m2 = function(n_pos, n_all){
  # browser()
  if(is.null(nrow(n_pos))){
    n_pos = matrix(n_pos, nrow = 1)
    n_all = matrix(n_all, nrow = 1)
  }
  n_c = n_all[,2, drop = F] + n_all[,3, drop = F]
  s_c = n_pos[,2, drop = F] + n_pos[,3, drop = F]
  
  mu_exp = (n_pos[,1]+1)/(n_all[,1]+2)
  mu_ctl = (s_c + 1)/(n_c + 2)
  
  var_exp = (n_pos[,1]+1)*(n_all[,1] - n_pos[,1] + 1)/(n_all[,1]+2)^2/(n_all[,1]+3)
  var_ctl = (s_c + 1)*(n_c - s_c + 1)/(n_c + 2)^2/(n_c+3)
  
  mu_diff = mu_exp - mu_ctl
  var_diff = var_exp + var_ctl
  
  # regret with normal approx
  rho = rowSums(n_all[,-ncol(n_all)])
  reg_ind = c((mu_diff + dnorm(mu_diff/sqrt(var_diff))/pnorm(mu_diff/sqrt(var_diff))*sqrt(var_diff))*pnorm(mu_diff/sqrt(var_diff)))
  sum(rho*reg_ind)/sum(rho)
}

library(snowfall)
sfInit(parallel = T, cpus = 8)
sfExport("binary_sim_multi", "binary_marg_multi",
         "binary_marg_multi_m1", "binary_marg_multi_m2", "permutation_binary_multi")
sfLibrary(dqrng)

all_res = lapply(seq(-0.3,0.3,0.05), function(x){
  print(x)
  sfSapply(1:1e3, function(rep){
    sim_res = binary_sim_multi(n_tot = c(75, 75), r_all = c(1/2, 1/2), p_ctrl_all = c(0.5, 0.5),
                               p_delta_all = c(0,x), n_e = c(375, 375), p_e = c(0.5, 0.5))
    
    p_org = permutation_binary_multi(sim_res$df, e_pos = sim_res$e_pos, e_all = sim_res$e_all, n_perm = 1e4)
    p_m1 = permutation_binary_multi(sim_res$df, e_pos = sim_res$e_pos, e_all = sim_res$e_all, 
                                    m = binary_marg_multi_m1, n_perm = 1e4)
    p_m2 = permutation_binary_multi(sim_res$df, e_pos = sim_res$e_pos, e_all = sim_res$e_all, 
                                    m = binary_marg_multi_m2, n_perm = 1e4)
    c(p_org, p_m1, p_m2)
  })
})

plt_data = data.frame(power = c(sapply(all_res, function(x) rowMeans(x<0.05))),
                      method = factor(c("EP-PT-m", "ED-PT-m1", "ED-PT-m2"), 
                                      levels = c("EP-PT-m", "ED-PT-m1", "ED-PT-m2")),
                      gamma1 = rep(seq(-0.3,0.3,0.05), each = 3))
library(ggplot2)
saveRDS(all_res, file = "binary_multiple_modified_S1.rds")

p1 = ggplot(plt_data, aes(x = gamma1, y = power, color = method, linetype = method)) +
  geom_smooth(se = F, linewidth = 0.5, method = "gam") +
  geom_hline(yintercept = 0.05, linetype = 4) + theme_bw() +
  xlab(expression(gamma[1])) + ylab("Frequency of rejection") + 
  scale_color_manual(values = c("black", "#F8766D", "#619CFF")) +
  scale_linetype_manual(values = c(2, 1, 1)) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, color = "black"))
ggsave("figures/modified_m_rejection_prob_updated_binary_multi.pdf", p1, width = 6, height = 4)
