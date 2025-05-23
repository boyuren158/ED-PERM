# binary case
# generic test
binary_marg = function(n_pos, n_all){
  lgamma(n_pos[1]+1) + lgamma(n_all[1] - n_pos[1] + 1) +
    lgamma(n_pos[2]+n_pos[3]+1) + lgamma(n_all[2] + n_all[3] - n_pos[2] - n_pos[3]+1)
}

# n_pos: # of success in exp, in control and in ED
# n_all: sample size of exp, control and ED
permutation_binary = function(n_pos, n_all, marg_fn = binary_marg, exact = T, ...){
  # browser()
  # r is the total success in ID
  r = n_pos[1] + n_pos[2]
  # compute log marginal likelihood
  m_obs = marg_fn(n_pos, n_all, ...)
  # compute all possible marginal likelhiood after permutation
  # note the range of s1 is max(0, s1 + s2 - n2), min(s1+s2, n1)
  s1_val_all = max(0, r - n_all[2]):min(r, n_all[1])
  m_perm_all = sapply(s1_val_all, function(x){
    marg_fn(c(x, r-x, n_pos[3]), n_all, ...)
  })
  s1_lower = s1_val_all[m_perm_all<m_obs]
  s1_eq = s1_val_all[m_perm_all==m_obs]
  # 
  # # compute the p-value with hypergeometric distribution
  if( exact ){
    1 - sum(dhyper(s1_lower, m = n_all[1], n = n_all[2],  k = r)) - runif(1)*sum(dhyper(s1_eq, m = n_all[1], n = n_all[2],  k = r))
  }else{
    1 - sum(dhyper(s1_lower, m = n_all[1], n = n_all[2],  k = r))
  }
}

# binary case m1 and m2 variants of m(D)
m1_binary = function(n_pos, n_all){
  # theta1 = rbeta(n_rep, shape1 = n_pos[1] + 1, shape2 = n_all[1] - n_pos[1] + 1)
  # theta0 = rbeta(n_rep, shape1 = n_pos[2] + n_pos[3] + 1, shape2 = n_all[2] + n_all[3] - n_pos[2] - n_pos[3] + 1)
  # mean(theta1 > theta0)
  # normal approximation
  alpha1 = n_pos[1] + 1
  beta1 = n_all[1] - n_pos[1] + 1
  alpha2 = n_pos[2] + n_pos[3] + 1
  beta2 = n_all[2] + n_all[3] - alpha2 + 2
  mu_1 = alpha1/(alpha1 + beta1)
  var_1 = alpha1*beta1/(alpha1 + beta1)^2/(alpha1 + beta1)
  mu_2 = alpha2/(alpha2 + beta2)
  var_2 = alpha2*beta2/(alpha2 + beta2)^2/(alpha2 + beta2)
  
  pnorm(0, mean = mu_1 - mu_2, sd = sqrt(var_1 + var_2), lower.tail = F)
}

m2_binary = function(n_pos, n_all){
  alpha1 = n_pos[1] + 1
  beta1 = n_all[1] - n_pos[1] + 1
  alpha2 = n_pos[2] + n_pos[3] + 1
  beta2 = n_all[2] + n_all[3] - alpha2 + 2
  mu_1 = alpha1/(alpha1 + beta1)
  var_1 = alpha1*beta1/(alpha1 + beta1)^2/(alpha1 + beta1)
  mu_2 = alpha2/(alpha2 + beta2)
  var_2 = alpha2*beta2/(alpha2 + beta2)^2/(alpha2 + beta2)
  
  mu = mu_1 - mu_2
  sigma = sqrt(var_1 + var_2)
  ptail = pnorm(0, mean = mu_1 - mu_2, sd = sqrt(var_1 + var_2), lower.tail = F)
  
  (mu + dnorm(mu/sigma)/pnorm(mu/sigma)*sigma)*ptail
}

p_approx = function(n_pos, r_all, n1){
  s = n_pos[1] + n_pos[2]#sum(n_pos)
  sprime = 2*(s+n_pos[3])/(sum(r_all) + 1)- n_pos[1]
  mu = s/(r_all[1]+1)
  # sd = sqrt((n_pos[1] + n_pos[2])*r_all[1]*(1-p0)/(r_all[1]+1)^2)
  var = s*r_all[1]*((r_all[1]+1)*n1-s)/((r_all[1]+1)*n1 - 1)/(r_all[1]+1)^2
  
  1 - pnorm((max(n_pos[1], sprime) - mu)/sqrt(var)) + pnorm((min(n_pos[1], sprime) - mu)/sqrt(var))
}

# plotting the approx p vs. true p
all_params = expand.grid(pi0 = seq(0.2,0.8, length = 10), 
                         a = seq(0,2,length = 10),
                         b = seq(-3,3, length = 10))

p_val_plot = as.data.frame(t(apply(all_params, 1, function(x){
  n1 = 1e4
  n_all = c(n1, n1/2, n1*5)
  n_pos = rbinom(3, size = n_all, prob = c(x[1] + x[2]/sqrt(n1), x[1], x[1] + x[3]/sqrt(n1)))
  p_e = permutation_binary(n_pos, n_all, exact = T)
  p_a = p_approx(n_pos, c(1/2, 5), n1)
  c(p_e, p_a)
})))

colnames(p_val_plot) = c("Actual", "Approximation")

library(ggplot2)
p_val_approx = ggplot(p_val_plot, aes(x = log10(Actual), y = log10(Approximation))) + geom_point() +
  theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = 2) + xlab(expression(log[10]("Actual p-value"))) + 
  ylab(expression(log[10]("Approximated p-value"))) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12))
p_val_approx

# rejection region
library(dplyr)
se = rbinom(n = 1, size = 5e4, prob = 0.5)
all_values = expand.grid(a = seq(-2,2,0.01), b = seq(-2,2,0.01))
p_region = apply(all_values, 1, function(x){
  p1 = 0.5 + x[1]/sqrt(1e4)
  p0 = 0.5 + x[2]/sqrt(5e3)
  p_approx(n_pos = c(round(p1*1e4), round(p0*5e3), se), r_all = c(1/2, 5), n1 = 1e4)
})
p_region_prop = apply(all_values, 1, function(x){
  p1 = 0.5 + x[1]/sqrt(1e4)
  p0 = 0.5 + x[2]/sqrt(5e3)
  prop.test(x = c(round(p1*1e4), round(p0*5e3)), n = c(1e4, 5e3), correct = F)$p.value
})
all_values$z = factor(ifelse(p_region<=0.05, "reject", "not reject"), level = c("reject", "not reject"))
all_values$z_prop = factor(ifelse(p_region_prop<=0.05, "reject", "not reject"), level = c("reject", "not reject"))
all_values = all_values %>% mutate( z_use = factor(case_when( z == "reject" & z_prop =="reject" ~ "Both",
                                                       z == "not reject" & z_prop =="not reject" ~ "Neither",
                                          z == "reject" & z_prop =="not reject" ~ "ED-PT only",
                                          z == "not reject" & z_prop =="reject" ~ "Test-B only"),
                                          levels = c('ED-PT only',
                                                     'Test-B only',
                                                     'Both',
                                                     'Neither')) )
library(ggpattern)
p_rej_region = ggplot(all_values %>% filter(!is.na(z_use)), aes(y = a, x = b, fill = as.factor(z_use))) + geom_tile(color = NA) +
  ylab(expression(sqrt(n[1])*(s[1]/n[1]-0.5))) + xlab(expression(sqrt(n[0])*(s[0]/n[0]-0.5))) + coord_fixed() +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() + scale_fill_manual(values = c("#F8766D", "#00BFC4", "#FFCCCC", "gray")) +
  guides(fill=guide_legend(title="")) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12))
  
p_rej_region

# comparing theoretical vs. approximation of limiting power function
limiting_power_true = function(r, re, a, b, p0){
  mu = c(sqrt(r)*a/sqrt((r+1)*p0*(1-p0)),
         (2*(r+1)*b*re - (r*(r + re + 1) + 2*re)*a)/sqrt(r*(r+1)*p0*(1-p0))/(r + re + 1))
  Sigma = matrix(c(1, -1, -1, 1+4*re/r/(r + re + 1)), 
                 nrow = 2)
  
  yall = MASS::mvrnorm(n = 1e5, mu = mu, Sigma = Sigma)
  mean(pnorm(pmax(yall[,1], yall[,2])) - pnorm(pmin(yall[,1], yall[,2])) > 0.95)
}

limiting_power_approx = function(r, re, a, b, p0){
  mu = sqrt(r)*a/sqrt((r+1)*p0*(1-p0)) - qnorm(0.95)
  pnorm(mu)
}

truth = sapply(c(1/2,1,2), function(r){
  sapply(seq(0,2,0.1), function(a){
    limiting_power_true(r, 5, a, b = 0, p0 = 0.5)
  })
})

approx = sapply(c(1/2,1,2), function(r){
  sapply(seq(0,2,0.1), function(a){
    limiting_power_approx(r, 5, a, b = 0, p0 = 0.5)
    })
})

ztest = sapply(c(1/2,1,2), function(r){
  sapply(seq(0,2,0.1), function(a){
    p0 = 0.5
    mu1 = sqrt(r)*a/sqrt((r+1)*p0*(1-p0)) - qnorm(0.975)
    mu2 = -sqrt(r)*a/sqrt((r+1)*p0*(1-p0)) - qnorm(0.975)
    pnorm(mu1) + pnorm(mu2)
  })
})

library(tidyr)
library(dplyr)
library(ggplot2)
asym_pwr = data.frame(a = seq(0,2,0.1),
                      truth = c(truth),
                      approx = c(approx),
                      ztest = c(ztest),
                      r = rep(c(1/2,1,2), each = 21)) %>% 
  pivot_longer(cols = truth:ztest, names_to = "type", values_to = "Power") %>%
  mutate(type = factor(type, levels = c("truth", "approx", "ztest"), labels = c("Actual", "Approximation", "Test-B")),
         r = as.factor(r))

library(ggsci)

p_asym_pwr = ggplot(asym_pwr, aes(x = a, y = Power, color = r)) + 
  geom_line(aes(linetype = "Approximation"), data = filter(asym_pwr, type =="Approximation")) + 
  geom_line(aes(linetype = "Test-B"), data = filter(asym_pwr, type =="Test-B")) +
  geom_point(aes(size = "Actual"), data = filter(asym_pwr, type == "Actual")) + 
  theme_bw() +
  scale_color_jco() +
  scale_size_manual(values = c("Actual" = 1.5, "Approximation" = .5, "Test-B" = .5), 
                    limits = c("Actual", "Approximation", "Test-B")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Approximation" = "solid", "Test-B" = "dashed"), 
                         limits = c("Actual", "Approximation", "Test-B")) +
  labs(size = NULL, linetype = NULL) +
  guides(
    size = guide_legend(
      override.aes = list(linetype = c("blank", "solid", "dashed"), shape = c(16, NA, NA))
    )
  ) +
  xlab(expression(a)) + ylab("Limiting Power") +
  theme(legend.position = c(0.3,0.65), legend.background = element_blank(), axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12))

pdf("figures/two_binary_example_asym_updated.pdf", width = 8, height = 4)
# gridExtra::grid.arrange(p_val_approx, p_asym_pwr, p_rej_region, widths = c(4,4.2,7), heights = 4)
gridExtra::grid.arrange(p_val_approx, p_asym_pwr, widths = c(4,4), heights = 4)
dev.off()
  
# simulations
example1_driver = function(n_all, p_all, n_rep){
  # browser()
  sapply(1:n_rep, function(rep){
    # generate data first
    # p_all = p1, p0, pe
    # n_all = n1, n0, ne
    n_pos = rbinom(3, size = n_all, prob = p_all)
    p_perm_in = permutation_binary(c(n_pos[1:2], 0), c(n_all[1:2], 0))
    p_perm_all = permutation_binary(n_pos, n_all)
    p_z_in = prop.test(x = n_pos[1:2], n = n_all[1:2], correct = F)$p.value
    p_z_all = prop.test(x = c(n_pos[1], n_pos[2] + n_pos[3]), n = c(n_all[1], n_all[2] + n_all[3]), correct = F)$p.value
    p_oracle = prop.test(x = n_pos[1], n = n_all[1], p = p_all[2], alternative = "two.sided")$p.value
    c(p_perm_in, p_perm_all, p_z_in, p_z_all, p_oracle)
  })
}

example1_driver_neg = function(n_all, p_all, n_rep, marg_fn, ...){
  # browser()
  sapply(1:n_rep, function(rep){
    # generate data first
    # p_all = p1, p0, pe
    # n_all = n1, n0, ne
    n_pos = rbinom(3, size = n_all, prob = p_all)
    p_perm_all = permutation_binary(n_pos, n_all, marg_fn = binary_marg, ...)
    p_m1_all = permutation_binary(n_pos, n_all, marg_fn = m1_binary, ...)
    p_m2_all = permutation_binary(n_pos, n_all, marg_fn = m2_binary, ...)
    c(p_perm_all, p_m1_all, p_m2_all)
  })
}

system.time(example1_type1e <- lapply(seq(-0.1, 0.1, 0.025), function(bias){
  example1_driver(n_all = c(66, 34, 1e2), p_all = c(0.5, 0.5, 0.5 + bias), n_rep = 1e4)}))

system.time(example1_power <- lapply(seq(-0.1, 0.1, 0.025), function(bias){
  example1_driver(n_all = c(66, 34, 5e2), p_all = c(0.75, 0.5, 0.5 + bias), n_rep = 5e4)}))

# negative treatment effects
system.time(example1_power_neg <- lapply(seq(-0.25, 0.25, 0.025), function(gamma){
  print(gamma)
  example1_driver_neg(n_all = c(66, 34, 5e2), p_all = c(0.5 + gamma, 0.5, 0.5), n_rep = 1e3)}))

# system.time(example1_power <- lapply(0, function(bias){
#   example1_driver(n_all = c(1e4, 5e3, 5e4), p_all = c(0.5+0.5/sqrt(1e4), 0.5, 0.5), n_rep = 10000)}))

res_type1e1 = sapply(example1_type1e, function(x){
  rowMeans(x<0.01)
})

res_type1e5 = sapply(example1_type1e, function(x){
  rowMeans(x<0.05)
})
res_type1e5[1,] = mean(res_type1e5[1,])
res_type1e5[3,] = mean(res_type1e5[3,])
res_type1e5[5,] = mean(res_type1e5[5,])

res_power = sapply(example1_power, function(x){
  rowMeans(x<0.05)
})
res_power[1,] = mean(res_power[1,])
res_power[3,] = mean(res_power[3,])
res_power[5,] = mean(res_power[5,])

# plot for the example
df_type1e = data.frame(type1 = c(c(res_type1e1), c(res_type1e5)),
                       Method = c("Permutation", "Permutation", "Wald", "Wald", "Oracle"),
                       Data = c("ID", "ID + ED", "ID", "ID + ED", "ID + ED"),
                       beta_0 = rep(seq(-0.1,0.1,0.025), each = 5),
                       alpha = rep(c("alpha-level = 0.01", "alpha-level = 0.05"), each = length(res_type1e1)),
                       level = rep(c(0.01, 0.05), each = length(res_type1e1)))
df_type1e = df_type1e %>% mutate(Data = factor(Data, levels = c("ID + ED", "ID")))
library(ggplot2)
p1 = ggplot(df_type1e %>% filter(level == 0.01), aes(x = beta_0, y = type1, linetype = Data, color = Method)) + geom_smooth(se = F, linewidth = 1) + 
  geom_hline(aes(yintercept = level), linetype = 4) +
  facet_grid(.~alpha) + theme_bw() + ylab("Type I Error Rate") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.position = "none") + ylim(c(0,0.2))
p2 = ggplot(df_type1e %>% filter(level == 0.05), aes(x = beta_0, y = type1, linetype = Data, color = Method)) + geom_smooth(se = F, linewidth = 1) + 
  geom_hline(aes(yintercept = level), linetype = 4) +
  facet_grid(.~alpha) + theme_bw() + ylab("Type I Error Rate") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.position = "none") + ylim(c(0,0.4))

df_power = data.frame(Power = c(res_power),
                      Method = c("Permutation", "Permutation", "Wald", "Wald", "Oracle"),
                      Data = c("ID", "ID + ED", "ID", "ID + ED", "ID + ED"),
                      beta_0 = rep(seq(-0.1,0.1,0.025), each = 5))
df_power = df_power %>% mutate(Data = factor(Data, levels = c("ID + ED", "ID")))
p3 = ggplot(df_power, aes(x = beta_0, y = Power, linetype = Data, color = Method)) + geom_smooth(se = F, linewidth = 1) + 
  theme_bw() + ylab("Power") + xlab(expression(beta[0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12) ) + ylim(c(0.6,1))
p3

pdf("figures/two_binary_example_all.pdf", width = 14, height = 8)
gridExtra::grid.arrange(p1, p2, p3, p_val_approx, p_asym_pwr, p_rej_region, widths = c(4,4,6), heights = c(4,4) )
dev.off()

# negative trt results
pwr_neg = sapply(example1_power_neg, function(x){
  rowMeans(x < 0.05)
})

plt_neg = data.frame( rej = c(pwr_neg[1,], pwr_neg[2,], pwr_neg[3,]),
                     gamma = seq(-0.25,0.25, 0.025),
                     method = rep(c("ED-PT-m", "ED-PT-m1", "ED-PT-m2"), each = 21))

p1 = ggplot(plt_neg, aes(x = gamma, y = rej, color = method, linetype = method)) + 
  geom_smooth(se = F, linewidth = 0.5, method = "gam") +
  geom_hline(yintercept = 0.05, linetype = 4) + theme_bw() +
  xlab(expression(gamma[1])) + ylab("Frequency of rejection") + 
  scale_color_manual(values = c("black", "#F8766D", "#619CFF")) +
  scale_linetype_manual(values = c(2, 1, 1)) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, color = "black"))
ggsave("figures/modified_m_rejection_prob_updated_binary.pdf", width = 6, height = 4)
