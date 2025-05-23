library(dplyr)

# two groups
library(snowfall)
sfInit(parallel = T, cpus = 30)
sfSource("example_2_helper.R")
sfLibrary(dplyr)
system.time(type1e_sg <- lapply(seq(-0.1, 0.1, 0.025), function(mu0){
  print(mu0)
  sfSapply(1:1e4, function(rep){
    data = data_sim(n_all = c(50,50), r_all = c(1/2, 1/2), re_all = c(7.5, 7.5), 
                    mu0 = 0, gamma = 0, gamma_1 = 0, beta_1 = 0.5, mu0e = mu0)
    permute_test_sg(data, 10, 1e3)
})
  }))

power_sg = lapply(seq(-0.1, 0.1, 0.025), function(mu0){
  print(mu0)
  sfSapply(1:1e4, function(rep){
    data = data_sim(n_all = c(50,50), r_all = c(1/2, 1/2), re_all = c(7.5, 7.5), 
                    mu0 = 0, gamma = 0.5, gamma_1 = -0.2, beta_1 = 0.5, mu0e = mu0)
    permute_test_sg(data, 10, 1e3)
  })
})

res_all = list(type1e_sg, power_sg, type1e_cont, power_cont)
res_all = readRDS("example_2.rds")
type1e_sg = res_all[[1]]
power_sg = res_all[[2]]

res_type1e1 = sapply(type1e_sg, function(x){
  rowMeans(x<0.01)
})
res_type1e1[2,] = mean(res_type1e1[2,])
res_type1e1[3,] = mean(res_type1e1[3,])
res_type1e1[5,] = mean(res_type1e1[5,])

res_type1e5 = sapply(type1e_sg, function(x){
  rowMeans(x<0.05)
})
res_type1e5[2,] = mean(res_type1e5[2,])
res_type1e5[3,] = mean(res_type1e5[3,])
res_type1e5[5,] = mean(res_type1e5[5,])

res_power = sapply(power_sg, function(x){
  rowMeans(x<0.05)
})
res_power[2,] = mean(res_power[2,])
res_power[3,] = mean(res_power[3,])
res_power[5,] = mean(res_power[5,])

# plot for the example
df_type1e = data.frame(type1 = c(c(res_type1e1), c(res_type1e5)),
                       Method = c("Permutation", "Permutation", "Wald", "Wald", "Oracle"),
                       Data = c("ID + ED", "ID", "ID", "ID + ED", "ID + ED"),
                       beta_0 = rep(seq(-0.1,0.1,0.025), each = 5),
                       alpha = rep(c("alpha-level = 0.01", "alpha-level = 0.05"), each = length(res_type1e1)),
                       level = rep(c(0.01, 0.05), each = length(res_type1e1)))
df_type1e = df_type1e %>% mutate(Data = factor(Data, levels = c("ID + ED", "ID")),
                                 Method = factor(Method, levels = c("Oracle", "Permutation", "Wald")))
library(ggplot2)
p1 = ggplot(df_type1e %>% filter(level == 0.05), aes(x = beta_0, y = type1, linetype = Data, color = Method)) +
  geom_line() + geom_point() + geom_hline(aes(yintercept = level), linetype = 4) +
  facet_grid(.~alpha) + theme_bw() + ylab("Type I Error Rate") + xlab(expression(beta[0][","][0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.position = "none") + ylim(c(0,0.12))

df_power = data.frame(Power = c(res_power),
                      Method = c("Permutation", "Permutation", "Wald", "Wald", "Oracle"),
                      Data = c("ID + ED", "ID", "ID", "ID + ED", "ID + ED"),
                      beta_0 = rep(seq(-0.1,0.1,0.025), each = 5))
df_power = df_power %>% mutate(Data = factor(Data, levels = c("ID + ED", "ID")),
                               Method = factor(Method, levels = c("Oracle", "Permutation", "Wald")))
p2 = ggplot(df_power, aes(x = beta_0, y = Power, linetype = Data, color = Method)) + geom_line() + geom_point() +
  theme_bw() + ylab("Power") + xlab(expression(beta[0][","][0])) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12) ) +
  ylim(c(0.4,1))

# now investigate the accuracy of the asymptotic results
F_func = function(t, v, r, re, a, b, rho){
  # v is a matrix
  # K*n_mc dimension
  K = length(a)
  scale = sqrt(r*(1+re+r)/re)
  c0 = v + sqrt(rho*re/(1+re+r))*((r+re)/re*a - b)
  noncenter = colSums((c0/scale)^2)
  tnew = t/scale^2
  
  pchisq(tnew, df = K, ncp = noncenter)
}

Ftilde_inv_func = function(alpha, v, r, re, a, b, rho){
  # browser()
  # v is a matrix
  # K*n_mc dimension
  K = length(a)
  scale = sqrt(r*(1+re+r)/re)
  c0 = v + sqrt(rho*re/(1+re+r))*(a/(1+r) - b)
  noncenter = colSums((c0/scale)^2)
  q_a = qchisq(alpha, df = K, ncp = noncenter)
  q_a*scale^2
}

g_norm = function(r, re, a, b, rho, alpha = 0.05, n_mc = 1e4){
  K = length(a)
  v_all = t(matrix(rnorm(K*n_mc), ncol = K))
  s = F_func(Ftilde_inv_func(1-alpha, v_all, r, re, a, b, rho), v_all, r, re, a, b, rho)
  1 - mean(s)
}

# three curves, different b
sfExport("g_norm", "Ftilde_inv_func", "F_func")

# own simulation funcs
sfSource("example_2_helper.R")
sfLibrary(dplyr)

power_asym_test = lapply(seq(0,10,0.5), function(a){
  print(a)
  sfSapply(1:1e4, function(rep){
    data <- sim_asym(n = 1e3, r = 1/2, re = 7.5, p_all = c(0.5,0.5), 
                     a_all = c(a/sqrt(1e3),0), b_all = c(0,0), beta = c(1,1))
    p_perm = perm_asym_sim(data, 10, 1e3)
    
    ID = data$ID
    data$ID$group = as.factor(data$ID$group)
    m1_id = lm(y ~ group + a + group:a, data = ID)
    S_id = solve(summary(m1_id)$cov.unscaled[3:4, 3:4])
    chi2_id = coef(m1_id)[3:4]%*%(S_id%*%coef(m1_id)[3:4])
    p_f_id = 1 - pchisq(chi2_id, df = 2)
    c(p_perm, p_f_id)
  })
})

power_asym_test_1 = lapply(seq(0,10,0.5), function(a){
  print(a)
  sfSapply(1:1e4, function(rep){
    data <- sim_asym(n = 1e3, r = 1, re = 7.5, p_all = c(0.5,0.5), 
                     a_all = c(a/sqrt(1e3),0), b_all = c(0,0), beta = c(1,1))
    p_perm = perm_asym_sim(data, 10, 1e3)
    
    ID = data$ID
    data$ID$group = as.factor(data$ID$group)
    m1_id = lm(y ~ group + a + group:a, data = ID)
    S_id = solve(summary(m1_id)$cov.unscaled[3:4, 3:4])
    chi2_id = coef(m1_id)[3:4]%*%(S_id%*%coef(m1_id)[3:4])
    p_f_id = 1 - pchisq(chi2_id, df = 2)
    c(p_perm, p_f_id)
  })
})

power_asym_test_2 = lapply(seq(0,10,0.5), function(a){
  print(a)
  sfSapply(1:1e4, function(rep){
    data <- sim_asym(n = 1e3, r = 2, re = 7.5, p_all = c(0.5,0.5), 
                     a_all = c(a/sqrt(1e3),0), b_all = c(0,0), beta = c(1,1))
    p_perm = perm_asym_sim(data, 10, 1e3)
    
    ID = data$ID
    data$ID$group = as.factor(data$ID$group)
    m1_id = lm(y ~ group + a + group:a, data = ID)
    S_id = solve(summary(m1_id)$cov.unscaled[3:4, 3:4])
    chi2_id = coef(m1_id)[3:4]%*%(S_id%*%coef(m1_id)[3:4])
    p_f_id = 1 - pchisq(chi2_id, df = 2)
    c(p_perm, p_f_id)
  })
})

power_simulations = t(cbind(sapply(power_asym_test, function(x) rowMeans(x<0.05)), 
                          sapply(power_asym_test_2, function(x) rowMeans(x<0.05)),
                          sapply(power_asym_test_1, function(x) rowMeans(x<0.05))))
simp_all = data.frame(r = rep(c(1/2, 2, 1), each = 21),
                      a1 = seq(0,10,0.5),
                      permute = power_simulations[,1],
                      wald = power_simulations[,2])

# theoretical results
all_params = expand.grid( a1 = seq(0,10,0.5), r = c(0.5, 1, 2))
ap_all = sfApply(all_params, 1, function(x){
  g_norm(r = x[2], re = 7.5, a = c(x[1]*sqrt(1+x[2]),0), b = c(0, 0), rho = c(0.5,0.5))
})
ap_plt = cbind(all_params, ap_all)

asym_plt_all = left_join(ap_plt, simp_all) %>% pivot_longer(cols = ap_all:wald, names_to = "Test", values_to = "Power")
asym_plt_all$r = as.factor(asym_plt_all$r)

p_asym_pwr = ggplot(asym_plt_all, aes(x = a1, y = Power, color = r)) + 
  geom_line(aes(linetype = "Large Sample"), data = filter(asym_plt_all, Test =="permute")) + 
  geom_line(aes(linetype = "Test-B"), data = filter(asym_plt_all, Test =="wald")) +
  geom_point(aes(size = "Exact"), data = filter(asym_plt_all, Test == "ap_all")) + 
  theme_bw() +
  scale_color_jco() +
  scale_size_manual(values = c("Exact" = 1.5, "Large Sample" = .5, "Test-B" = .5), 
                    limits = c("Exact", "Large Sample", "Test-B")) +
  scale_linetype_manual(values = c("Exact" = "solid", "Large Sample" = "solid", "Test-B" = "dashed"), 
                        limits = c("Exact", "Large Sample", "Test-B")) +
  labs(size = NULL, linetype = NULL) +
  guides(
    size = guide_legend(
      override.aes = list(linetype = c("blank", "solid", "dashed"), shape = c(16, NA, NA))
    )
  ) +
  xlab(expression(a[1])) + ylab("Limiting Power") +
  theme(legend.position = c(0.25,0.65), legend.background = element_blank(), axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12))

pdf("figures/two_normal_example_all_updated.pdf", width = 14, height = 4)
gridExtra::grid.arrange(p1, p2, p_asym_pwr, widths = c(4,5.5,4.5), heights = c(4))
dev.off()
