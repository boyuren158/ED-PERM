library(dplyr)
library(snowfall)
library(ggplot2)
sfInit(parallel = T, cpus = 8)

# two groups simulation, ED-PT-m and ED-PT-mtilde1
sfSource("example_2_helper.R")
sfLibrary(dplyr)

res2 <- sapply(seq(-1,1,length.out = 21), function(gamma2){
  print(gamma2)
  power_negative <- sfSapply(1:1e3, function(rep){
    data <- sim_asym(n = 50, r = 1/2, re = 15/2, p_all = c(0.5,0.5), 
                     a_all = c(0,gamma2), b_all = c(0,0), beta = c(0,0))
    p1 = perm_asym_sim_pos(data, 10, 1e3)
    p2 = perm_asym_sim(data, 10, 1e3)
    c(p1,p2)
  })
  rowMeans(power_negative < 0.05)
})

sfExport("dec_m", "dec_perm")
  res3 <- sapply(seq(-1,1,length.out = 21), function(gamma2){
  print(gamma2)
  power_negative <- sfSapply(1:1e3, function(rep){
    data <- sim_asym(n = 50, r = 1/2, re = 15/2, p_all = c(0.5,0.5), 
                     a_all = c(0,gamma2), b_all = c(0,0), beta = c(0,0))
    p1 = dec_perm(data, n_perm = 1e3)
    p2 = perm_asym_sim(data, 10, 1e3)
    c(p1,p2)
  })
  rowMeans(power_negative < 0.05)
})

# plotting
plt_df = data.frame( rej = c(res2[2,], res2[1,], res3[1,]),
                     gamma = seq(1,-1, length.out = 21),
                     method = rep(c("ED-PT-m", "ED-PT-m1", "ED-PT-m2"), each = 21))

p1 = ggplot(plt_df, aes(x = gamma, y = rej, color = method, linetype = method)) + 
  geom_smooth(se = F, linewidth = 0.5, method = "gam") +
  geom_hline(yintercept = 0.05, linetype = 4) + theme_bw() +
  xlab(expression(gamma[1])) + ylab("Frequency of rejection") + 
  scale_color_manual(values = c("black", "#F8766D", "#619CFF")) +
  scale_linetype_manual(values = c(2, 1, 1)) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, color = "black"))
ggsave("figures/modified_m_rejection_prob_updated_neg.pdf", width = 6, height = 4)

# additional setting
power_pp0 <- sfSapply(1:1e3, function(rep){
  data <- sim_asym(n = 75, r = 1/2, re = 15/2, p_all = c(1/3,1/3,1/3), 
                   a_all = c(0.75, 0.75, 0), b_all = c(0,0,0), beta = c(0,0,0))
  p1 = dec_perm(data, n_perm = 1e3)
  p2 = perm_asym_sim_pos(data, 10, 1e3)
  p3 = perm_asym_sim(data, 10, 1e3)
  c(p1,p2,p3)
})
rowMeans(power_pp0 < 0.05)

power_ppn <- sfSapply(1:1e3, function(rep){
  data <- sim_asym(n = 75, r = 1/2, re = 15/2, p_all = c(1/3,1/3,1/3), 
                   a_all = c(0.75, 0.75, -0.75), b_all = c(0,0,0), beta = c(0,0,0))
  p1 = dec_perm(data, n_perm = 1e3)
  p2 = perm_asym_sim_pos(data, 10, 1e3)
  p3 = perm_asym_sim(data, 10, 1e3)
  c(p1,p2,p3)
})
rowMeans(power_ppn < 0.05)

# negative effects
