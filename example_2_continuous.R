library(snowfall)
sfInit(cpus = 8, parallel = T)
sfSource("example_2_helper.R")
sfLibrary(dplyr)

system.time(power_cont <- lapply(c(1,seq(5,40,5)), function(p){
  print(p)
  sfSapply(1:1e4, function(rep){
    data = data_sim_continuous_subg(n = 100, r = 1/2, re = 7.5, p = p, rho = c(0.5,0.5),
                               mu0 = 0, gamma = 0.6, gamma_1 = c(-0.2, rep(0,p)), beta_1 = c(0.5, rep(1/sqrt(p), p)), mu0e = 0)
    permute_test_continuous_subg(data, 10, 1e3)
  })
}))

system.time(power_cont_mis <- lapply(c(1), function(p){
  print(p)
  sfSapply(1:1e4, function(rep){
    data = data_sim_continuous_subg(n = 100, r = 1/2, re = 7.5, p = p, rho = c(0.5,0.5),
                                    mu0 = 0, gamma = 0.6, gamma_1 = c(-0.2, rep(0,p)), 
                                    beta_1 = c(0.5, rep(1/sqrt(p), p)), mu0e = 0.2)
    permute_test_continuous_subg(data, 10, 1e3)
  })
}))

z_test_cont <- lapply(c(1,seq(5,40,5)), function(p){
  sapply(1:1e4, function(rep){
    data = data_sim_continuous_subg(n = 100, r = 1/2, re = 7.5, p = p, rho = c(0.5,0.5),
                                  mu0 = 0, gamma = 0.6, gamma_1 = c(-0.2, rep(0,p)), beta_1 = c(0.5, rep(1/sqrt(p), p)), mu0e = 0)
    res = lm( y ~ a, data = data$D_in )
    summary(res)$coefficients[2,4]
  })
})

system.time(power_cont_single <- lapply(c(1,seq(5,40,5)), function(p){
  print(p)
  sfSapply(1:1e4, function(rep){
    data = data_sim_continuous_subg(n = 100, r = 1/2, re = 7.5, p = p, rho = c(0.5,0.5),
                                    mu0 = 0, gamma = 0.75, gamma_1 = c(-0.75, rep(0,p)), beta_1 = c(0.5, rep(1/sqrt(p), p)), mu0e = 0)
    permute_test_continuous_subg(data, 10, 1e3)
  })
}))

z_test_cont_single <- lapply(c(1, seq(5,40,5)), function(p){
  sapply(1:1e4, function(rep){
    data = data_sim_continuous_subg(n = 100, r = 1/2, re = 7.5, p = p, rho = c(0.5,0.5),
                                    mu0 = 0, gamma = 0.75, gamma_1 = c(-0.75, rep(0,p)), beta_1 = c(0.5, rep(1/sqrt(p), p)), mu0e = 0)
    res = lm( y ~ a, data = data$D_in )
    summary(res)$coefficients[2,4]
  })
})

res_power_both = sapply(power_cont, function(x){
  rowMeans(x<0.05)
})

res_power_both[4,] = sapply(z_test_cont, function(x) mean(x<0.05))

res_power_single = sapply(power_cont_single, function(x){
  rowMeans(x<0.05)
})
res_power_single[4,] = sapply(z_test_cont_single, function(x) mean(x<0.05))

plt_data = data.frame(power = c(res_power_both[1,],
                                res_power_both[3,],
                                res_power_both[4,]),
                      Method = rep(c("ED-PT", "Test-B", "Z-test"), each = ncol(res_power_both)),
                      d = c(1, seq(5,40,5))+1)

plt_data_s = data.frame(power = c(res_power_single[1,],
                                res_power_single[3,],
                                res_power_single[4,]),
                      Method = rep(c("ED-PT", "Test-B", "Z-test"), each = ncol(res_power_both)),
                      d = c(1, seq(5,40,5))+1)


library(ggplot2)
plt_data = plt_data %>% mutate(Data = factor(ifelse(Method == "ED-PT", "ID + ED", "ID"), levels = c("ID + ED", "ID")))
plt_data_s = plt_data_s %>% mutate(Data = factor(ifelse(Method == "ED-PT", "ID + ED", "ID"), levels = c("ID + ED", "ID")))
p1 = ggplot(data = plt_data, aes(x = d, y = power, color = Method, linetype = Data)) + geom_line() + geom_point() +
  scale_color_manual(values = c("#00BA38", "#619CFF", "#C77CFF")) +
  ylim(c(0.4,1)) + theme_bw() + ylab("Power") + ggtitle("Treatment effects = (0.6, 0.4)") +
  theme( title = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12), legend.position = "none")
p2 = ggplot(data = plt_data_s, aes(x = d, y = power, color = Method, linetype = Data)) + geom_line() + geom_point() +
  scale_color_manual(values = c("#00BA38", "#619CFF", "#C77CFF")) +
  ylim(c(0.2,1)) + theme_bw() + ylab("Power") + ggtitle("Treatment effects = (0.75, 0)") +
  theme(title = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12))

saveRDS(list(res1 = plt_data, res2 = plt_data_s), "example_2_continuous.rds")

pdf("figures/example2_continuous_new.pdf", width = 9, height = 4)
gridExtra::grid.arrange(p1, p2, widths = c(3.9,5.1), heights = c(4))
dev.off()
