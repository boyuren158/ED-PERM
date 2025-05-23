library(dplyr)
library(survival)
library(gtools)
set.seed(1)
source("permutation_src.R")

# load data
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

cbind( coef(glm( outcome ~ Age + Sex + eor + kps90_100*mgmt, data = data_in_use)),
       coef(glm( outcome ~ Age + Sex + eor + kps90_100*mgmt, data = data_out_use)) )
# seems like eor and kps:mgmt has the highest difference

### odds ratio testing code
n_out = 100
p_all = sapply(1:1e3, function(rep){
  data1 = data_generate_module(data_in_use = data_in_use, data_out_use = data_out_use,
                               n_in = 500, n_out = n_out, ratio = 2/3, alt = T, 
                               lor = c(1,0,0,0), out = "ext")
  data1$data_in$subtype = as.factor(data1$data_in$subtype)
  res = glm( outcome ~ Age + Sex + eor + subtype*casectrl, family = "binomial", data = data1$data_in)
  res_null = glm( outcome ~ Age + Sex + eor + subtype, family = "binomial", data = data1$data_in)
  anova(res, res_null, test = "Chisq")$`Pr(>Chi)`[2]
})
mean(p_all < 0.05)

###

library(snowfall)
sfInit(parallel = T, cpus = 8)
sfExport("data_in_use", "data_out_use")
sfSource("permutation_src.R")
sfLibrary(dplyr)
sfLibrary(gtools)

system.time(type1e <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, out = "ext")
    c( permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "marginal"),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "prob"),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "regret") )
    } )
}))

sapply(type1e, function(x) rowMeans(x<0.05))

system.time(type1e_in <- lapply(seq(5,5,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, out = "ext")
    c( permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "marginal", ID = T),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "prob", ID = T),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "regret", ID = T) )
  } )
}))

lapply(type1e_in, function(x) rowMeans(x<0.05))

df1 = data.frame(type1e = c(c(sapply(type1e, function(x) rowMeans(x<0.05))), 
                 rep(c(sapply(type1e_in, function(x) rowMeans(x<0.05))),5)),
                 method = c("ED-PT-m", "ED-PT-m1", "ED-PT-m2"),
                 ne = rep(seq(50,250,50), each = 3),
                 source = rep(c("ID + ED", "ID"), each = 15))
type1e_lrz = sapply(type1e_other, function(x) rowMeans(x<0.05))
type1e_lrz[1,] = mean(type1e_lrz[1,])
type1e_lrz[3,] = mean(type1e_lrz[3,])

df2 = data.frame(type1e = c(type1e_lrz),
                 method = rep(c("Test-LR and LR-ED", 
                                "Test-B and C"), each = 2),
                 source = c("ID", "ID + ED"),
                 ne = rep(seq(50,250,50), each = 4))

df3 = data.frame(type1e = c(sapply(p_all_dfci, function(x) mean(x<0.05)),
                            colMeans(p_all_dfci_ipw<0.05)),
                 method = rep(c("Test-Matching", "Test-IPW"), each = 5),
                 ne = seq(50,250,50),
                 source = "ID + ED")

df = rbind(df1, df2, df3)
df = df %>% mutate(method = factor(method, levels = c("ED-PT-m", "ED-PT-m1",
                                                      "ED-PT-m2", "Test-LR and LR-ED", "Test-B and C",
                                                      "Test-Matching", "Test-IPW")),
                   source = factor(source, levels = c("ID + ED", "ID")))

saveRDS(df, "Robj/type1e_plt.rds")

# scale_color_manual(name = "Method", labels = c("ED-PT", "Test-A", "Test-LR", "Test-B", "ED-PT-Inf"),
#                    values = c("#F8766D", "#F8766D", "#00BA38", "#619CFF","#F8766D"))

library(ggplot2)
p1 = ggplot(df, aes(x = ne, y = type1e, color = method, linetype = source)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0.05, linetype = 4) +
  theme_bw() + ylab("Type I Error Rate") + xlab(expression(n[E])) +
  scale_color_manual(values = c(
    "#F8766D", "#CC3366", "#660066", "#00BA38",
    "#619CFF", "#FF9933", "#FFCC00")
  ) +
  theme(legend.position = c(0.2,0.9), legend.background = element_blank()) +
  theme(axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12))
ggsave("figures/application_example_type1e.pdf", p1, width = 6, height = 6)

system.time(type1e_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e5, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, out = "self")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
       # permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
       other_test(data_all) #)
  } )
}))

sapply(type1e_self, function(x) rowMeans(x<0.05))
sapply(type1e_self, function(x){
  rowMeans(x<0.05)
})

lrt_thresh = quantile( c(sapply(type1e, function(x) x[1,])), probs = 0.05 )

system.time(pwr_negative_id <- lapply(seq(5,5,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(-2,0,0,0), out = "ext")
    c( permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "marginal", ID = T),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "prob", ID = T),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "regret", ID = T) )
  })
}))

sapply(pwr_negative, function(x) mean(x<0.05))

system.time(pwr_negative_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(-2,0,0,0), out = "self")
    c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
       permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
       other_test(data_all) )
  })
}))

sapply(pwr_negative_self, function(x) rowMeans(x<0.05))

system.time(pwr_global <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = (1:4)/2, out = "ext")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
    permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "regret")
    })
}))

sapply(pwr_global, function(x) mean(x<0.05))

system.time(pwr_global_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = (1:4)/2, out = "self")
    permute_test(data_all, prior_sd = 10, self = T, n_perm = 1e2, prob = T)
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
  })
}))

sapply(pwr_global_self, function(x) rowMeans(x<0.05))

system.time(pwr_single <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(0, 0, 3, 0), out = "ext")
    c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
       permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
       other_test(data_all) )
  })
}))

system.time(pwr_single_self <- lapply(seq(50,50,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(0, 0, 3, 0), out = "self")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = F),
    c( permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = T, ID = T, threshold = 5),
       permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = T, self = T, threshold = 5) )#,
       # other_test(data_all) )
  })
}))

data_single = lapply(1:1e3, function(rep){
  data_generate_module(data_in_use = data_in_use, 
                       data_out_use = data_out_use,
                       n_in = 150, n_out = 250, ratio = 2/3, alt = T, 
                       lor = c(0, 0, 3, 0), out = "ext")
})
sfExport("data_single")

system.time(pwr_single_self_regret <- lapply(seq(50,250,50), function(ne){
  print(ne)
  sfSapply( 1:1e3, function(rep){
    permute_test(list(data_in = data_single[[rep]]$data_in,
                      data_out = data_single[[rep]]$data_out[1:ne,]), 
                 prior_sd = 10, n_perm = 1e3, prob = "regret", threshold = 2)
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
  })
}))

sapply(pwr_single_self_regret, function(x) mean(x<0.05))        
              
mean(pwr_single_self_in < 0.05)

system.time(pwr_single_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_raw = data_single_self[[rep]]
    ED = data_raw$data_out[1:n_out,]
    permute_test(list(data_in = data_raw$data_in, data_out = ED), prior_sd = 10,
                 self = T, n_perm = 1e2, prob = T, ID = F)
  })
}))

# was 0.621 0.612 0.612 0.608 0.589
mean(pwr_single_self_in < 0.05)
sapply(pwr_single_self, function(x) mean(x<0.05))

plot( pwr_single_self[[1]][2,], pwr_single_self[[2]][2,] )
abline(a = 0, b = 1)
abline(h = 0.05, v = 0.05)
which(( pwr_single_self[[1]][2,] < pwr_single_self[[2]][2,] ))


system.time(pwr_mgmt <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(5, 5, 0, 0), out = "ext")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
    permute_test(data_all, prior_sd = 10, n_perm = 1e3, prob = "regret")
  })
}))

sapply(pwr_mgmt, function(x) rowMeans(x<0.05))

system.time(pwr_mgmt_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(5, 5, 0, 0), out = "self")
    permute_test(data_all, prior_sd = 10, self = T, n_perm = 1e2, prob = T)
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
  })
}))

rowMeans(pwr_mgmt_self[[1]] < 0.05)

system.time(pwr_kps <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(2, 0, 2, 0), out = "ext")
    c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
       permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
       other_test(data_all) )
  })
}))
sapply(pwr_kps, function(x) rowMeans(x<0.05))

system.time(pwr_kps_self <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e3, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(2, 0, 2, 0), out = "self")
    c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
       permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
       other_test(data_all) )
  })
}))