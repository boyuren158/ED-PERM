system.time(type1e_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, out = "ext")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    #    permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    #    other_test(data_all) )
    other_test(data_all)
  } )
}))

sapply(type1e_other, function(x) rowMeans(x<0.05))

system.time(type1e_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, out = "self")
    # c( permute_test(data_all, prior_sd = 10, n_perm = 100, prob = F),
    # permute_test(data_all, prior_sd = 10, n_perm = 100, prob = T),
    other_test(data_all) #)
  } )
}))

system.time(pwr_global_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = (1:4)/2, out = "ext")
    other_test(data_all)
  })
}))

system.time(pwr_global_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = (1:4)/2, out = "self")
    other_test(data_all)
  })
}))

system.time(pwr_single_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(0, 0, 3, 0), out = "ext")
    other_test(data_all)
  })
}))

system.time(pwr_single_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(0, 0, 3, 0), out = "self")
    other_test(data_all)
  })
}))

system.time(pwr_mgmt_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(5, 5, 0, 0), out = "ext")
    other_test(data_all)
  })
}))

system.time(pwr_mgmt_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(5, 5, 0, 0), out = "self")
    other_test(data_all)
  })
}))

system.time(pwr_kps_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(2, 0, 2, 0), out = "ext")
    other_test(data_all)
  })
}))

system.time(pwr_kps_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(2, 0, 2, 0), out = "self")
    other_test(data_all)
  })
}))

system.time(pwr_neg_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(-2,0,0,0), out = "ext")
    other_test(data_all)
  })
}))

sapply(pwr_neg_other, function(x) rowMeans(x<0.05))

system.time(pwr_neg_self_other <- lapply(seq(50,250,50), function(n_out){
  print(n_out)
  sfSapply( 1:1e4, function(rep){
    data_all = data_generate_module(data_in_use = data_in_use, 
                                    data_out_use = data_out_use,
                                    n_in = 150, n_out = n_out, ratio = 2/3, alt = T, 
                                    lor = c(-2,0,0,0), out = "self")
    other_test(data_all)
  })
}))

# lrt threshold
lrt_thresh = quantile(type1e_other[[1]][1,], prob = 0.05)

test = lapply(list(pwr_global_other, pwr_global_self_other, pwr_single_other, pwr_single_self_other,
     pwr_mgmt_other, pwr_mgmt_self_other, pwr_kps_other, pwr_kps_self_other), function(ind){
       raw = sapply(ind, function(ind_n) rowMeans(ind_n < c(rep(lrt_thresh, 2), 0.05, 0.05)))
       raw[1,] = mean(raw[1,])
       raw[3,] = mean(raw[3,])
       raw
     })
other_pwr = data.frame( power_new = unlist(lapply(test, c)),
                        method = c("LR", "LR", "Z", "Z"),
                        source = c("ID", "ID + ED"),
                        ne = rep(seq(50,250,50), each = 4),
                        scenario = rep(c("Global", "Single", "MGMT", "KPS"), each = 40),
                        ED = rep(c("DFCI", "AVAGLIO"), each = 20) )

saveRDS(other_pwr, "Robj/other_tests.rds")

#negative other results
test = lapply(list(pwr_neg_other, pwr_neg_self_other), function(ind){
                     raw = sapply(ind, function(ind_n) rowMeans(ind_n < c(rep(lrt_thresh, 2), 0.05, 0.05)))
                     raw[1,] = mean(raw[1,])
                     raw[3,] = mean(raw[3,])
                     raw
                   })
neg_pwr = data.frame( power = unlist(lapply(test, c)),
                      method = c("Test-LR and LR-ED", "Test-LR and LR-ED", "Test-B and C", "Test-B and C"),
                      source = c("ID", "ID + ED"),
                      ne = rep(seq(50,250,50), each = 4),
                      scenario = "6. Negative (-2,0,0,0)",
                      ED = rep(c("DFCI", "AVAGLIO"), each = 20)) 
