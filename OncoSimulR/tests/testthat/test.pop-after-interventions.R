test_that("Population after one intervention Exp",{
  popSize <- 5000
  fractionPopSize <- 0.15
  sampleEvery <- 0.01
  modelChanges <- list( Change1 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(pop_after_interv <= fraction_pop) #all.equal(pop_after_interv, fraction_pop, tol=0)
  
})

test_that("Population after one intervention McFL",{
  popSize <- 5000
  fractionPopSize <- 0.15
  sampleEvery <- 0.01
  modelChanges <- list( Change1 = list( trigger = list(popSize = popSize),
                                        action = list(fractionPopSize = fractionPopSize)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "McFL", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(pop_after_interv <= fraction_pop)
  
})
#Test from here
test_that("Population after 3 repeated interventions Exp",{
  popSize <- 6000
  fractionPopSize <- 0.10
  sampleEvery <- 0.01
  
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)),
                       Change2 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)),
                       Change3 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))
})

test_that("Population after 3 ordered interventions Exp",{
  popSize1 <- 1000
  fractionPopSize1 <- 0.15
  popSize2 <- 3000
  fractionPopSize2 <- 0.1
  popSize3 <- 5000
  fractionPopSize3 <- 0.01
  
  
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize1),
                                       action = list(fractionPopSize = fractionPopSize1)),
                       Change2 = list( trigger = list(popSize = popSize2),
                                       action = list(fractionPopSize = fractionPopSize2)),
                       Change3 = list( trigger = list(popSize = popSize3),
                                       action = list(fractionPopSize = fractionPopSize3)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))
})

test_that("Population after 3 disordered interventions Exp",{
  popSize1 <- 3000
  fractionPopSize1 <- 0.15
  popSize2 <- 5000
  fractionPopSize2 <- 0.1
  popSize3 <- 1000
  fractionPopSize3 <- 0.15
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize1),
                                       action = list(fractionPopSize = fractionPopSize1)),
                       Change2 = list( trigger = list(popSize = popSize2),
                                       action = list(fractionPopSize = fractionPopSize2)),
                       Change3 = list( trigger = list(popSize = popSize3),
                                       action = list(fractionPopSize = fractionPopSize3)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))
})

test_that("Population after 3 repeated interventions McFL",{
  popSize <- 6000
  fractionPopSize <- 0.10
  sampleEvery <- 0.01
  
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)),
                       Change2 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)),
                       Change3 = list( trigger = list(popSize = popSize),
                                       action = list(fractionPopSize = fractionPopSize)))
  
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))
})

test_that("Population after 3 ordered interventions McFL",{
  popSize1 <- 1000
  fractionPopSize1 <- 0.15
  popSize2 <- 3000
  fractionPopSize2 <- 0.1
  popSize3 <- 5000
  fractionPopSize3 <- 0.01
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize1),
                                       action = list(fractionPopSize = fractionPopSize1)),
                       Change2 = list( trigger = list(popSize = popSize2),
                                       action = list(fractionPopSize = fractionPopSize2)),
                       Change3 = list( trigger = list(popSize = popSize3),
                                       action = list(fractionPopSize = fractionPopSize3)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))

})

test_that("Population after 3 disordered interventions McFL",{
  popSize1 <- 3000
  fractionPopSize1 <- 0.15
  popSize2 <- 5000
  fractionPopSize2 <- 0.1
  popSize3 <- 1000
  fractionPopSize3 <- 0.15
  modelChanges = list( Change1 = list( trigger = list(popSize = popSize1),
                                       action = list(fractionPopSize = fractionPopSize1)),
                       Change2 = list( trigger = list(popSize = popSize2),
                                       action = list(fractionPopSize = fractionPopSize2)),
                       Change3 = list( trigger = list(popSize = popSize3),
                                       action = list(fractionPopSize = fractionPopSize3)))
  nd <- 70  
  np <- 0 
  s <- 0.1  
  sp <- 1e-3 
  spp <- -sp/(1 + sp)
  mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                            drvNames = seq.int(nd))
  s1 <-  oncoSimulIndiv(mcf1,
                        model = "Exp", 
                        mu = 1e-5,
                        detectionProb = NA,
                        detectionSize = 1e4, 
                        detectionDrivers = NA,
                        sampleEvery = sampleEvery,
                        keepEvery = 0.01,
                        initSize = 1000,
                        finalTime = 2000,
                        onlyCancer = FALSE, 
                        modelChanges = modelChanges)
  plot(s1, addtot = TRUE, lwdClone = 0.9, log = "", 
       thinData = TRUE, thinData.keep = 0.5) # delete
  
  p_time <- s1$pops.by.time
  # Take the pop size before intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time))
  pop_before_interv <- rowSums(df[,2:ncol(df)])
  
  # Take the pop size after intervention
  df <- as.data.frame(p_time) %>% 
    filter(V1 %in% (s1$other$interv_time + sampleEvery))
  pop_after_interv <- rowSums(df[,2:ncol(df)])
  
  # Obtain the expected value of pop
  fraction_pop <- pop_before_interv * fractionPopSize
  
  expect_true(all(pop_after_interv <= fraction_pop))
})

