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
       thinData = TRUE, thinData.keep = 0.5)
  
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
  
  salida <- testthat:::expect_true(pop_after_interv <= fraction_pop) #all.equal(pop_after_interv, fraction_pop, tol=0)
  
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
       thinData = TRUE, thinData.keep = 0.5)
  
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
  
  salida <- testthat:::expect_true(pop_after_interv <= fraction_pop)
  
})
