myTestBarLevels = function(){
  mybys = c(-1,-2,-3,-4)
  res = NULL
  index = 1
  for(myby in mybys){
    barLevels = paste0(unlist(lapply(paste0("E-",seq(8,1,myby)),function(x){ paste0(seq(9,1,-2),x) })),collapse=",")
    partial = selectLearnEval(barLevels=barLevels,forcePRSice = T,outcome="disc")
    res[[index]] = partial
    print(barLevels)
    cat("Number of SNPs",partial$prsicesummary$Num_SNP,"\n")
    print(partial$confMat$byClass)
    print(partial$evalConfMat)
    index = index + 1
  }
  res
}

myTestIncs = function(){
  slower = 0
  #sincs = c(0.05,0.025,0.016,0.0125,0.01)
  sincs = c(0.05,0.025,0.016)
  supper = 0.5

  res = NULL
  index = 1
  for(myinc in sincs){
    partial = selectLearnEval(barLevels=NULL,
                              slower=slower,
                              sinc=myinc,
                              supper=supper,
                              forcePRSice = T,
                              outcome="disc",
                              gwas="small.tab")
    res[[index]] = partial
    cat("Number of SNPs",partial$prsicesummary$Num_SNP,"\n")
    print(partial$confMat$byClass)
    print(partial$evalConfMat)
    index = index + 1
  }
  res
}

myTestClumpingKBs = function(){
  options(scipen=999)
  kbs = c(850,500,250,5)

  res = NULL
  index = 1
  for(kb in kbs){
    partial = selectLearnEval(clumpkb=kb,
                              prsiceseed = 123,
                              forcePRSice = T,
                              barLevels=paste0(unlist(lapply(paste0("E-",seq(8,1,-1)),function(x){ paste0(seq(4,1,-2),x) })),collapse=","),
                              outcome="disc",
                              gwas="small.tab")
    res[[index]] = partial
    cat("Number of SNPs",partial$prsicesummary$Num_SNP,"\n")
    print(partial$confMat$byClass)
    print(partial$evalConfMat)
    index = index + 1
  }
  res
}

myTestClumpingR2 = function(r2s = c(0.01,0.05,0.1,0.2)){
  options(scipen=999)


  res = NULL
  index = 1
  for(r2 in r2s){
    partial = selectLearnEval(clumpr2=r2,
                              prsiceseed = 123,
                              forcePRSice = T,
                              barLevels=paste0(unlist(lapply(paste0("E-",seq(8,1,-1)),function(x){ paste0(seq(4,1,-2),x) })),collapse=","),
                              outcome="disc",
                              gwas="small.tab")
    res[[index]] = partial
    cat("Number of SNPs",partial$prsicesummary$Num_SNP,"\n")
    print(partial$confMat$byClass)
    print(partial$evalConfMat)
    index = index + 1
  }
  res
}
