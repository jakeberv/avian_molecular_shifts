######################################################################################################1
#Bayou allometry shifting regimes model pipeline 
#written by: Sean McHugh, Josef Uyeda, Jacob Berv

#estimates evoltuion of an evolutionary allometric slope and intercept as multiregime OU
#process using the R package bayou (Uyeda and Harmon 2014; Uyeda et al. 2017), with shifts 
#either provided a priori or estimated via reversible jump MCMC

# includes marginal likelihood estimation using the stepping stone sampling method (Xie et al. 2011)
# to compare model support for fixed, reversible jump, and global regime null models   
#-- but this currently does not work, so, skipped for now
########################################################################################################1


#load packages
library(ape)
library(coda)
library(geiger)
library(ips)
library(phytools)
library(lmtest)
library(treeplyr)
library(bayou)
library(doParallel)
require(viridis)
require(pbmcapply)

#this next section will run the analyses -- 
#but it is separated with brackets { } for organization
#skip to line 459 to load pre-analyzed RDS file and proceed with
#violin plot generation, as in the main text

{
#install.packages("devtools")
#require(devtools)
#install_github("uyedaj/bayou")
#require(bayou)

#################
# Prepare Trees #
################
{
  
  
  # put your directory here
  setwd("SET DIRECTORY TO TOP LEVEL")
  setwd("./BMR")
  
  #make sure you have these destinations directories set to your machines path and 
  #ending as here (ex fixed and saved_obejcts must be the terminal directories)
  dir.create("./fixed")
  dir.create("./saved_objects")
  
  #import the tree RDS file
  simmap.janus.nuc.alldata<-readRDS(file="simmap.janus.nuc.alldata.RDS")
 
  #get regime mapping across tips
  plot(simmap.janus.nuc.alldata)
  edgelabels(cex=.4, frame="c")
  
  #simmap.janus.nuc.alldata$mapped.edge
  #getStates(simmap.janus.nuc.alldata, type='tips')
  
  #data (WITH missing data) in log space
  #saveRDS(scaling_data, file="scaling_data.RDS")
  #scaling_data<-readRDS(file="scaling_data.RDS")
  
  #data with imputed values for missing data in log space
  #data imputed under mvBM, multi-regime
  #nomiss<-mvbm.imputed$estimates[,c(which(colnames(mvbm.imputed$estimates)==Y),which(colnames(mvbm.imputed$estimates)==X)),drop=F]
  #saveRDS(nomiss, file="nomiss.RDS")
  nomiss<-readRDS(file="nomiss.RDS")
  
  #squared standard errors in log space
  #saveRDS(scaling_data_errors, file="standarderror_sq.RDS")
  #standarderror_sq<-readRDS(file="standarderror_sq.RDS")
  
  ##t<-sample(tortoise_chronogram,size=100)
  #
  #nodelabels()
  
  #transform simmap to normal phylo object
  tree <-  ape::reorder.phylo(as.phylo(simmap.janus.nuc.alldata), order="postorder")
  

  #add standard error to dataframe with traits, kinda not necessary in current implementation but thought it may be useful down the road
  birb_traits<- cbind(nomiss)#,  "se_bmr"=standarderror_sq[,2], "se_mass" = standarderror_sq[,1])
  

  #SSDI_no_error        <-  TortTraits$dat$SSDI
  bmr<-    birb_traits[,2]
  mass<-  birb_traits[,1]
  

  #se_bmr<-birb_traits[,3]
  
  #se_mass<-birb_traits[,4]
  
  
  MEvar <- 0.01
  #MEvar <- 0.1
  
}

{

#fixed_branches<-identifyBranches(tree, 7, fixed.loc = TRUE, plot.simmap = TRUE)
#I identified branches with identifyBranches function, this is manual branch selection so goodluckhavefun
fixed_branches<-list(uncex.merged.mtdnas=c(13, 14,  15, 199, 265, 370, 374, 387, 388, 389, 390, 393),
                     uncex.merged       =c(13, 14,  15, 199, 265, 370, 374, 387, 388, 389, 390),
                     uncex.fixed_rjMCMC   =c(13, 14, 199, 265, 374, 387, 388, 390),
                     uncex.exons        =c(13,  14 ,199 ,265 ,370 ,374, 389),
                     uncex.introns      =c(13, 265, 387, 388, 390),
                     uncex.utrs         =c(15, 265, 387),
                     uncex.mtdnas.all   =c(14, 393),
                     mtDNAs.proteins    =c(393)
)
                     
  
#t2 is the regimes I identified 2:7 regimes (1 is the starting regime), one unique regime for each shift, we could also set this up to have convergent regimes as shown in the molecular shifts
#
t2 <-pbmclapply(1:length(fixed_branches), function(shifts) c(2:(length(fixed_branches[[shifts]])+1)), mc.cores = 50)

# parameter for where on the branch eahc shift occurs locations shifts a location of 0 approximates the node the branch follows from
loc<-pbmclapply(1:length(fixed_branches), function(shifts) sample(0:0, length(fixed_branches[[shifts]]), replace=T), mc.cores = 50)

#sample(0:0, length(fixed_branches), replace=T)
  #names(SSDI[[rep]])<-TortTraits$phy$tip.label
  names(bmr )<-tree$tip.label
  names(mass)<-tree$tip.label
  
  


plot.phylo(tree, type= "phylogram", use.edge.length = T, show.node.label = T) #plot incase extracting clades interactively, plus just to check it ou
nodelabels()


#pred_list<-list(mass,bmr, mass)

# mass is the allometry predictor
pred_list<-mass

#Dat_list<-list(SSDI, SSDI, bmr)

#bmr is the data
Dat_list<-bmr




mcmc.list<-list()

#i=1
reps<-3

}


#setup 4 reps to setup the directories and models for running the 5 allometric models:


# F = fixed model global slope and fixed shifts on intercept
# FNN = fixed model shifts on slope and intercept
# 11 = global intercept and slope
# N1 = Global slope and shifting intercept
# NN = Shifting slope and Intercept

#each run the number of times specified by reps (I usually do 4 but up to you)

for(i in 1:reps){
  
  
  # this is usually where data is scaled hence the object names, but since it came in transformed we are good and just transform trait object into a matrix
  
  Pred_unscale<-as.matrix(as.numeric(pred_list), row.names = names(pred_list))
  

  Pred<-Pred_unscale
  
  names(Pred)<-names(mass)
  
  

  Dat<-Dat_list
  
  names(Dat)<-names(bmr)
  

  
  length(Dat)
  length(Pred)
  
#shifts=4
  
prior.F=list()
model.F=list()
model.FNN=list()
mcmc.F=list()
mcmc.FNN=list()
  ######setup fixed models first######
  
  for (shifts in 1:length(fixed_branches)){
  
    # fixed model prior (same for FNN and F)
    #comments on this prior should apply to the rest of the priors blocks
    
    prior.F[[shifts]] <- make.prior(tree, plot.prior = F,  #you can plot distributions but it gets busy very quickly
                          #setup distributions
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_Pred="dnorm",
                                     dsb="fixed", dk="fixed", dtheta="dnorm"), 
                          #setup pars for distributions
                          param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                     dbeta_Pred=list(mean=0.7, sd=0.1),
                                     dtheta=list(mean=-3.5, sd=1.75)),
                          #what have you fixed in the model, this is the fixed model so provide your branch numbers, the number of shifts, and the number of thetas in the model
                          #other models without fixed shifts assumed will have zeros here as seen below
                          fixed=list(k=length(fixed_branches[[shifts]]), sb=fixed_branches[[shifts]], ntheta=length(fixed_branches[[shifts]])+1, t2=t2[[shifts]])
    )
    
    names(fixed_branches)[[shifts]]
    
    
    DF=   list(alpha=2, sig2=2, beta_Pred=0.1, k=1, theta=0.5, slide=1)
    #DNN = list(alpha=2, sig2=2, beta_Pred=0.3, k=c(1,1), theta=2, slide=1)
    DNN = list(alpha=2, sig2=2, beta_Pred=0.3, k=c(1,1), theta=2, slide=1)
    

    file_dir<-paste("/fixed/",names(fixed_branches)[[shifts]],i, sep="")
    
    
    
    model.F[[shifts]]<- makeBayouModel(Dat ~ Pred, rjpars = c(), #rjpars indicate what pars we assume will have shifts w/o a priori assumptions, in the fixed model this is none but in others it may be theta or beta
                             tree=tree, dat=Dat, pred=Pred, SE=MEvar, prior=prior.F[[shifts]], D=DF,
                             startpar = list(alpha=1, sig2=1, k=length(fixed_branches[[shifts]]), ntheta=length(fixed_branches[[shifts]])+1, beta_Pred= 0,theta=sample(0:0, length(fixed_branches[[shifts]])+1, replace=T), sb=fixed_branches[[shifts]], t2=t2[[shifts]], loc=loc[[shifts]]),
                             slopechange = "alphaWeighted"
    )
    
    
    
    
    
    model.FNN[[shifts]]<- makeBayouModel(Dat ~ Pred, rjpars = c("theta", "Pred"), 
                               tree=tree, dat=Dat, pred=Pred, SE=MEvar, prior=prior.F[[shifts]], D=DNN,
                               startpar = list(alpha=1, sig2=1, k=length(fixed_branches[[shifts]]), ntheta=length(fixed_branches[[shifts]])+1, beta_Pred= sample(0:0, length(fixed_branches[[shifts]])+1, replace=T),theta=sample(0:0, length(fixed_branches[[shifts]])+1, replace=T), sb=fixed_branches[[shifts]], t2=t2[[shifts]], loc=loc[[shifts]]), 
                               slopechange = "alphaWeighted"
    )
    
    
    mcmc.F[[shifts]]   <- bayou.makeMCMC(tree, Dat, pred=Pred, SE=MEvar, model=model.F[[shifts]]$model, prior=prior.F[[shifts]], startpar=model.F[[shifts]]$startpar,    file.dir = file_dir,  outname="modelF_r001", plot.freq=NULL,samp=samp)
    
    mcmc.FNN[[shifts]] <- bayou.makeMCMC(tree, Dat, pred=Pred, SE=MEvar, model=model.FNN[[shifts]]$model, prior=prior.F[[shifts]], startpar=model.FNN[[shifts]]$startpar, file.dir = file_dir, outname="modelFNN_r001", plot.freq=NULL,samp=samp)
    
    
    
  
  }
  
names(mcmc.F)<-names(fixed_branches)
names(mcmc.FNN)<-names(fixed_branches)




  #######edited_allo#############
  ###allommeetric model###
  #set up priors
  
  
  
  #prior list 
  
  
  # alpha: OU adaptive pull
  # sig2:  the OU process step variance (how much change occurs at each step in the OU process) 
  #        note: this does not equate to total trait displacement like BM as it is constrained towards a optima over all steps 
  # beta_Pred: for allometric slope
  # sb:  shift locations (branch the shift is on)
  # k:  number of shifts 
  # theta: dallometric intercept
  
  
  
  
  
  #fixed=list(k=(number of shifts), sb=(branch location of the shift)
  #everything global

  
  
#models labelled the same across prior, model, and mcmc objects
# bmax = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
#          0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0,
#          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
#          1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
#          1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
#          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
#          1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0,
#          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

bmax = rep(0,394)
#next restrict edges to candidate edges only
bmax[c(13,  14,  15, 199, 265, 370, 374, 387, 388, 389, 390, 393)]=1

edges<- tree$edge.length

edges[bmax==0]=0
edges[bmax==1]=1 #uncomment to make probability of shifts equal across branches


#Global null (Global intercept AND slope)
prior.11 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_Pred="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_Pred=list(mean=0.7, sd=0.1),
                                  dtheta=list(mean=-3.5, sd=1.75)), #list(mean=-3.5, sd=1.75)
                       fixed=list(k=0, sb=numeric(0))
)
#Global null (Global slope ONLY) standard model only intercept changes (theta) 
prior.N1 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_Pred="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_Pred=list(mean=0.7, sd=0.1),
                                  dsb=list(bmax=Inf,prob=edges),
                                  dk=list(lambda=8, kmax=24),
                                  dtheta=list(mean=-3.5, sd=1.75))
)

# full reversible jump slope and intercept changing (everything changes, slope and intercept)
prior.NN <- make.prior(tree, plot.prior = T, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_Pred="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_Pred=list(mean=0.7, sd=0.1),
                                  dsb=list(bmax=Inf,prob=edges),
                                  dk=list(lambda=8, kmax=24),
                                  dtheta=list(mean=-3.5, sd=1.75))
)


  
  #set tuning pars for MCMC######
  #fiddle with these as needed to get good convergence (see gelmans R)
  D11 = list(alpha=2, sig2=2, beta_Pred=0.1, k=1, theta=0.5, slide=1)
  DN1 = list(alpha=2, sig2=2, beta_Pred=0.1, k=1, theta=2, slide=1)
  #DNN = list(alpha=2, sig2=2, beta_Pred=0.3, k=c(1,1), theta=2, slide=1)
  
  
  
  #setup model###########
  
  #setup model by plugging in the priors you made before, your data, and the starting points for each parameter

  
  model.11 <- makeBayouModel(Dat ~ Pred, rjpars = c(), 
                             tree=tree, dat=Dat, pred=Pred, SE=MEvar, prior=prior.11, D=D11, slopechange = "alphaWeighted")
  
  model.N1 <- makeBayouModel(Dat ~ Pred, rjpars = c("theta"),  
                             tree=tree, dat=Dat, pred=Pred, SE=MEvar, prior=prior.N1, D=DN1, slopechange = "alphaWeighted")
  
  model.NN <- makeBayouModel(Dat ~ Pred, rjpars = c("theta", "Pred"),  
                             tree=tree, dat=Dat, pred=Pred, SE=MEvar, prior=prior.NN, D=DNN, slopechange = "alphaWeighted")
  
  #bayou.makeMCMC
  
  model.F$startpar
  
  
  # change names to whatever directory you want
  
  file_dir<-paste("","global",i, sep="")
  
  
  # chose MCMC chain length here and how many iterations you keep in your posterior
  iterations=10000000
  samp=1000#iterations/5000

  
  #make the mcmc objec tthat you will execute for each model
  
  mcmc.11  <- bayou.makeMCMC(tree, Dat, pred=Pred, SE=MEvar, model=model.11$model, prior=prior.11, startpar=model.11$startpar, file.dir = file_dir,  outname="model11_r001", plot.freq=NULL,samp=samp)
  
  mcmc.N1  <- bayou.makeMCMC(tree, Dat, pred=Pred, SE=MEvar, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, file.dir = file_dir,  outname="modelN1_r001", plot.freq=NULL, samp=samp)
  
  mcmc.NN  <- bayou.makeMCMC(tree, Dat, pred=Pred, SE=MEvar, model=model.NN$model, prior=prior.NN, startpar=model.NN$startpar, file.dir = file_dir,  outname="modelNN_r001", plot.freq=NULL, samp=samp)
  
  
  # create the list to pass into doParallel
  
  
  
  mcmc.list[[i]]<-c(mcmc.F,  
                    mcmc.FNN,
                    list(mcmc.11,
                         mcmc.N1, 
                         mcmc.NN ))
  
  name_vec<-c(paste( "F_", names(fixed_branches), sep=""),  paste( "FNN_", names(fixed_branches), sep=""), "mcmc.11", "mcmc.N1", "mcmc.NN")
  
  
  names(mcmc.list[[i]]) = paste(i, "_", name_vec, sep="")
  
      
  
}

}

#names(mcmc.list[[i]])
#saveRDS(object = mcmc.list, file = "saved_objects/mcmcs_full_batch.rds")
#mcmc.list<-readRDS(file="saved_objects/mcmcs_full_batch.rds")

#new location
mcmc.list <- readRDS(file = "./saved_objects/mcmcs_full_batch.rds")
mcmc.full.list<- unlist(mcmc.list, recursive=F)

#View(unlist(mcmc.list, recursive=F))

{

#View(mcmc.full.list)

{
  registerDoParallel(cores=length(mcmc.full.list))
  foreach(i=1:length(mcmc.full.list)) %dopar% {
    mcmc.full.list[[i]]$run(iterations)
  }
}


  
}  


########check run#########



chain.list<- pbmclapply(mcmc.full.list, function(mcmc) mcmc$load(), mc.cores=50)
saveRDS(object = chain.list, file = "saved_objects/chains_full_batch.rds")

#chain.list<-readRDS(file="saved_objects/chains_full_batch.rds")
chain.list <- readRDS(file = "./saved_objects/chains_full_batch.rds")

#names(mcmc.full.list)
names(chain.list) <-names(mcmc.full.list)


#group chains into objects by model
F_uncex.merged.mtdnas.chains  <-chain.list[grepl("F_uncex.merged.mtdnas" ,  names(chain.list)) ]
F_uncex.merged.chains         <-chain.list[grepl("F_uncex.merged"       ,  names(chain.list)) ]
F_uncex.fixed_rjMCMC.chains     <-chain.list[grepl("F_uncex.fixed_rjMCMC"   ,  names(chain.list)) ]
F_uncex.exons.chains          <-chain.list[grepl("F_uncex.exons"        ,  names(chain.list)) ]
F_uncex.introns.chains        <-chain.list[grepl("F_uncex.introns"      ,  names(chain.list)) ]
F_uncex.utrs.chains           <-chain.list[grepl("F_uncex.utrs"         ,  names(chain.list)) ]
F_uncex.mtdnas.all.chains     <-chain.list[grepl("F_uncex.mtdnas.all"   ,  names(chain.list)) ]
F_mtDNAs.proteins.chains      <-chain.list[grepl("F_mtDNAs.proteins"    ,  names(chain.list)) ]

FNN_uncex.merged.mtdnas.chains<-chain.list[grepl("FNN_uncex.merged.mtdnas" ,  names(chain.list)) ]
FNN_uncex.merged.chains       <-chain.list[grepl("FNN_uncex.merged"       ,  names(chain.list)) ]
FNN_uncex.fixed_rjMCMC.chains   <-chain.list[grepl("FNN_uncex.fixed_rjMCMC"   ,  names(chain.list)) ]
FNN_uncex.exons.chains        <-chain.list[grepl("FNN_uncex.exons"        ,  names(chain.list)) ]
FNN_uncex.introns.chains      <-chain.list[grepl("FNN_uncex.introns"      ,  names(chain.list)) ]
FNN_uncex.utrs.chains         <-chain.list[grepl("FNN_uncex.utrs"         ,  names(chain.list)) ]
FNN_uncex.mtdnas.all.chains   <-chain.list[grepl("FNN_uncex.mtdnas.all"   ,  names(chain.list)) ]
FNN_mtDNAs.proteins.chains    <-chain.list[grepl("FNN_mtDNAs.proteins"    ,  names(chain.list)) ]


NN.chains                     <-chain.list[grepl("c.NN",  names(chain.list)) ]
N1.chains                     <-chain.list[grepl("c.N1",  names(chain.list)) ]
global.chains                 <-chain.list[grepl("c.11",  names(chain.list)) ]

#as.mcmc.list(lapply(F.chains, function(x) as.mcmc(x$sig2) ) )

#create new list of chains from the grouped by model chain objects
chains<-list(F_uncex.merged.mtdnas  =F_uncex.merged.mtdnas.chains,
             F_uncex.merged         =F_uncex.merged.chains,       
             F_uncex.fixed_rjMCMC     =F_uncex.fixed_rjMCMC.chains,   
             F_uncex.exons          =F_uncex.exons.chains,
             F_uncex.introns        =F_uncex.introns.chains,
             F_uncex.utrs           =F_uncex.utrs.chains,
             F_uncex.mtdnas.all     =F_uncex.mtdnas.all.chains,
             F_mtDNAs.proteins      =F_mtDNAs.proteins.chains,
             FNN_uncex.merged.mtdnas=FNN_uncex.merged.mtdnas.chains,
             FNN_uncex.merged       =FNN_uncex.merged.chains,       
             FNN_uncex.fixed_rjMCMC   =FNN_uncex.fixed_rjMCMC.chains,   
             FNN_uncex.exons        =FNN_uncex.exons.chains,
             FNN_uncex.introns      =FNN_uncex.introns.chains,
             FNN_uncex.utrs         =FNN_uncex.utrs.chains,
             FNN_uncex.mtdnas.all   =FNN_uncex.mtdnas.all.chains,
             FNN_mtDNAs.proteins    =FNN_mtDNAs.proteins.chains,
             NN                     =NN.chains, 
             N1                     =N1.chains, 
             global                 =global.chains                 )




#####testing convergence of chains using Gelman's R########

# using this a a measure to quickly asses convergance either to choose a burnin point or see if the run was successful, 
#the closer to 1 statistic over the chain the the better, ideally it is >1.005 but that can be a bit hopeful, even if 
#you dont have convergance combining chains can still be recommended

#Note: RJ pars cant be tested since it is difficult to isolate separate chains for each possible shift. Maybe its possible 
#but I dont know how and again this isnt a "word of god" test as much as a diagnostic as implemented here


# pars helper object for labeling figures
pars=c("sig2",
       "alpha",
       "beta_Pred",
       "ntheta",
       "k")

for (model in 1:length(chains)){
  
  pdf(file= paste(names(chains)[[model]],"_gelmans_plots.pdf", sep=""))    
  
  
  
  
  gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$sig2) ) ), autoburnin = T);
  title("sig2")
  gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$alpha) ) ), autoburnin = T);
  title("alpha")
  
  if(length(chains[[model]][[1]]$beta_Pred[[1]])>1){
    gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc( unlist( lapply(x$beta_Pred, function(x) x[[1]]) ) ) ) ), autoburnin = T);
  } else{
    gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$beta_Pred) ) ), autoburnin = T);
    
  }
  title("beta root")
  
  if(length(chains[[model]][[1]]$theta[[1]])>1){
    
    gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc( unlist( lapply(x$theta, function(x) x[[1]]) ) ) ) ), autoburnin = T);
    
  }else{
    gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$ntheta) ) ), autoburnin = T);
  }
  title("theta root")
  
  try(gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$ntheta) ) ), autoburnin = T));
  title("ntheta")
  
  try(gelman.plot(as.mcmc.list(lapply(chains[[model]], function(x) as.mcmc(x$k) ) ), autoburnin = T));
  title("k")
  
  dev.off()  
  
}

F_uncex.merged.mtdnas.combined.chains  <-bayou::combine.chains(F_uncex.merged.mtdnas.chains  , burnin.prop=0.4)
F_uncex.merged.combined.chains         <-bayou::combine.chains(F_uncex.merged.chains         , burnin.prop=0.4)
F_uncex.fixed_rjMCMC.combined.chains     <-bayou::combine.chains(F_uncex.fixed_rjMCMC.chains     , burnin.prop=0.4)
F_uncex.exons.combined.chains          <-bayou::combine.chains(F_uncex.exons.chains          , burnin.prop=0.4)
F_uncex.introns.combined.chains        <-bayou::combine.chains(F_uncex.introns.chains        , burnin.prop=0.4)
F_uncex.utrs.combined.chains           <-bayou::combine.chains(F_uncex.utrs.chains           , burnin.prop=0.4)
F_uncex.mtdnas.all.combined.chains     <-bayou::combine.chains(F_uncex.mtdnas.all.chains     , burnin.prop=0.4)
F_mtDNAs.proteins.combined.chains      <-bayou::combine.chains(F_mtDNAs.proteins.chains      , burnin.prop=0.4)
FNN_uncex.merged.mtdnas.combined.chains<-bayou::combine.chains(FNN_uncex.merged.mtdnas.chains, burnin.prop=0.4)
FNN_uncex.merged.combined.chains       <-bayou::combine.chains(FNN_uncex.merged.chains       , burnin.prop=0.4)
FNN_uncex.fixed_rjMCMC.combined.chains   <-bayou::combine.chains(FNN_uncex.fixed_rjMCMC.chains   , burnin.prop=0.4)
FNN_uncex.exons.combined.chains        <-bayou::combine.chains(FNN_uncex.exons.chains        , burnin.prop=0.4)
FNN_uncex.introns.combined.chains      <-bayou::combine.chains(FNN_uncex.introns.chains      , burnin.prop=0.4)
FNN_uncex.utrs.combined.chains         <-bayou::combine.chains(FNN_uncex.utrs.chains         , burnin.prop=0.4)
FNN_uncex.mtdnas.all.combined.chains   <-bayou::combine.chains(FNN_uncex.mtdnas.all.chains   , burnin.prop=0.4)
FNN_mtDNAs.proteins.combined.chains    <-bayou::combine.chains(FNN_mtDNAs.proteins.chains    , burnin.prop=0.4)
NN.combined.chains                     <-bayou::combine.chains(NN.chains                     , burnin.prop=0.4)
N1.combined.chains                     <-bayou::combine.chains(N1.chains                     , burnin.prop=0.4)
global.combined.chains                 <-bayou::combine.chains(global.chains                 , burnin.prop=0.4)

#create list of combined chains
combined.chainlist<-list(F_uncex.merged.mtdnas  =F_uncex.merged.mtdnas.combined.chains,
             F_uncex.merged         =F_uncex.merged.combined.chains,       
             F_uncex.fixed_rjMCMC     =F_uncex.fixed_rjMCMC.combined.chains,   
             F_uncex.exons          =F_uncex.exons.combined.chains,
             F_uncex.introns        =F_uncex.introns.combined.chains,
             F_uncex.utrs           =F_uncex.utrs.combined.chains,
             F_uncex.mtdnas.all     =F_uncex.mtdnas.all.combined.chains,
             F_mtDNAs.proteins      =F_mtDNAs.proteins.combined.chains,
             FNN_uncex.merged.mtdnas=FNN_uncex.merged.mtdnas.combined.chains,
             FNN_uncex.merged       =FNN_uncex.merged.combined.chains,       
             FNN_uncex.fixed_rjMCMC   =FNN_uncex.fixed_rjMCMC.combined.chains,   
             FNN_uncex.exons        =FNN_uncex.exons.combined.chains,
             FNN_uncex.introns      =FNN_uncex.introns.combined.chains,
             FNN_uncex.utrs         =FNN_uncex.utrs.combined.chains,
             FNN_uncex.mtdnas.all   =FNN_uncex.mtdnas.all.combined.chains,
             FNN_mtDNAs.proteins    =FNN_mtDNAs.proteins.combined.chains,
             NN                     =NN.combined.chains, 
             N1                     =N1.combined.chains, 
             global                 =global.combined.chains                 )



# set burning over all chains
burnin.combined.chain.list<- pbmclapply(combined.chainlist, function(chain) set.burnin(chain, 0.4), mc.cores=10)

#summarize all MCMC chains
summary.list<- pbmclapply(burnin.combined.chain.list, function(chain) summary(chain), mc.cores=10)

#shift summaries (used to make the 4 panel plot)
shiftSummaries.list<-pbmclapply(1:length(burnin.combined.chain.list), function(model)  try(shiftSummaries(burnin.combined.chain.list[[model]], mcmc =mcmc.list[[1]][[model]] , pp.cutoff=.001) ), mc.cores=10)

names(shiftSummaries.list)<-names(burnin.combined.chain.list)


#making 4 panel plot of shifting allometries
{
  pdf(file= paste("shift_summaries_combined_full-separate.pdf", sep=""))    
  par(mfrow = c(2, 2))
  
  lapply(1:length(shiftSummaries.list), function(SS){ if(length(shiftSummaries.list[[SS]])>1){ plotShiftSummaries(shiftSummaries.list[[SS]], lwd=2, single.plot=T, label.pts=F );
    mtext(names(shiftSummaries.list)[[SS]], side = 3, line = -3, outer = TRUE)}
    
  })
  
  dev.off()
}


cols<-function(n){
  set.seed(3)
  head(c("black", sample(rcartocolor:::carto_pal(12, "Safe"))), n)
}


#' A function to plot a list produced by \code{shiftSummaries}
plot_shift_sum <- function(summaries, pal=rainbow, ask=FALSE, single.plot=FALSE, label.pts=TRUE, ...){
  #oldpar <- graphics::par(no.readonly = TRUE)    # code line i
  #on.exit(graphics::par(oldpar))            # code line i + 1 
  #px <- par()
  ndens <- length(summaries$cladesummaries[[1]]$densities)
  #par(mfrow=c(2,max(ndens,2)), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white", ...)
  blank.panels <- prod(par()$mfrow) - (2+ndens)
  #par(ask=ask)
  regressions <- summaries$regressions
  if(ncol(regressions)==1){ regressions <- data.frame(regressions, "slope"=0)}
  dat <- summaries$dat
  tree <- summaries$tree
  sumpars <- summaries$pars
  descendents <- summaries$descendents
  PP <- c("Root",round(summaries$PP,2))
  xlimits <- apply(do.call(rbind, 
                           lapply(summaries$cladesummaries, function(x) 
                             sapply(x$densities, function(y) range(y$x))
                           )), 2, range)
  xlimits[1,] <- xlimits[1,]-0.1*apply(xlimits, 2, diff)
  xlimits[2,] <- xlimits[2,]+0.1*apply(xlimits, 2, diff)
  #print(xlimits)
  if(ndens > 1){
    xint <- setNames(data.frame(summaries$pred)[[1]], names(dat))
    xlimits2 <- range(xint)
    #print(xlimits2)
  } else {
    xint <- jitter(.tipregime(sumpars, tree))
    xlimits2 <- c(-2, sumpars$ntheta+3)
  }
  if(!single.plot){
    for(i in (1:nrow(regressions))){
      plotBayoupars(sumpars, tree, col=setNames(c(pal(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), cex=0.2)
      plot(xint, dat, pch=21, xlim=xlimits2, bg=makeTransparent("gray80", 100), col =makeTransparent("gray80", 10), main=paste("Posterior prob: ", PP[i], sep=""))
      if(length(descendents[[i]] > 0)){
        if(label.pts) text(xint[descendents[[i]]], dat[descendents[[i]]], labels=names(dat[descendents[[i]]]), col="white", cex=0.4, pos = 2)
        points(xint[descendents[[i]]], dat[descendents[[i]]], pch=21, bg=makeTransparent(pal(nrow(regressions))[i], 100), col =makeTransparent(pal(nrow(regressions))[i], 10))
      } else{
        warnings("No descendents for this shift")
      }
      abline(a=regressions[i,1], b=regressions[i,2], col=pal(nrow(regressions))[i], lwd=2, lty=2)
      dens <- summaries$cladesummaries[[i]]$densities
      gbg <- lapply(1:length(dens), function(y)plot(dens[[y]], col=pal(nrow(regressions))[i], main=names(dens)[y],xlim=c(xlimits[,y])))
      
      if(blank.panels >0){lapply(1:blank.panels,function(x) plot.new())}
    }
  } else {
    plotBayoupars(sumpars, tree, col=setNames(pal(sumpars$ntheta), 1:sumpars$ntheta), cex=0.2, tip.col="white", no.margin=T)
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    plot(xint, dat, pch=21, xlim=xlimits2, bg=makeTransparent("gray80", 100), col = makeTransparent("gray80", 10), bty='n')
    for(i in 1:length(descendents)){
      abline(a=regressions[i,1], b=regressions[i,2], col=pal(nrow(regressions))[i], lwd=0.5, lty=1)
      if(length(descendents[[i]] > 0)){
        if(label.pts) text(xint[descendents[[i]]], dat[descendents[[i]]], labels=names(dat[descendents[[i]]]), col="white", cex=0.4, pos = 2)
        par(lwd=0.05)
        points(xint[descendents[[i]]], dat[descendents[[i]]], pch=21, bg=makeTransparent(pal(nrow(regressions))[i], 100), 
               col = "black", cex=1)
      } else{
        warnings("No descendents for this shift")
      }
      
    }
    varN <- length(summaries$cladesummaries[[1]]$densities)
    for(j in 1:varN){
      dens <- lapply(1:length(summaries$cladesummaries), function(x) summaries$cladesummaries[[x]]$densities[[j]])
      xrange <- quantile((do.call(c, lapply(dens, function(q) q$x))), c(0.3, 0.8))
      #print(xrange)
      ymax <- max(do.call(c, lapply(dens, function(q) q$y)))
      plot(xrange, c(0, ymax*1.1), type="n", main=names(summaries$cladesummaries[[1]]$densities)[j], xlab="", ylab="Kernel Density Estimate", bty='n')
      gbg <- lapply(1:length(dens), function(y) {polygon(dens[[y]], col=adjustcolor(pal(nrow(regressions))[y], alpha.f=1), 
                                                        lwd=0.0001, border=pal(nrow(regressions))[y]); 
        abline(lty=2, lwd=1, lend=2, v=dens[[y]]$x[which.max(dens[[y]]$y)], col=pal(nrow(regressions))[y]) })
    }
    
  }
  #px <- px[!(names(px) %in% c("cin", "cra", "cxy", "csi", "din", "page"))]
  #suppressWarnings(par(px))
  #return(summaries)
}

#tmp<-plot_shift_sum(shiftSummaries.list$F_uncex.fixed_rjMCMC)

#making 4 panel plot of shifting allometries
{
  pdf(file= paste("shift_summaries_full_custom_test.pdf", sep=""))    
  par(mfrow = c(2, 2))
  
  lapply(1:length(shiftSummaries.list), function(SS){ if(length(shiftSummaries.list[[SS]])>1){ plot_shift_sum(shiftSummaries.list[[SS]], lwd=2, single.plot=T, label.pts=F );
    mtext(names(shiftSummaries.list)[[SS]], side = 3, line = -3, outer = TRUE)}
    
  })
  
  dev.off()
}

dev.off()

#par(mfrow = c(2, 2))
#plot_shift_sum(shiftSummaries.list$FNN, lwd=2, single.plot=T, label.pts = F)

save.image('workspace.RDS')



cbind(shiftSummaries.list$NN$regressions, c(0,shiftSummaries.list$NN$PP))

# 
# 
# # MArginal likelihood estimation
# 
# 
# {
#   #running only 10 cores for a length.out=50 stepping stone sampler with 100k iterations. which means you will have 50 jobs of 100k gen mcmcs and 10 will be done in parallel at a time
#   # it is even better to run for more iterations (like 500k) but it will take a full day maybe more
#   registerDoParallel(cores=62)
#   
#   Bk <- qbeta(seq(0,1, length.out=100), 0.3,1)
#   #set.seed(1)
#   iterations=100000
#   
#   
#   #this is to compare the combined chains from each model tested, the mcmc.list object must be used to call 
#   #the stepping stone function, but we dont have a specific mcmc object for the combine chain so we use one of the mcmc object (mcmc.F, mcmc.NN, etc made at lines 315 to 325)
#   # we still provide the full combined chain for comparisson and it is still compare between chains just wante to make that clear to reduce confusion
#   
#   #additionall "mcmc.list[[3]]$`~/SSDI/birbs/nomiss_bmr_dat_mass_pred3_mcmc.FNN`" probably wont work for you, you need to chose what name you have for each mcmc object that you assigned above with " file_dir<-paste("~/SSDI/birbs/","nomiss_bmr_dat_mass_pred",i, sep="")"
#   
#   names(mcmc.list[[1]])
#   
#   
#   ss.F_uncex.merged.mtdnas  <-mcmc.list[[1]]$`1_F_uncex.merged.mtdnas`  $steppingstone( iterations, burnin.combined.chain.list$F_uncex.merged.mtdnas  , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.merged         <-mcmc.list[[1]]$`1_F_uncex.merged`         $steppingstone( iterations, burnin.combined.chain.list$F_uncex.merged         , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.fixed_rjMCMC     <-mcmc.list[[1]]$`1_F_uncex.fixed_rjMCMC`     $steppingstone( iterations, burnin.combined.chain.list$F_uncex.fixed_rjMCMC     , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.exons          <-mcmc.list[[1]]$`1_F_uncex.exons`          $steppingstone( iterations, burnin.combined.chain.list$F_uncex.exons          , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.introns        <-mcmc.list[[1]]$`1_F_uncex.introns`        $steppingstone( iterations, burnin.combined.chain.list$F_uncex.introns        , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.utrs           <-mcmc.list[[1]]$`1_F_uncex.utrs`           $steppingstone( iterations, burnin.combined.chain.list$F_uncex.utrs           , Bk, burnin=0.3, plot=FALSE)
#   ss.F_uncex.mtdnas.all     <-mcmc.list[[1]]$`1_F_uncex.mtdnas.all`     $steppingstone( iterations, burnin.combined.chain.list$F_uncex.mtdnas.all     , Bk, burnin=0.3, plot=FALSE)
#   ss.F_mtDNAs.proteins      <-mcmc.list[[1]]$`1_F_mtDNAs.proteins`      $steppingstone( iterations, burnin.combined.chain.list$F_mtDNAs.proteins      , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.merged.mtdnas<-mcmc.list[[1]]$`1_FNN_uncex.merged.mtdnas`$steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.merged.mtdnas, Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.merged       <-mcmc.list[[1]]$`1_FNN_uncex.merged`       $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.merged       , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.fixed_rjMCMC   <-mcmc.list[[1]]$`1_FNN_uncex.fixed_rjMCMC`   $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.fixed_rjMCMC   , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.exons        <-mcmc.list[[1]]$`1_FNN_uncex.exons`        $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.exons        , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.introns      <-mcmc.list[[1]]$`1_FNN_uncex.introns`      $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.introns      , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.utrs         <-mcmc.list[[1]]$`1_FNN_uncex.utrs`         $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.utrs         , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_uncex.mtdnas.all   <-mcmc.list[[1]]$`1_FNN_uncex.mtdnas.all`   $steppingstone( iterations, burnin.combined.chain.list$FNN_uncex.mtdnas.all   , Bk, burnin=0.3, plot=FALSE)
#   ss.FNN_mtDNAs.proteins    <-mcmc.list[[1]]$`1_FNN_mtDNAs.proteins`    $steppingstone( iterations, burnin.combined.chain.list$FNN_mtDNAs.proteins    , Bk, burnin=0.3, plot=FALSE)
#   ss.mcmc.11                <-mcmc.list[[1]]$`1_mcmc.11`                $steppingstone( iterations, burnin.combined.chain.list$global               , Bk, burnin=0.3, plot=FALSE)
#   ss.mcmc.N1                <-mcmc.list[[1]]$`1_mcmc.N1`                $steppingstone( iterations, burnin.combined.chain.list$N1             , Bk, burnin=0.3, plot=FALSE)
#   ss.mcmc.NN                <-mcmc.list[[1]]$`1_mcmc.NN`$steppingstone( iterations, burnin.combined.chain.list$NN , Bk, burnin=0.3, plot=FALSE)
#   
# 
#   
#   mlnL <- c("F_uncex.merged.mtdnas"  =ss.F_uncex.merged.mtdnas$lnr,
#             "F_uncex.merged"         =ss.F_uncex.merged$lnr,
#             "F_uncex.fixed_rjMCMC"     =ss.F_uncex.fixed_rjMCMC$lnr,
#             "F_uncex.exons"          =ss.F_uncex.exons$lnr,
#             "F_uncex.introns"        =ss.F_uncex.introns$lnr,
#             "F_uncex.utrs"           =ss.F_uncex.utrs$lnr,
#             "F_uncex.mtdnas.all"     =ss.F_uncex.mtdnas.all$lnr,
#             "F_mtDNAs.proteins"      =ss.F_mtDNAs.proteins$lnr,
#             "FNN_uncex.merged.mtdnas" =ss.FNN_uncex.merged.mtdnas$lnr,
#             "FNN_uncex.merged"       =ss.FNN_uncex.merged$lnr,
#             "FNN_uncex.fixed_rjMCMC"   =ss.FNN_uncex.fixed_rjMCMC$lnr,
#             "FNN_uncex.exons"        =ss.FNN_uncex.exons$lnr,
#             "FNN_uncex.introns"      =ss.FNN_uncex.introns$lnr,
#             "FNN_uncex.utrs"         =ss.FNN_uncex.utrs$lnr,
#             "FNN_uncex.mtdnas.all"   =ss.FNN_uncex.mtdnas.all$lnr,
#             "FNN_mtDNAs.proteins"    =ss.FNN_mtDNAs.proteins$lnr,
#             "mcmc.11"                =ss.mcmc.11$lnr,
#             "mcmc.N1"                =ss.mcmc.N1$lnr,
#             "mcmc.NN"                =ss.mcmc.NN$lnr
#   )
#   
#   
# }
# 
# 
# #check marginal likelihoods
# mlnL
# 
# saveRDS(mlnL, file="mlnL.RDS")
# #mlnL<-readRDS(file="mlnL.RDS")
# 
# as.data.frame(sort(mlnL))
# 
# 
# #extract the chain data #example
# as.data.frame(do.call(rbind, FNN_uncex.introns.combined.chains$theta))
# 
# #extract the chain data #example
# as.data.frame(do.call(rbind, NN.combined.chains$theta))[,2]
# 
# 
#write a function to extract the mean, median, and 95% HPDs for each param for each regime
table_generator<-function(data, output, HPD=0.95) {
  tmp<-list()
  params<-names(data)

  for( i in 1:length(params)){

    if(class(data[[i]]) == "numeric"){
      tmp[[i]]<-data[[i]]
    }

    if(class(data[[i]]) == "list") {
      tmp[[i]] <- as.data.frame(do.call(rbind, data[[i]]))
    }

  }

  names(tmp)<-params
  subset<-tmp[output]
  result<-list()

  for(i in 1:length(subset)){
    result[[i]]<-  cbind(as.data.frame(colMeans(subset[[i]])),
                         as.data.frame(apply(subset[[i]],2,median)),
                         as.data.frame(coda::HPDinterval(coda::mcmc(subset[[i]]), prob = HPD)))
    colnames(result[[i]])<-c('mean', 'median', 'lower', 'upper')
  }

  names(result)<-names(subset)

  return(result)

  }

test<-table_generator(data=FNN_uncex.merged.mtdnas.combined.chains,
                      output=c("theta", "beta_Pred"), HPD=0.95)


estimate_mode <- function(x, ...) {
  d <- density(x)
  mode<-d$x[which.max(d$y)]
  #int<-HDInterval::hdi(d)
  return(c(mode))
}


require(vioplot)


theta<-do.call(rbind, FNN_uncex.merged.mtdnas.combined.chains$theta)
#colnames(theta)<-c("Reg0", "Reg1", "Reg2", "Reg3", "Reg4", "Reg5", "Reg6", "Reg7", "Reg8", "Reg9", "Reg10", "Reg11", "Reg12")
colnames(theta)<-c("Reg0", "Rheiformes, Casuariiformes, and Apterygiformes", "Tinamiformes", "Notopaleognathae", "Aequornithes", "Coraciimorphae", "Passeri", "Psittaciformes", "Reminader of Neoaves", "Otidae", "Passerea", "Columbea", "Neognathae")
beta<-do.call(rbind, FNN_uncex.merged.mtdnas.combined.chains$beta_Pred)
#colnames(beta)<-c("Reg0", "Reg1", "Reg2", "Reg3", "Reg4", "Reg5", "Reg6", "Reg7", "Reg8", "Reg9", "Reg10", "Reg11", "Reg12")
colnames(beta)<-c("Reg0", "Rheiformes, Casuariiformes, and Apterygiformes", "Tinamiformes", "Notopaleognathae", "Aequornithes", "Coraciimorphae", "Passeri", "Psittaciformes", "Reminader of Neoaves", "Otidae", "Passerea", "Columbea", "Neognathae")

#testing 2d density plot (skip this section to generate violin plots, only for testing)
{
 theta<-as.data.frame(theta)
  theta<-tidyr::pivot_longer(theta, cols = everything())
  colnames(theta)<-c("Reg", "value")
  theta$Reg<-as.factor(theta$Reg)
  #theta<-dplyr::sample_n(theta, 10000)
  #saveRDS(theta,file="test.RDS")
  beta<-as.data.frame(beta)
  beta<-tidyr::pivot_longer(beta, cols = everything())
  colnames(beta)<-c("Reg", "value")
  beta$Reg<-as.factor(beta$Reg)
  
  
  data<-cbind(beta=beta,theta=theta) #%>%
    # group_by(beta.Reg) %>%
    # mutate(quant = quantile(beta.value, 0.45), row = row_number()) %>%
    # filter(beta.value >= quant) %>% 
    # mutate(quant = quantile(beta.value, 0.55), row = row_number()) %>%
    # filter(beta.value <= quant) %>%
    # 
    # group_by(theta.Reg) %>%
    # mutate(quant = quantile(theta.value, 0.45), row = row_number()) %>%
    # filter(theta.value >= quant) %>%
    # mutate(quant = quantile(theta.value, 0.55), row = row_number()) %>%
    # filter(theta.value <= quant)

  #data<- data %>% filter(theta.Reg == "Reg0")
  
  # Bin size control + color palette
  ggplot(data, aes(x=theta.value, y=beta.value) ) +
    #stat_dens2d_filter_g(keep.fraction = 1/4) +
    geom_hex(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    xlim(-4.2, -3) + 
    ylim(0.6, 0.85) +
    geom_density_2d(mapping=aes(x=theta.value, y=beta.value, groups=theta.Reg), inherit.aes = F)
    
  
  
  ggplot(theta, aes(x=Reg, y=value)) +
    geom_boxplot(fill='#A4A4A4', color="black")+
    theme_classic()
}
#beta[,-c(4,11)]

#violin plots
pdf(file="BMR_vioplots.pdf", height=7, width=6)
par(mfrow=c(2,1), mar=c(0, 4.1, 8, 2.1))
vioplot(beta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)], las=2, names=NA, main="test", ylab="slope", ylim=c(0.45, 1))
#add prior
abline(h=0.70, lty=2)
#add maximum
abline(h=max(matrixStats::colMedians(beta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)])), lty=3)
#add minimum
abline(h=min(matrixStats::colMedians(beta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)])), lty=3)

par(mar=c(7, 4.1, 1, 2.1))
vioplot(theta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)], las=2, ylab="intercept", ylim=c(-6, -1.75))
#add prior
abline(h=-3.5, lty=2)
#add maximum
abline(h=max(matrixStats::colMedians(theta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)])), lty=3)
#add minimum
abline(h=min(matrixStats::colMedians(theta[,c(1, 3, 2, 13, 12, 10, 9, 5, 6, 8, 7)])), lty=3)

dev.off()


#vioplot(beta[,-c(4,11)][,order(matrixStats::colMedians(beta[,-c(4,11)]))], las=2)
#vioplot(theta[,-c(4,11)][,order(matrixStats::colMedians(beta[,-c(4,11)]))], las=2)


par(mfrow=c(1,2))
plotrix::plotCI(x=seq(1:length(test$theta$mean)),y=test$theta$mean, ui=test$theta$upper,
                li=test$theta$lower, pch=19, ylab="theta", xlab="regime")
plotrix::plotCI(x=seq(1:length(test$theta$mean)),y=test$beta_Pred$mean, ui=test$beta_Pred$upper,
                li=test$beta_Pred$lower, pch=19, ylab="betaPred", xlab="regime")


dev.off()



