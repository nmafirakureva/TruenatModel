## Pete version
## quick run through the PPD pathways tree

## flags for sensitivity analyses
shell <- FALSE # whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  # none,base/lo/hi,tptru,hicoprev
  if(SA == 'none'){
    SA <- ''
  } 
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'cdr' = making cdr higher for incidence
  ## 'txd' = making the completion influence tx/pt outcome
  sacases <- c('','lo','tptru','hicoprev', 'ctryeff','ugaattcsts', 'cdr', 'hivprev')
  SA <- sacases[1]
}

# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
source(here('R/truenat_pathways_tree1.R'))           #tree structure and namings: also tree functions & libraries
source(here('R/truenat_functions.R'))      #functions for tree parameters

## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
hivlevels <- c(0,1)
artlevels <- c(0,1)
tblevels <- c('TB', 'noTB') # Active TB disease, no TB
agelevels <- c('0-4','5-14')

isoz <- c('NGA') #relevant countries

## --- life years and other outputs NOTE needs to be set FALSE on first run thru
LYSdone <- FALSE
if(!LYSdone){
  ## make discounted life-years if they haven't been done
  LYKc <- GetLifeYears(isolist=isoz,discount.rate=0.03,yearfrom=2021)
  LYKc0 <- GetLifeYears(isolist=isoz,discount.rate=0.00,yearfrom=2021)
  LYKc5 <- GetLifeYears(isolist=isoz,discount.rate=0.05,yearfrom=2021)
  LYKc <- merge(LYKc,LYKc0[,.(iso3,age,LYS0=LYS)],by=c('iso3','age'))
  LYKc <- merge(LYKc,LYKc5[,.(iso3,age,LYS5=LYS)],by=c('iso3','age'))
  LYK <- LYKc[,.(LYS=mean(LYS),LYS0=mean(LYS0),LYS5=mean(LYS5)),by=.(age)] #averaged life-years 4 generic tests
  save(LYKc,file=here('indata/LYKc.Rdata'))
  save(LYK,file=here('indata/LYK.Rdata'))
} else {
  load(file=here('indata/LYKc.Rdata'))
  load(file=here('indata/LYK.Rdata'))
}

# # Sensitivity analysis: 0% & 5% discount rates
# if(SA %in% c('hi','lo')){
#   LYKc[,LYS:=ifelse(SA=='lo', LYS0, 
#                     ifelse(SA=='hi',LYS5, LYS))]  
# }

## prior parameters
PD <- read.csv(here('indata/parms.csv')) #read in main parameters
PDI <- read.csv(here('indata/CASCI.csv')) #read in presentation parameters
PDA <- read.csv(here('indata/CASCA.csv')) #read in age splits 
PDC <- read.csv(here('indata/CASCP.csv')) #read in cascade parameters
PDO <- read.csv(here('indata/CASCPO.csv')) #read in outcome parameters
# AD <- read.csv(here('indata/DiagnosticAccuracy.csv')) #read in accuracy parameters
# RD <- fread(gh('indata/RUParms.csv'))    #read resource use data
CD <- fread(gh('indata/costs.csv'))    #read cost data

PDCO <- rbind(PDC,PDO, PDI, PDA) # combine cascade and outcome parameters
names(PD)
names(PDCO)

# Quick check & fix
unique(PDCO$parm)
PDCO$parm <- gsub('prusumed', 'presumed', PDCO$parm)

## parameters to be determined from cascade data
PD0 <- PD |> 
  mutate(DISTRIBUTION = as.numeric(DISTRIBUTION)) |> 
  filter(!is.na(DISTRIBUTION)) 

## the rest
PD1 <- PD |> 
  mutate(dist = as.numeric(DISTRIBUTION)) |> 
  filter(is.na(dist)) |>
  select(-dist)

names(PDCO)
names(PD0)

unique(PDCO$parm)
PD2 <- rbind(
  PDCO |> 
    rename(NAME = parm, DISTRIBUTION = frac) |> 
    filter(!grepl('clinbac.assess|TrueNat|TBLamp', NAME)) |> # dropping things not needed for now
    select(NAME, DISTRIBUTION),
  PD0 |> 
    select(NAME, DISTRIBUTION)) |> 
  distinct(NAME, .keep_all = TRUE) |> 
  pivot_wider(names_from = NAME, values_from = DISTRIBUTION)

names(PD2)

# convert into parameter object
P <- parse.parmtable(PD1)             

## make base PSA dataset
set.seed(1234) #random number seed

D0 <- makePSA(nreps,P,dbls = list(c('cfrhivor','cfrartor')))


# merge in fixed parameter
names(PD2)
D0 <- cbind(D0, PD2)

## use these parameters to construct input data by attribute
D0 <- makeAttributes(D0)
D0[,sum(value),by=.(id, tb)] #CHECK

## read and make cost data
rcsts <- CD

names(CD)

rcsts <- setDT(rcsts)

## turn cost data into PSA
rcsts[is.na(rcsts)] <- 0 # some quick fix >> setting NA to 0
rcsts[cost.sd==0,cost.sd:=cost.m/40]        #SD such that 95% UI ~ 10% of mean

allcosts <- rcsts[,.(cost=parm, cost.m, cost.sd)]

C <- MakeCostData(allcosts,nreps)               # make cost PSA NOTE using CMR cost data

## NOTE can re-run from here to implement changes to MakeTreeParms
## add cost data
D <- merge(D0,C,by=c('id'),all.x=TRUE)       # merge into PSA (differentiated D and D0 to facilitate rerunning)

## compute other parameters (adds by side-effect)
MakeTreeParms(D,P)


# names(D)[grepl('cost', names(D))]

## checks
D[,sum(value),by=.(id)] #CHECK
D[,sum(value),by=.(id, age)] #CHECK
D[,sum(value),by=.(id,age, tb)] #CHECK


## check for leaks
head(SOC.F$checkfun(D)) #SOC arm
head(INT.F$checkfun(D)) #INT arm

names(SOC.F)

## === RUN MODEL
arms <- c('SOC','INT')
D <- runallfuns(D,arm=arms)                      #appends anwers

## restricted trees:
D[['soc_att_check']] <- SOC.att.F$checkfun(D)
D[['soc_att_cost']] <- SOC.att.F$costfun(D)
D[['soc_att_ppd']] <- SOC.att.F$cost.ppdfun(D)
D[['soc_att_nhs']] <- SOC.att.F$cost.nhsfun(D)
D[['soc_tpt_check']] <- SOC.tpt.F$checkfun(D)
D[['soc_tpt_cost']] <- SOC.tpt.F$costfun(D)
D[['soc_tpt_ppd']] <- SOC.tpt.F$cost.ppdfun(D)
D[['soc_tpt_nhs']] <- SOC.tpt.F$cost.nhsfun(D)
D[['soc_notx_check']] <- SOC.notx.F$checkfun(D)
D[['soc_notx_cost']] <- SOC.notx.F$costfun(D)

D[['int_att_check']] <- INT.att.F$checkfun(D)
D[['int_att_cost']] <- INT.att.F$costfun(D)
D[['int_att_ppd']] <- INT.att.F$cost.ppdfun(D)
D[['int_att_nhs']] <- INT.att.F$cost.nhsfun(D)
D[['int_tpt_check']] <- INT.tpt.F$checkfun(D)
D[['int_tpt_cost']] <- INT.tpt.F$costfun(D)
D[['int_tpt_ppd']] <- INT.tpt.F$cost.ppdfun(D)
D[['int_tpt_nhs']] <- INT.att.F$cost.nhsfun(D)
D[['int_notx_check']] <- INT.notx.F$checkfun(D)
D[['int_notx_cost']] <- INT.notx.F$costfun(D)

## cross checks: compute probability of endpoints from whole tree vs pruned trees
head(D[,int_att_check])
head(D[, attend.int])

## NOTE OK
all(D[, attend.int] == D[, int_att_check])
all(D[, attend.soc] == D[, soc_att_check])
all(D[, tptend.int] == D[, int_tpt_check])
all(D[, tptend.soc] == D[, soc_tpt_check])
## NOTE confusingly there is also a notxend variable, which is different
all(D[, notx.int] == D[, int_notx_check])
all(D[, notx.soc] == D[, soc_notx_check])


D[,table(tb)]

## create restricted PSA
DR <- D[,.(id,tb,
           soc_att_check,
           soc_att_cost,
           soc_att_ppdcost=soc_att_ppd/soc_att_cost,
           soc_tpt_check,
           soc_tpt_cost,
           soc_tpt_ppdcost=soc_tpt_ppd/soc_tpt_cost,
           soc_notx_check,
           soc_notx_cost,
           int_att_check,
           int_att_cost,
           int_att_ppdcost=int_att_ppd/int_att_cost,
           int_tpt_check,
           int_tpt_cost,
           int_tpt_ppdcost=int_tpt_ppd/int_tpt_cost,
           int_notx_check,
           int_notx_cost
           )]

## condition costs on outcome
DR[,c('soc_att_cost',
      'soc_tpt_cost',
      'soc_notx_cost',
      'int_att_cost',
      'int_tpt_cost',
      'int_notx_cost'):=.(
        soc_att_cost/soc_att_check,
        soc_tpt_cost/soc_tpt_check,
        soc_notx_cost/soc_notx_check,
        int_att_cost/int_att_check,
        int_tpt_cost/int_tpt_check,
        int_notx_cost/int_notx_check
      )]

save(DR,file=here('outdata/DR.Rdata'))

## summary
DRS <- DR[,lapply(.SD,mean),.SDcols=names(DR)[-c(1,2)],by=tb]
DRS <- melt(DRS,id='tb')
DRS[,c('arm','outcome','quantity'):=tstrsplit(variable,split='_')]
options(scipen=999)
(DRS <- dcast(data=DRS,formula=arm+quantity+outcome~tb,value.var='value'))
fwrite(DRS,file=here('outdata/DRS.csv'))

# D[tb=='noTB',.(soc.prop.prev.tb.dx,
#     soc.prop.no.prev.tb.dx.symp,
#     soc.prop.no.prev.tb.dx.symp.gp.assess,
#     soc.prop.no.prev.tb.dx.symp.tb.suspicion,
#     soc.prop.no.prev.tb.dx.symp.nhs.referral,
#     soc.prop.no.prev.tb.dx.symp.tb.dx,
#     soc.prop.starting.att,
#     soc.prop.completing.att,
#     sens.any.abn.xray,spec.any.abn.xray)]

## soc.prop.no.prev.tb.dx.symp.tb.dx??

## cost of getting ATT (from CSV output)
D[,mean((pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + 
           (1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visit + cost.prison.escort)) + 
          (pDSTB*DurDSTB*cost.dsatt.drugs + (1-pDSTB)*DurMDRTB*cost.mdratt.drugs) + cost.dots*(pDSTB*DurDSTB + (1-pDSTB)*DurMDRTB) +
          cost.inpatient),
  by=tb]
D[,mean(IncompDurDSTB/DurDSTB*(pDSTB*dstb.visits*(cost.dstb.opd.visit + cost.prison.escort) + pDSTB*DurDSTB*(cost.dsatt.drugs + cost.dots)) + 
          IncompDurMDRTB/DurMDRTB*((1-pDSTB)*mdrtb.visits*(cost.mdrtb.opd.visits + cost.prison.escort)  + (1-pDSTB)*DurMDRTB*(cost.mdratt.drugs + cost.dots)) + 
          cost.inpatient
),by=tb]
D[,mean(durTPT*(cost.ltbi.drugs + cost.dots)),by=tb]
D[,mean((durTPT*(cost.ltbi.drugs + cost.dots) + TPT.visits*(cost.tpt.opd.visit + cost.prison.escort))),by=tb]
D[,mean((IncompDurTPT*(durTPT*(cost.ltbi.drugs + cost.dots) + TPT.visits*(cost.tpt.opd.visit + cost.prison.escort)))),by=tb] # BUG fixed
D[,mean((IncompDurTPT/durTPT*(durTPT*(cost.ltbi.drugs + cost.dots) + TPT.visits*(cost.tpt.opd.visit + cost.prison.escort)))),by=tb]
