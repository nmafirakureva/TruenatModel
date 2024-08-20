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
PDCO <- PDCO |> 
  # filter(grepl('att',parm)) |>
  mutate(parm = case_when(
    age=='0-4'& grepl('att',parm) ~ paste0('u5.',parm),
    age=='5-14'& grepl('att',parm) ~ paste0('o5.',parm),
    .default = parm))

unique(PDCO$parm)
PDCO |> 
  filter(grepl('att',parm)) 

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
D[['soc_notx_check']] <- SOC.notx.F$checkfun(D)
D[['soc_notx_cost']] <- SOC.notx.F$costfun(D)

D[['int_att_check']] <- INT.att.F$checkfun(D)
D[['int_att_cost']] <- INT.att.F$costfun(D)
D[['int_notx_check']] <- INT.notx.F$checkfun(D)
D[['int_notx_cost']] <- INT.notx.F$costfun(D)

## cross checks: compute probability of endpoints from whole tree vs pruned trees
head(D[,int_att_check])
head(D[, attend.int])

## NOTE OK
all(D[, attend.int] == D[, int_att_check])
all(D[, attend.soc] == D[, soc_att_check])

## NOTE confusingly there is also a notxend variable, which is different
all(D[, notx.int] == D[, int_notx_check])
all(D[, notx.soc] == D[, soc_notx_check])

D[,table(tb)]

## create restricted PSA
DR <- D[,.(id,tb,
           soc_att_check,
           soc_att_cost,
           soc_notx_check,
           soc_notx_cost,
           int_att_check,
           int_att_cost,
           int_notx_check,
           int_notx_cost
           )]

## condition costs on outcome
DR[,c('soc_att_cost',
      'soc_notx_cost',
      'int_att_cost',
      'int_notx_cost'):=.(
        soc_att_cost/soc_att_check,
        soc_notx_cost/soc_notx_check,
        int_att_cost/int_att_check,
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

(dxb <- names(D)[grepl('dxb|dxc|att.soc|att.int',names(D))])
summary(D[tb=='TB',..dxb])
summary(D[tb=='noTB',..dxb])

## Cost-effectiveness analysis
## --- run over different countries
cnmz <- names(C)
cnmz <- cnmz[cnmz!='id']

names(labz)[grepl('cost$',names(labz))]
costsbystg <- c('DH.presumptive.cost.soc','DH.evaluated.cost.soc','DH.treated.cost.soc',
                'DH.presumptive.cost.int','DH.evaluated.cost.int','DH.treated.cost.int',
                'PHC.presumptive.cost.soc','PHC.evaluated.cost.soc','PHC.treated.cost.soc',
                'PHC.presumptive.cost.int','PHC.evaluated.cost.int','PHC.treated.cost.int')
toget <- c('id',
           'cost.soc','cost.int',
           costsbystg,
           'att.soc','att.int',
           'deaths.soc','deaths.int',
           'LYS','LYS0','value'
)
notwt <- c('id','LYS','LYS0','value') #variables not to weight against value
lyarm <- c('LYL.soc','LYL.int')
lyarm <- c(lyarm,gsub('\\.','0\\.',lyarm)) #include undiscounted
tosum <- c(setdiff(toget,notwt),lyarm)
## heuristic to scale top value for thresholds:
heur <- c('id','value','deaths.int','deaths.soc')
out <- D[,..heur]
out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=c('deaths.int','deaths.soc'),by=id] #sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.iph)]
topl <- 500 #100
lz <- seq(from = 0,to=topl,length.out = 1000) #threshold vector for CEACs
## staged costs by arm
soc.sc <- grep('soc',costsbystg,value=TRUE); psoc.sc <- paste0('perATT.',soc.sc)
int.sc <- grep('int',costsbystg,value=TRUE); pint.sc <- paste0('perATT.',int.sc)


## containers & loop
allout <- allpout <- allscout <- list() #tabular outputs
ceacl <- NMB <- list()             #CEAC outputs etc

## cn <- isoz[1]
for(cn in isoz){
  cat('running model for:',cn,'\n')
  ## --- costs
  ## drop previous costs
  D[,c(cnmz):=NULL]
  ## add cost data
  C <- MakeCostData(allcosts,nreps) #make cost PSA
  D <- merge(D,C,by='id',all.x = TRUE)        #merge into PSA
  ## --- DALYs
  ## drop any that are there
  if('LYS' %in% names(D)) D[,c('LYS','LYS0'):=NULL]
  D <- merge(D,LYKc[,.(age,LYS,LYS0)],by='age',all.x = TRUE)        #merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(D <- runallfuns(D,arm=arms)))
  ## --- grather outcomes
  out <- D[,..toget]
  out[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
                   LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  out <- out[,lapply(.SD,function(x) sum(x*value)),.SDcols=tosum,by=id] #sum against popn
  ## non-incremental cost per ATT
  out[,costperATT.soc:=cost.soc/att.soc];
  out[,(psoc.sc):=lapply(.SD,function(x) x/att.soc),.SDcols=soc.sc]
  out[,costperATT.int:=cost.int/att.int];
  out[,(pint.sc):=lapply(.SD,function(x) x/att.int),.SDcols=int.sc]
  ## increments wrt SOC (per child presenting at either DH/PHC)
  out[,Dcost.int:=cost.int-cost.soc] #inc costs
  out[,Datt.int:=att.int-att.soc] #inc atts
  out[,attPC.int:=1e2*att.int/att.soc] #rel inc atts
  out[,Ddeaths.int:=deaths.int-deaths.soc] #inc deaths
  out[,DLYL0.int:=LYL0.int-LYL0.soc] #inc LYLs w/o discount
  out[,DLYL.int:=LYL.int-LYL.soc] #inc LYLs
  ## per whatever
  out[,DcostperATT.int:=cost.int/att.int-cost.soc/att.soc];
  out[,Dcostperdeaths.int:=-cost.int/deaths.int+cost.soc/deaths.soc]
  out[,DcostperLYS0.int:=-cost.int/LYL0.int+cost.soc/LYL0.soc]
  out[,DcostperLYS.int:=-cost.int/LYL.int+cost.soc/LYL.soc]
  ## D/D
  out[,DcostperDATT.int:=Dcost.int/Datt.int];
  out[,DcostperDdeaths.int:=-Dcost.int/Ddeaths.int]
  out[,DcostperDLYS0.int:=-Dcost.int/DLYL0.int]
  out[,DcostperDLYS.int:=-Dcost.int/DLYL.int]
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs; pouts <- smy$pouts; scouts <- smy$scouts
  outs[,iso3:=cn]; pouts[,iso3:=cn]; scouts[,iso3:=cn]
  ## capture tabular
  allout[[cn]] <- outs; allpout[[cn]] <- pouts; allscout[[cn]] <- scouts
  ## capture data for NMB
  NMB[[cn]] <- out[,.(iso3=cn,DLYL.int,Dcost.int)]
  ## ceac data
  ceacl[[cn]] <- data.table(iso3=cn,
                            int=make.ceac(out[,.(Q=-DLYL.int,P=Dcost.int)],lz),
                            threshold=lz)
}

allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
allscout <- rbindlist(allscout)
ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)

fwrite(allout,file=gh('outdata/allout') + SA + '.csv')
fwrite(allpout,file=gh('outdata/allpout') + SA + '.csv')
fwrite(allscout,file=gh('outdata/allscout') + SA + '.csv')
save(ceacl,file=gh('outdata/ceacl') + SA + '.Rdata')
save(NMB,file=gh('outdata/NMB') + SA + '.Rdata')

## checks
out[,.(att.int/att.soc)]
out[,.(Datt.int/att.soc)]
out[,.(att.soc,att.int)]
out[,.(Ddeaths.int/Datt.int)] #OK
out[,.(DLYL.int/Ddeaths.int)] #OK
## check
allout[,.(costperATT.int.mid-costperATT.soc.mid,DcostperATT.int.mid)]

## CEAC plot
cbPalette <- c("#009E73")
ceaclm <- melt(ceacl,id=c('iso3','threshold'))
ceaclm[,Intervention:=ifelse(variable=='int','Intervention','Standard of care')]
## name key
ckey <- data.table(iso3=c('NGA'),
                   country=c('Nigeria'))

ceaclm <- merge(ceaclm,ckey,by='iso3',all.x=TRUE)

## plot: int only
GP <- ggplot(ceaclm[variable=='int' &
                      iso3 %in% c('NGA')],
             aes(threshold,value,
                 col=country,lty=Intervention)) +
  geom_line(show.legend = FALSE) +
  theme_classic() +
  theme(legend.position = 'top',legend.title = element_blank())+
  ggpubr::grids()+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY averted)')+
  scale_colour_manual(values=cbPalette)
GP

ggsave(GP,file=gh('plots/CEAC') + SA + '.png',w=7,h=5)

allout$ICER.int

