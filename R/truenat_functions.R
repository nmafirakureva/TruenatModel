## TODO questions:
##
## - stool/sputum
## - flag assumption = groups for SA or pending data
##

## ========= UTILITIES ===============
logit <- function(x) log(odds(x))
ilogit <- function(x) iodds(exp(x))
AOR <- function(base,OR) ilogit(log(OR) + logit(base))
AOR2 <- function(base,OR) OR*base/(1-base+(base*OR))
odds <- function(x) x/(1-x)
iodds <- function(x) x/(1+x)
lo <- function(x) quantile(x,probs = 0.025, na.rm=TRUE)
hi <- function(x) quantile(x,probs = 1-0.025, na.rm=TRUE)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

brkt <- function(M,L,H,ndp=0) paste0(round(M,ndp),' (',
                                     round(L,ndp),' - ',
                                     round(H,ndp),')')
gm <- function(x) exp(mean(log(x))) #geometric mean
gh <- function(x) glue(here(x))

## ========= OUTCOMES ===============
## TODO - remove excess RNG here
CFRtxY <- function(age,hiv=0,art=0){#NB optimized for clarity not speed
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$ontxY$r(length(age))
  tmp[age=='5-14'] <- P$ontxO$r(sum(age=='5-14'))  #NB this could be achieved in the tree model
  ## hivartOR
  Z <- P$hivartOR$r(length(age))
  hor <- rep(1,length(age))
  tmp <- logit(tmp)                     #transform
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
## CFRtxY(1:10,P)                            #test
## summary(CFRtxY(1:1e3,P))
## summary(CFRtxY(1:1e3,hiv=1))
## summary(CFRtxY(1:1e3,hiv=1,art=1,P))


## == CFR off tx
CFRtxN <- function(age,hiv=0,art=0){
  if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
  if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
  tmp <- P$notxY$r(length(age))          #default a<5 and hiv=art=0
  tmp[age!='5-14' & hiv>0 & art==0] <- P$notxY$r(sum(age!='5-14' & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[age!='5-14' & hiv>0 & art>0] <- P$notxHAY$r(sum(age!='5-14' & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[age=='5-14'] <- P$notxO$r(sum(age=='5-14'))    #o5, HIV-ve
  tmp[age=='5-14' & hiv>0 & art==0] <- P$notxHAO$r(sum(age=='5-14' & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[age=='5-14' & hiv>0 & art>0] <- P$notxHAO$r(sum(age=='5-14' & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
## CFRtxN(1:10,P)                            #test
## summary(CFRtxN(1:1e3,P))
## summary(CFRtxN(1:1e3,hiv=1,P))
## summary(CFRtxN(1:1e3,hiv=1,art=1,P))

## add CFRs to data by side-effect
AddCFRs <- function(D,P){
  ## d.cfr.notx & d.cfr.tx
  D[,c('cfr.notx','cfr.tx'):=0] #NOTE neglect non-TB mortality
  ## CFR on  ATT
  D[,cfr.tx:=CFRtxY(age,hiv,art)]
  ## CFR w/o ATT
  D[,cfr.notx:=CFRtxN(age,hiv,art)]
}


## new parameters as part of reach work
AddDetectionLabels <- function(D){
  # D[,CDR:=ifelse(isoz=='CMR', 0.1878, 0.6190)] # fixed to mean background TB detection
  # D[,CDRi:=CDR*1.5] # TB detection if household visited - upscale background mean by a factor of 1.5
  D[,CDR:=cdr] #TODO country-specific background TB detection sampled from a beta distribution
  # D[,CDRi:=pmin(cdr*1.5,1)] #TODO TB detection if household visited based on adjusted country-specific background TB detection
  D[,CDRi:=cdri] #TODO TB detection if household visited based on adjusted country-specific background TB detection
  tmp <- eval(parse(text=INTtbprev),envir=D) #TODO unclear why 0?!
  tmp <- rep(0.03,length(tmp))               #TODO for testing - made up
  # D[,int.tbprev.symptomatic:=tmp] #TB prev in symptomatics, based on INT
  # D[,soc.tbprev.symptomatic := tmp]
  # D[,soc.tbprev.symptomatic := int.tbprev.symptomatic]
  D[,soc.tbprev.symptomatic := int.tbprev/soc.frac.symp]
  D[,int.tbprev.symptomatic := int.tbprev/int.frac.symp]
}

## == case detection
aCDR <- function(mn,ab){
  mn <- mn*(1 + runif(length(mn))) #CDR adjustment 2
  mn <- pmin(mn,1)
  a <- mn*ab
  b <- (1-mn)*ab
  rbeta(n=length(mn),shape1 = a,shape2 = b)
  ## 0.4
}

## test <- eval(parse(text=INTtbprev),envir=D)
## tail(test)

## ========= DIAGNOSIS ===============

## function for combining sample modality with bacteriological test
## to calculate the probablility bac+ TB is diagnosed
TBbacsampletest <- function(samplepossible,testpos){
  samplepossible * testpos
}

## function to add in the Sample/Test combined probabilities of TB dx
## works by side effect
AddSampleTests <- function(D){

  ## Bacteriological test availability
  # D[,soc.dh.test:=rbeta(nrow(D),8.251046,5.500698)]
  # D[,int.dh.test:=rbeta(nrow(D),8.251046,5.500698)]
  # D[,soc.dh.test:=1]
  # D[,int.dh.test:=1]

  # PHC
  # D[,soc.phc.test:=rbeta(nrow(D),5.335947,172.5289)]
  # D[,int.phc.test:=rbeta(nrow(D),8.251046,5.500698)]
  # D[,int.phc.test:=0.6] # to increase this with Truenat introduction
  # summary(D$int.phc.test)

  # Type of test available
  D[,soc.dh.fracUltra:=ifelse(age=='5-14',1,1)] # assume only Xpert is available in DH
  D[,int.dh.fracUltra:=ifelse(age=='5-14',1,1)]
  summary(D$soc.dh.fracUltra)

  D[,soc.phc.fracUltra:=ifelse(age=='5-14',0,0)] # assume no Xpert is not available in PHC
  D[,int.phc.fracUltra:=ifelse(age=='5-14',0,0)]

  # Bacteriological test possibility
  # Test is only possible for children who can provide a sample
  D[,soc.dh.test:=ifelse(age=='5-14',soc.dh.test*soc.dh.exp.sputum.o5,soc.dh.test*soc.dh.exp.sputum.u5)] # Actual children tested
  D[,int.dh.test:=ifelse(age=='5-14',int.dh.test*soc.dh.exp.sputum.o5,int.dh.test*soc.dh.exp.sputum.u5)]

  D[,soc.phc.fracUltra:=soc.phc.fracUltra*ifelse(age=='5-14',soc.phc.exp.sputum.o5,soc.phc.exp.sputum.u5)] # Actual children tested
  D[,int.phc.fracUltra:=int.phc.fracUltra*ifelse(age=='5-14',soc.phc.exp.sputum.o5,soc.phc.exp.sputum.u5)]

  ## ------- X on sputum/GA -------
  ## People receiving Xpert Ultra testing [either sputum or GA], in those identified as having presumptive TB

  ## TB dx bac+ on Xpert Ultra on sputum
  D[,soc.dh.ptbxUltra:=ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum)]
  D[,int.dh.ptbxUltra:=ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum)]

  D[,soc.phc.ptbxUltra:=ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum)]
  D[,int.phc.ptbxUltra:=ifelse(tb=='TB+',sens.xsputum,1-spec.xsputum)]

  # TB dx bac+ on TrueNat on stool
  D[,soc.dh.ptbxTrueNat:=ifelse(tb=='TB+',sens.truenat.sputum,1-spec.truenat.sputum)]
  D[,int.dh.ptbxTrueNat:=ifelse(tb=='TB+',sens.truenat.sputum,1-spec.truenat.sputum)]

  D[,soc.phc.ptbxTrueNat:=ifelse(tb=='TB+',sens.truenat.sputum,1-spec.truenat.sputum)]
  D[,int.phc.ptbxTrueNat:=ifelse(tb=='TB+',sens.truenat.sputum,1-spec.truenat.sputum)]

  ## ------- clinical -------
  ## TB dx clinical, in bac- people identified as having presumptive TB
  D[,soc.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,soc.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

  D[,int.dh.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,int.phc.test.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

  ## TB dx clinical, in bac- people identified as having presumptive TB
  D[,soc.dh.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,soc.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

  D[,int.dh.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,int.phc.notest.ptbc:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

  # reassessment
  # D[,soc.dh.bact.tbdx:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,soc.dh.bact.tbdx:=ifelse(tb!='noTB',sens.clin,1-spec.clin)] #Assume no bacteriological testing @ reassessment
  D[,int.dh.bact.tbdx:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

  D[,soc.phc.bact.tbdx:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]
  D[,int.phc.bact.tbdx:=ifelse(tb!='noTB',sens.clin,1-spec.clin)]

}

## ======= COMBINED LABELLER ===========

## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D){

  # generate some parameters
  (check <- names(D)[grep('presumptive|presented',names(D))])
  summary(D[,..check])

  # Initial presentation

  # TB SPEED Decentralization assumed 90 (80-100) # https://doi.org/10.1016/j.eclinm.2024.102528
  # (p <- getAB(0.90, ((80-100)^2)/392^2))
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  # D[,phc.presented :=rbeta(nrow(D),30.21696,3.35744)]
  # D[,phc.presented :=0.75]

  ## TB more likely to go to DH
  D[tb!='noTB' & age=='0-4',phc.presented:=1-iodds(OR.dh.if.TB.u5*odds(1-phc.presented))]
  D[tb!='noTB' & age!='0-4',phc.presented:=1-iodds(OR.dh.if.TB.o5*odds(1-phc.presented))]

  # summary(D[,.(phc.presented),])
  # (D[,.(mean(phc.presented)),by=tb])
  D[,dh.presented :=1-phc.presented]

  # presumptive TB
  # (p <- getAB(0.89, ((98-60)^2)/392^2))  #sensitivity:89% (95% CI 52% to 98%) based on https://doi.org/10.1002/14651858.CD013693.pub2
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  # D[,sens.sympt.screen:=rbeta(nrow(D),8.38209,1.035989)]

  # (p <- getAB(0.69, ((51-83)^2)/392^2))  #specificity: 69% (95% CI 51% to 83%) based on https://doi.org/10.1002/14651858.CD013693.pub2
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  # D[,spec.sympt.screen:=rbeta(nrow(D),21.45787,9.640494)]

  # PD2 <- setDT(PD2)
  # summary(PD2[,dh.presumed/phc.presumed])

  D[,fac:=dh.presumed/phc.presumed] # factor to scale down sensitivity of presuming at PHC
  ## specificity of presuming
  D[,dh.presumed:=ifelse(tb!='noTB',1,1-spec.sympt.screen)] #TODO:placeholder for now
  D[,phc.presumed:=ifelse(tb!='noTB',1,1-spec.sympt.screen)] #TODO:placeholder for now

  # # TODO: check if this is OKAY
  # # Approach to normalizing everything to presumptive TB
  # D[,phc.presented:=phc.presented/phc.presumed]
  # D[,dh.presented:=dh.presented/dh.presumed]
  # D[,phc.presumed:=phc.presumed/phc.presumed]
  # D[,dh.presumed:=dh.presumed/dh.presumed]
  #
  # D |> select(tb, phc.presented,dh.presented,phc.presumed,dh.presumed) |>
  #   group_by(tb) |> summarise_all(mean)
  # hospital referral loss to follow-up
  # presumptive TB
  # (p <- getAB(0.5, ((0-80)^2)/392^2))
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  # names(PD)
  # PD |> filter(NAME=='soc.phc.rltfu')
  # summary(D$soc.phc.rltfu)
  # D[,soc.phc.rltfu:=1]
  D[,int.phc.rltfu:=soc.phc.rltfu]
  summary(D$int.phc.rltfu)

  # RIF resistance
  # summary(D[,dh.prr]*100)
  # summary(D[,phc.prr]*100)
  # 26/2414
  D[,prr:=dh.prr]

  # Bacteriological confirmation
  # D[,soc.dh.bact.tbdx:=1-dh.att]
  # D[,int.dh.bact.tbdx:=1-phc.att]

  # pre-treatment loss to follow-up
  # (check <- names(D)[grep('u5.dh.att|o5.dh.att|u5.phc.att|o5.phc.att',names(D))])
  # summary(D[,..check])
  D[,dh.ptltfu:=ifelse(age=='0-4',1-u5.dh.att,1-o5.dh.att)]
  D[,phc.ptltfu:=ifelse(age=='0-4',1-u5.phc.att,1-o5.phc.att)]
  summary(D$dh.ptltfu)

  D[,phc.14phcltfu:=dh.14dhltfu]
  #
  D[,hivprev.u5:=frac.hiv]
  D[,hivprev.o5:=frac.hiv]

  # costs
  D[,cost.phc.rsATT:=ifelse(age=='0-4',cost.phc.rsATT.04,cost.phc.rsATT.514)]
  D[,cost.dh.rsATT:=ifelse(age=='0-4',cost.dh.rsATT.04,cost.dh.rsATT.514)]

  D[,cost.phc.rrATT:=ifelse(age=='0-4',cost.phc.rrATT.04,cost.phc.rrATT.514)]
  D[,cost.dh.rrATT:=ifelse(age=='0-4',cost.dh.rrATT.04,cost.dh.rrATT.514)]

  # TODO: need % receiving X-ray, % samples referred for testing by facility too add these costs

  # X-ray
  D[,cost.phc.evaluation:=cost.phc.evaluation+ifelse(age=='0-4',u5.phc.p.xray*cost.cxr.exam,o5.phc.p.xray*cost.cxr.exam)]
  D[,cost.dh.evaluation:=cost.dh.evaluation+ifelse(age=='0-4',u5.dh.p.xray*cost.cxr.exam,o5.dh.p.xray*cost.cxr.exam)]

  # Currently assuming no cost of patient referral
  D[,cost.phc.refer:=0]
  # TODO: Check if ATT initiation and follow up costs are included

}

## combined function to add the labels to the tree prior to calculations
MakeTreeParms <- function(D,P){
  ## -- use of other functions
  AddSampleTests(D) #samples/tests
  AddCFRs(D,P) #outcomes
  # AddTPTrr(D,P)
  # AddProgProb(D, P)
  ## new labels from data
  AddDataDrivenLabels(D)
  # AddDetectionLabels(D)
}

## ======= EPIDEMIOLOGY ===========

makeAttributes <- function(D){
    nrep <- nrow(D)
    D[,id:=1:nrep]
    fx <- list(age=agelevels,tb=tblevels,hiv=hivlevels,art=artlevels, isoz = isoz)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    ## --- age
    D[,value:=ifelse(age=='5-14',1-dh.fracO,dh.fracO)] #NOTE value first set here
    ## --- HIV/ART
    D[,hivprev.u5:=HHhivprev04]
    D[,hivprev.o5:=HHhivprev514]
    D[,artcov:=1]
    D[,h01:=0]
    D[age!='5-14',h10:=hivprev.u5*(1-artcov)]
    D[age=='5-14',h10:=hivprev.o5*(1-artcov)]
    D[age!='5-14',h00:=1-hivprev.u5]
    D[age=='5-14',h00:=1-hivprev.o5]
    D[age!='5-14',h11:=hivprev.u5*artcov]
    D[age=='5-14',h11:=hivprev.o5*artcov]
    D[hiv==0 & art==0,value:=value*h00]
    D[hiv==0 & art==1,value:=value*h01]
    D[hiv==1 & art==0,value:=value*h10]
    D[hiv==1 & art==1,value:=value*h11]
    D[,c('h00','h01','h10','h11'):=NULL]
    ## --- TB
    ## ## (old version) calculate true TB prev among initial care seeking as:
    ## ## tbi = f x tbp + (1-f) x tbd
    ## ## where: f=fraction initially seeking care at PHC; tbp=prev @ phc; tbd=prev @ dh
    ## ## NOTE the 'underlying' TB prev in care-seekers in controlled by soc parms
    D[,tbi:= phc.presented * phc.tbprev + (1-phc.presented) * dh.tbprev]
    D[tb=='noTB',value:=value*(1-tbi)]
    D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    D[tb=='TB+',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]

    # D[tb=='noTB',value:=value*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    # D[tb=='TB',value:=value*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    # D[,tbi:=NULL]                            #remove temporary variable
    return(D)
}


## function for generating random sample of costs
MakeCostData <- function(csts,          #base data table of cost data
                         nrep,          #number of replicates being used in PSA
                         anmz=NULL     #attribute names (if any)
                         ){
  if(nrow(csts[cost.sd>0 & cost.m==0])>0) warning(paste0('Some cost input variables have zero mean & SD>0. These will be treated as fixed variables:\n',paste0(csts[cost.sd>0 & cost.m==0,cost],collapse='\n')))
  if(is.null(anmz)& any(csts[,table(cost)]>1)) warning('Some cost names occur >1 times, but no attributes have been specified! This is unlikely to do what you want.')
  csts[cost.m>0,gmsc:=cost.sd^2/cost.m]
  csts[!is.na(gmsc) & gmsc > 0, gmk:=cost.m/gmsc]
  NR <- nrow(csts)
  csts <- csts[rep(1:NR,nrep)]
  csts[,id:=rep(1:nrep,each=NR)]
  csts[,rnd:=!is.na(gmsc) & !is.na(gmk) & gmk>0 & gmsc > 0]
  csts[rnd==TRUE,value:=rgamma(sum(rnd),shape=gmk,scale = gmsc)] #random sample from gamma distribution
  csts[rnd!=TRUE,value:=cost.m]                                  #fixed values
  ## csts[,cnms:=paste0('c_',cost)]
  csts[,cnms:=paste0(cost)]
  F <- 'id '
  if(!is.null(anmz)) F <- paste0(F,'+ ',paste(anmz,collapse='+')) #split out by attributes if included
  F <- paste0(F, ' ~ cnms')
  dcast(csts,as.formula(F),value.var = 'value')      #id ~ cnms
}
## NOTE
## if attributes are included, all costs need to be specified by them even if this means duplicating those without dependence


## making life years
GetLifeYears <- function(isolist,discount.rate,yearfrom){
    ## template:
    LYT <- data.table(age=0:14,
                      age_group=c(rep('0-4',5),rep('5-14',10)),
                      LYS=0.0)
    ## make country/age key
    LYK <- list()
    for(iso in isolist){
        ## iso <- cn
        tmp <- copy(LYT)
        tmp[,iso3:=iso]
        for(ag in tmp$age)
            tmp[age==ag,LYS:=discly::discly(iso3=iso,
                                            age=ag,
                                            yearnow=yearfrom,
                                            sex='Total',
                                            endyear = 2098,
                                            HR=1,
                                            dr=discount.rate,
                                            hiv='both'
                                            )]
        LYK[[iso]] <- tmp
    }
    LYK <- rbindlist(LYK)
    ## assume unweighted & collapse
    LYK <- LYK[,.(LYS=mean(LYS)),by=.(iso3,age=age_group)]
    setkey(LYK,age)
    LYK
}

## ## scraps for development of below fn
## data <- merge(D,LYK,by='age') #add age
## Kmax <- 1e3
## file.id <- 'test'
## wtp <- 500

## NOTE this is more illustrative for now
## NOTE needs a folder called plots/ creating (which is currently excluded from the repo)
## some automatic CEA outputs
file.id='';Kmax=5e3;wtp=5e3;
MakeCEAoutputs <- function(data,LY,
                           file.id='',Kmax=5e3,wtp=5e3,
                           arms=c('SOC','INT')){
  data <- merge(data,LY,by='age') #add age
  DS <- D[,.(cost.SOC=sum(cost.soc*value),
                cost.INT=sum(cost.int*value),
                lyl.SOC=sum(deaths.soc*value*LYS),
                lyl.INT=sum(deaths.int*value*LYS)),
             by=id] #PSA summary

  ## prep for BCEA
  LYS <- CST <- matrix(nrow=nreps,ncol=2)
  LYS[,1] <- 1-DS$lyl.SOC #NOTE this is life years lost
  LYS[,2] <- 1-DS$lyl.INT
  CST[,1] <- DS$cost.SOC
  CST[,2] <- DS$cost.INT
  ## BCEA outputs
  M <- bcea(e=LYS,c=CST,ref=1,interventions = arms,Kmax=Kmax)
  print(summary(M))

  fn <- paste0(here('outdata/kstar_'),file.id,'.txt')
  cat(M$kstar,file = fn)
  fn <- paste0(here('outdata/ICER_'),file.id,'.txt')
  cat(M$ICER,file = fn)

  ## NOTE may need more configuration
  ceac.plot(M,graph='ggplot2') +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('plots/CEAC_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  ceplane.plot(M,graph='ggplot2',wtp=wtp)+
    scale_x_continuous(label=comma) +
    theme_classic() +
    theme(legend.position = 'top') + ggpubr::grids()
  fn <- paste0(here('plots/CE_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  eib.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('plots/EIB_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  evi.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('plots/EVI_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

}

## --- for reformatting costs
reformatCosts <- function(rcsts){
  iextra <- outer(isoz,c('.lo','.hi','drop'),paste0)
  iextra <- c(t(iextra)); iextra <- rev(rev(iextra)[-1])
  nnmz <- c('drop','DESCRIPTION','NAME',iextra)
  names(rcsts)[1:length(nnmz)] <- nnmz
  drop <- grep('drop',names(rcsts),value=TRUE)
  rcsts[,c(drop):=NULL]
  rcsts[is.na(rcsts)] <- 0.1 #dummy
  rcsts <- melt(rcsts,id=c('NAME','DESCRIPTION'))
  rcsts[,DESCRIPTION:=NULL]
  rcsts[,c('iso3','hilo'):=tstrsplit(variable,split="\\.")]
  rcsts <- dcast(rcsts,iso3 + NAME ~ hilo,value.var = 'value')
  rcsts[,c('cost.m','cost.sd'):=.((lo+hi)/2,(hi-lo)/3.92)]
  rcsts <- rcsts[,.(iso3,cost=NAME,cost.m,cost.sd)]
  rcsts
}

## calculate mid/lo/hi
MLH <- function(dat){
  nnmz <- names(dat)
  lnmz <- paste0(nnmz,'.lo')
  hnmz <- paste0(nnmz,'.hi')
  mnmz <- paste0(nnmz,'.mid')
  L <- dat[,lapply(.SD,lo),.SDcols=nnmz]
  M <- dat[,lapply(.SD,mean),.SDcols=nnmz]
  H <- dat[,lapply(.SD,hi),.SDcols=nnmz]
  setnames(L,nnmz,lnmz); setnames(M,nnmz,mnmz); setnames(H,nnmz,hnmz);
  list(L=L,M=M,H=H)
}
## MLH(out[,.(DcostperATT,DcostperATT.soc)]) #test


## =========== output formatters
outsummary <- function(out){

  keep <- c('costperATT.soc','costperATT.int',
            'DcostperATT',
            'Ddeaths',
            'DLYL',
            'DLYL0',
            'Dcost',
            'attPC',
            'DcostperLYS0',
            'DcostperLYS',
            'Dcostperdeaths',
            ## D/D
            'DcostperDATT',
            'DcostperDLYS0',
            'DcostperDLYS',
            'DcostperDdeaths')
  scr <- c(psoc.sc,pint.sc)
  scrm <- paste0(scr,'.mid')
  keep <- c(keep,scr)

  ## mid/lo/hi
  outa <- MLH(out[,..keep])

  ## more bespoke statistics
  outi <- out[,.(ICER= -mean(Dcost) / mean(DLYL))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  pouts <- outs[,.(costperATT.soc = brkt(costperATT.soc.mid,costperATT.soc.lo,costperATT.soc.hi),
                   costperATT.int = brkt(costperATT.int.mid,costperATT.int.lo,costperATT.int.hi),
                   DcostperATT = brkt(DcostperATT.mid,DcostperATT.lo,DcostperATT.hi),
                   DcostperLYS0 = brkt(DcostperLYS0.mid,
                                           DcostperLYS0.lo,DcostperLYS0.hi),
                   DcostperLYS = brkt(DcostperLYS.mid,DcostperLYS.lo,DcostperLYS.hi),
                   Dcostperdeaths = brkt(Dcostperdeaths.mid,
                                             Dcostperdeaths.lo,Dcostperdeaths.hi),
                   ## D/D
                   DcostperDATT = brkt(DcostperDATT.mid,DcostperDATT.lo,DcostperDATT.hi),
                   DcostperDLYS0 = brkt(DcostperDLYS0.mid,
                                            DcostperDLYS0.lo,DcostperDLYS0.hi),
                   DcostperDLYS = brkt(DcostperDLYS.mid,DcostperDLYS.lo,DcostperDLYS.hi),
                   DcostperDdeaths = brkt(DcostperDdeaths.mid,
                                              DcostperDdeaths.lo,DcostperDdeaths.hi),
                   ## end D/D
                   DcostPerOPD = brkt(Dcost.mid,Dcost.lo,Dcost.hi),
                   DdeathsPer100kOPD = brkt(-1e5*Ddeaths.mid,
                                                -1e5*Ddeaths.hi,-1e5*Ddeaths.lo),
                   DLYS0Per100kOPD = brkt(-1e5*DLYL0.mid,
                                              -1e5*DLYL0.hi,-1e5*DLYL0.lo),
                   DLYSPer100kOPD = brkt(-1e5*DLYL.mid,
                                             -1e5*DLYL.hi,-1e5*DLYL.lo),
                   attPC = brkt(attPC.mid,attPC.lo,attPC.hi),
                   DattPC = brkt(attPC.mid-1e2,attPC.lo-1e2,attPC.hi-1e2),
                   ICER=round(ICER,0))]

  ## staged costs
  scouts <- outa$M[,..scrm]

  ## return value
  list(outs=outs,pouts=pouts,scouts=scouts)
}


## ---- utilities for making CEACs
make.ceac <- function(CEA,lamz){
    crv <- lamz
    for(i in 1:length(crv)) crv[i] <- CEA[,mean(lamz[i]*Q-P>0)]
    crv
}


## additional table 2
## --- cols:
## CMR, UGA x SOC, INT
## --- rows:
## contacts = value //
## TPT courses = tpt //
## ATT courses = att //
## prev TB = prevtb
## inc TB = inctb
## prev deaths = deaths-incdeaths
## inc tb deaths = incdeaths
## discounted LYL TODO //
## ATT cost TODO
## TPT cost TODO
## total cost TODO //
## ICER TODO



## =========== output formatters
Table2 <- function(dat){

  ## mid/lo/hi
  outa <- MLH(dat[,.(
    Dpttb,pttb.int,pttb.soc,
    DPHC.presumptive,PHC.presumptive.int,PHC.presumptive.soc,
    DDH.presumptive,DH.presumptive.int,DH.presumptive.soc,
    Dpresumptive,presumptive.int,presumptive.soc,
    DPHC.evaluated,PHC.evaluated.int,PHC.evaluated.soc,
    DDH.evaluated,DH.evaluated.int,DH.evaluated.soc,
    Dassessments,assessments.int,assessments.soc,
    Dpassessphc,passessphc.int,passessphc.soc,
    Dbacassess,bacassess.int,bacassess.soc,
    Dpbacassessphc,pbacassessphc.int,pbacassessphc.soc,
    Drefers,refers.int,refers.soc,
    Ddx,dx.int,dx.soc,
    Dpdxphc,pdxphc.int,pdxphc.soc,
    Dpdxb,pdxb.int,pdxb.soc,
    DPHC.treated,PHC.treated.int,PHC.treated.soc,
    DDH.treated,DH.treated.int,DH.treated.soc,
    Datt,att.int,att.soc,
    Dpatt.phc,patt.phc.int,patt.phc.soc,
    Dpatt.bac,patt.bac.int,patt.bac.soc,
    Dpttbtx, pttbtx.int, pttbtx.soc,
    Dpftbtx, pftbtx.int, pftbtx.soc,
    Ddeaths,deaths.int,deaths.soc,
    DLYL0,LYL0.int,LYL0.soc,
    DLYL,LYL.int,LYL.soc,
    DPHC.evaluated.cost,PHC.evaluated.cost.int,PHC.evaluated.cost.soc,
    DDH.evaluated.cost,DH.evaluated.cost.int,DH.evaluated.cost.soc,
    Dcost.assessments,cost.assessments.int,cost.assessments.soc,
    DPHC.treated.cost,PHC.treated.cost.int,PHC.treated.cost.soc,
    DDH.treated.cost,DH.treated.cost.int,DH.treated.cost.soc,
    Dcost,cost.int,cost.soc
  )])

  ## more bespoke statistics
  outi <- dat[,.(ICER= -mean(Dcost) / mean(DLYL),
                 DICER= -mean(Dcost) / mean(DLYL0))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  fac <- 1e2 #per  per 100 children with presumptive TB (actually presenting for now)
  pouts <- outs[,.(
    Dpttb = brkt(fac*Dpttb.mid,fac*Dpttb.lo,fac*Dpttb.hi),
    pttb.int = brkt(fac*pttb.int.mid,fac*pttb.int.lo,fac*pttb.int.hi),
    pttb.soc = brkt(fac*pttb.soc.mid,fac*pttb.soc.lo,fac*pttb.soc.hi),
    DPHC.presumptive = brkt(fac*DPHC.presumptive.mid,fac*DPHC.presumptive.lo,fac*DPHC.presumptive.hi),
    PHC.presumptive.int = brkt(fac*PHC.presumptive.int.mid,fac*PHC.presumptive.int.lo,fac*PHC.presumptive.int.hi),
    PHC.presumptive.soc = brkt(fac*PHC.presumptive.soc.mid,fac*PHC.presumptive.soc.lo,fac*PHC.presumptive.soc.hi),
    DDH.presumptive = brkt(fac*DDH.presumptive.mid,fac*DDH.presumptive.lo,fac*DDH.presumptive.hi),
    DH.presumptive.int = brkt(fac*DH.presumptive.int.mid,fac*DH.presumptive.int.lo,fac*DH.presumptive.int.hi),
    DH.presumptive.soc = brkt(fac*DH.presumptive.soc.mid,fac*DH.presumptive.soc.lo,fac*DH.presumptive.soc.hi),
    Dpresumptive = brkt(fac*Dpresumptive.mid,fac*Dpresumptive.lo,fac*Dpresumptive.hi),
    presumptive.int = brkt(fac*presumptive.int.mid,fac*presumptive.int.lo,fac*presumptive.int.hi),
    presumptive.soc = brkt(fac*presumptive.soc.mid,fac*presumptive.soc.lo,fac*presumptive.soc.hi),
    DPHC.evaluated = brkt(fac*DPHC.evaluated.mid,fac*DPHC.evaluated.lo,fac*DPHC.evaluated.hi),
    PHC.evaluated.int = brkt(fac*PHC.evaluated.int.mid,fac*PHC.evaluated.int.lo,fac*PHC.evaluated.int.hi),
    PHC.evaluated.soc = brkt(fac*PHC.evaluated.soc.mid,fac*PHC.evaluated.soc.lo,fac*PHC.evaluated.soc.hi),
    DDH.evaluated = brkt(fac*DDH.evaluated.mid,fac*DDH.evaluated.lo,fac*DDH.evaluated.hi),
    DH.evaluated.int = brkt(fac*DH.evaluated.int.mid,fac*DH.evaluated.int.lo,fac*DH.evaluated.int.hi),
    DH.evaluated.soc = brkt(fac*DH.evaluated.soc.mid,fac*DH.evaluated.soc.lo,fac*DH.evaluated.soc.hi),
    Dassessments = brkt(fac*Dassessments.mid,fac*Dassessments.lo,fac*Dassessments.hi),
    assessments.int = brkt(fac*assessments.int.mid,fac*assessments.int.lo,fac*assessments.int.hi),
    assessments.soc = brkt(fac*assessments.soc.mid,fac*assessments.soc.lo,fac*assessments.soc.hi),
    Dpassessphc = brkt(fac*Dpassessphc.mid,fac*Dpassessphc.lo,fac*Dpassessphc.hi),
    passessphc.int = brkt(fac*passessphc.int.mid,fac*passessphc.int.lo,fac*passessphc.int.hi),
    passessphc.soc = brkt(fac*passessphc.soc.mid,fac*passessphc.soc.lo,fac*passessphc.soc.hi),
    Dbacassess = brkt(fac*Dbacassess.mid,fac*Dbacassess.lo,fac*Dbacassess.hi),
    bacassess.int = brkt(fac*bacassess.int.mid,fac*bacassess.int.lo,fac*bacassess.int.hi),
    bacassess.soc = brkt(fac*bacassess.soc.mid,fac*bacassess.soc.lo,fac*bacassess.soc.hi),
    Dpbacassessphc = brkt(fac*Dpbacassessphc.mid,fac*Dpbacassessphc.lo,fac*Dpbacassessphc.hi),
    pbacassessphc.int = brkt(fac*pbacassessphc.int.mid,fac*pbacassessphc.int.lo,fac*pbacassessphc.int.hi),
    pbacassessphc.soc = brkt(fac*pbacassessphc.soc.mid,fac*pbacassessphc.soc.lo,fac*pbacassessphc.soc.hi),
    Drefers = brkt(fac*Drefers.mid,fac*Drefers.lo,fac*Drefers.hi),
    refers.int = brkt(fac*refers.int.mid,fac*refers.int.lo,fac*refers.int.hi),
    refers.soc = brkt(fac*refers.soc.mid,fac*refers.soc.lo,fac*refers.soc.hi),
    Ddx = brkt(fac*Ddx.mid,fac*Ddx.lo,fac*Ddx.hi),
    dx.int = brkt(fac*dx.int.mid,fac*dx.int.lo,fac*dx.int.hi),
    dx.soc = brkt(fac*dx.soc.mid,fac*dx.soc.lo,fac*dx.soc.hi),
    Dpdxphc = brkt(fac*Dpdxphc.mid,fac*Dpdxphc.lo,fac*Dpdxphc.hi),
    pdxphc.int = brkt(fac*pdxphc.int.mid,fac*pdxphc.int.lo,fac*pdxphc.int.hi),
    pdxphc.soc = brkt(fac*pdxphc.soc.mid,fac*pdxphc.soc.lo,fac*pdxphc.soc.hi),
    Dpdxb = brkt(fac*Dpdxb.mid,fac*Dpdxb.lo,fac*Dpdxb.hi),
    pdxb.int = brkt(fac*pdxb.int.mid,fac*pdxb.int.lo,fac*pdxb.int.hi),
    pdxb.soc = brkt(fac*pdxb.soc.mid,fac*pdxb.soc.lo,fac*pdxb.soc.hi),
    DPHC.treated = brkt(fac*DPHC.treated.mid,fac*DPHC.treated.lo,fac*DPHC.treated.hi),
    PHC.treated.int = brkt(fac*PHC.treated.int.mid,fac*PHC.treated.int.lo,fac*PHC.treated.int.hi),
    PHC.treated.soc = brkt(fac*PHC.treated.soc.mid,fac*PHC.treated.soc.lo,fac*PHC.treated.soc.hi),
    DDH.treated = brkt(fac*DDH.treated.mid,fac*DDH.treated.lo,fac*DDH.treated.hi),
    DH.treated.int = brkt(fac*DH.treated.int.mid,fac*DH.treated.int.lo,fac*DH.treated.int.hi),
    DH.treated.soc = brkt(fac*DH.treated.soc.mid,fac*DH.treated.soc.lo,fac*DH.treated.soc.hi),
    Datt = brkt(fac*Datt.mid,fac*Datt.lo,fac*Datt.hi),
    att.int = brkt(fac*att.int.mid,fac*att.int.lo,fac*att.int.hi),
    att.soc = brkt(fac*att.soc.mid,fac*att.soc.lo,fac*att.soc.hi),
    Dpatt.phc = brkt(fac*Dpatt.phc.mid,fac*Dpatt.phc.lo,fac*Dpatt.phc.hi),
    patt.phc.int = brkt(fac*patt.phc.int.mid,fac*patt.phc.int.lo,fac*patt.phc.int.hi),
    patt.phc.soc = brkt(fac*patt.phc.soc.mid,fac*patt.phc.soc.lo,fac*patt.phc.soc.hi),
    Dpatt.bac = brkt(fac*Dpatt.bac.mid,fac*Dpatt.bac.lo,fac*Dpatt.bac.hi, ndp=0),
    patt.bac.int = brkt(fac*patt.bac.int.mid,fac*patt.bac.int.lo,fac*patt.bac.int.hi, ndp=0),
    patt.bac.soc = brkt(fac*patt.bac.soc.mid,fac*patt.bac.soc.lo,fac*patt.bac.soc.hi, ndp=0),
    Dpttbtx = brkt(fac*Dpttbtx.mid,fac*Dpttbtx.lo,fac*Dpttbtx.hi),
    pttbtx.int = brkt(fac*pttbtx.int.mid,fac*pttbtx.int.lo,fac*pttbtx.int.hi),
    pttbtx.soc = brkt(fac*pttbtx.soc.mid,fac*pttbtx.soc.lo,fac*pttbtx.soc.hi),
    # Dpttbtx = brkt(fac*Dpttbtx.mid,fac*Dpttbtx.lo,fac*Dpttbtx.hi),
    # pttbtx.int = brkt(fac*pttbtx.int.mid,fac*pttbtx.int.lo,fac*pttbtx.int.hi),
    # pttbtx.soc = brkt(fac*pttbtx.soc.mid,fac*pttbtx.soc.lo,fac*pttbtx.soc.hi),
    Dpftbtx = brkt(fac*Dpftbtx.mid,fac*Dpftbtx.lo,fac*Dpftbtx.hi),
    pftbtx.int = brkt(fac*pftbtx.int.mid,fac*pftbtx.int.lo,fac*pftbtx.int.hi),
    pftbtx.soc = brkt(fac*pftbtx.soc.mid,fac*pftbtx.soc.lo,fac*pftbtx.soc.hi),
    Ddeaths = brkt(fac*Ddeaths.mid,fac*Ddeaths.lo,fac*Ddeaths.hi),
    deaths.int = brkt(fac*deaths.int.mid,fac*deaths.int.lo,fac*deaths.int.hi),
    deaths.soc = brkt(fac*deaths.soc.mid,fac*deaths.soc.lo,fac*deaths.soc.hi),
    DLYL = brkt(fac*DLYL.mid,fac*DLYL.lo,fac*DLYL.hi),
    LYL.int = brkt(fac*LYL.int.mid,fac*LYL.int.lo,fac*LYL.int.hi),
    LYL.soc = brkt(fac*LYL.soc.mid,fac*LYL.soc.lo,fac*LYL.soc.hi),
    DLYL0 = brkt(fac*DLYL0.mid,fac*DLYL0.lo,fac*DLYL0.hi),
    LYL0.int = brkt(fac*LYL0.int.mid,fac*LYL0.int.lo,fac*LYL0.int.hi),
    LYL0.soc = brkt(fac*LYL0.soc.mid,fac*LYL0.soc.lo,fac*LYL0.soc.hi),
    DPHC.evaluated.cost = brkt(fac*DPHC.evaluated.cost.mid,fac*DPHC.evaluated.cost.lo,fac*DPHC.evaluated.cost.hi),
    PHC.evaluated.cost.int = brkt(fac*PHC.evaluated.cost.int.mid,fac*PHC.evaluated.cost.int.lo,fac*PHC.evaluated.cost.int.hi),
    PHC.evaluated.cost.soc = brkt(fac*PHC.evaluated.cost.soc.mid,fac*PHC.evaluated.cost.soc.lo,fac*PHC.evaluated.cost.soc.hi),
    DDH.evaluated.cost = brkt(fac*DDH.evaluated.cost.mid,fac*DDH.evaluated.cost.lo,fac*DDH.evaluated.cost.hi),
    DH.evaluated.cost.int = brkt(fac*DH.evaluated.cost.int.mid,fac*DH.evaluated.cost.int.lo,fac*DH.evaluated.cost.int.hi),
    DH.evaluated.cost.soc = brkt(fac*DH.evaluated.cost.soc.mid,fac*DH.evaluated.cost.soc.lo,fac*DH.evaluated.cost.soc.hi),
    Dcost.assessments = brkt(fac*Dcost.assessments.mid,fac*Dcost.assessments.lo,fac*Dcost.assessments.hi),
    cost.assessments.int = brkt(fac*cost.assessments.int.mid,fac*cost.assessments.int.lo,fac*cost.assessments.int.hi),
    cost.assessments.soc = brkt(fac*cost.assessments.soc.mid,fac*cost.assessments.soc.lo,fac*cost.assessments.soc.hi),
    DPHC.treated.cost = brkt(fac*DPHC.treated.cost.mid,fac*DPHC.treated.cost.lo,fac*DPHC.treated.cost.hi),
    PHC.treated.cost.int = brkt(fac*PHC.treated.cost.int.mid,fac*PHC.treated.cost.int.lo,fac*PHC.treated.cost.int.hi),
    PHC.treated.cost.soc = brkt(fac*PHC.treated.cost.soc.mid,fac*PHC.treated.cost.soc.lo,fac*PHC.treated.cost.soc.hi),
    DDH.treated.cost = brkt(fac*DDH.treated.cost.mid,fac*DDH.treated.cost.lo,fac*DDH.treated.cost.hi),
    DH.treated.cost.int = brkt(fac*DH.treated.cost.int.mid,fac*DH.treated.cost.int.lo,fac*DH.treated.cost.int.hi),
    DH.treated.cost.soc = brkt(fac*DH.treated.cost.soc.mid,fac*DH.treated.cost.soc.lo,fac*DH.treated.cost.soc.hi),
    Dcost = brkt(fac*Dcost.mid,fac*Dcost.lo,fac*Dcost.hi),
    cost.int = brkt(fac*cost.int.mid,fac*cost.int.lo,fac*cost.int.hi),
    cost.soc = brkt(fac*cost.soc.mid,fac*cost.soc.lo,fac*cost.soc.hi),
    ICER=round(ICER,0),
    DICER=round(DICER,0)
  )]

## return value
  list(outs=outs,pouts=pouts)
}
