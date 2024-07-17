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

## == LTBI infection probability
#NB this is LTBI given not active: it is taken to be max(0,LTBI-coprev)
ltbi.prev <- function(age,coprev,hinco=FALSE){
        if(length(age)>1 & length(hinco)==1) hinco <- rep(hinco,length(age))
        tmp <- P$LTBI04$r(length(age))
        tmp[hinco] <- P$LTBI04hi$r(sum(hinco))
        tmp[age>=5] <- P$LTBI514$r(sum(age>=5))
        tmp[age>=5 & hinco] <- P$LTBI514hi$r(sum(age>=5 & hinco))
        tmp
        ## pmax(0,tmp - coprev) # already taken into account with decision tree
}
# ltbi.prev(1:10,0.1,hinco=TRUE)

## progression probability
progprob <- function(age,hiv=0,art=0){
        if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
        if(length(age)>1 & length(art)==1) art <- rep(art,length(age))
        ans <- P$prog04$r(length(age))
        ans[age>=5] <- P$prog514$r(sum(age>=5))
        if(any(hiv>0)){ #treat as IRR for escape
                hr <- P$hivpi$r(sum(hiv>0))
                ans[hiv>0] <- 1-(1-ans[hiv>0])^hr
        }
        if(any(art>0)){ #treat as IRR for escape
                hr <- P$artp$r(sum(art>0))
                ans[art>0] <- 1-(1-ans[art>0])^hr
        }
        ans
}

## add CFRs to data by side-effect
AddProgProb <- function(D, P){
  D[,p.tbdx.1yr:=progprob(age,hiv,art) * ltbi.prev(age,0)] #NOTE handling of coprev happens explicitly
}

## progprob(c(rep(3,5),rep(10,5)))


## === IPT efficacy
TPTrr <- function(age,hiv=0,
                  tst='none'           #a flag for: given to TST+ or not
){
        if(length(age)>1 & length(hiv)==1) hiv <- rep(hiv,length(age))
        if(tst=='none'){
                ans <- P$tptRR$r(length(age))
                ans[hiv>0] <- P$tptRRhivpos$r(sum(hiv>0))  #HIV+
        } else {
                ans <- P$tptRRhivpos$r(length(age))  #NOTE only applied to HIV-ves
        }
        ans
}

# TPTrr(1:10, P=P)
# summary(TPTrr(runif(1e3),hiv=0)) #0.37
# summary(TPTrr(runif(1e3),hiv=1)) #0.35
# summary(TPTrr(runif(1e3),tst='yes')) #0.09
AddTPTrr <- function(D,P){
        D[,tptRR0:=TPTrr(age,hiv)]  #base efficacy of IPT
        D[,tptRR1:=TPTrr(age,tst="+ve")]   #base efficacy of PT: TST+ve
        D[,tptRR:=tptRR0] 
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

## D[,tmp:=0+(int.frac.bac.assess)*(0+(int.frac.bac.dx)*(1+(1)*(0+(1)*(0))+(1-int.frac.bac.dx)*(0+(1)*(0+(int.frac.bac.clin.dx)*(1+(1)*(0))+(int.frac.bac.noclin.dx*(1-int.frac.bac.clin.dx))*(0+(1)*(0))+((1-int.frac.bac.clin.dx)*(1-int.frac.bac.noclin.dx))*(0+(int.frac.bac.7d.clin.dx)*(1+(1)*(0))+(int.frac.bac.7d.noclin.dx*(1-int.frac.bac.7d.clin.dx))*(0+(1)*(0))+((1-int.frac.bac.7d.clin.dx)*(1-int.frac.bac.7d.noclin.dx))*(0)))))+(1-int.frac.bac.assess)*(0+(int.frac.clin.dx)*(1+(1)*(0))+(int.frac.noclin.dx*(1-int.frac.clin.dx))*(0+(1)*(0))+((1-int.frac.clin.dx)*(1-int.frac.noclin.dx))*(0+(int.frac.clin.7d.clin.dx)*(1+(1)*(0))+(int.frac.clin.7d.noclin.dx*(1-int.frac.clin.7d.clin.dx))*(0+(1)*(0))+((1-int.frac.clin.7d.clin.dx)*(1-int.frac.clin.7d.noclin.dx))*(0))))+(1-int.frac.symp.attending)*(0)]

## ========= DIAGNOSIS ===============

## function for combining sample modality with bacteriological test
## to calculate the probablility bac+ TB is diagnosed
TBbacsampletest <- function(samplepossible,testpos){
  samplepossible * testpos
}


## OR function - with NPA/stool in mind
## positive if either sample/test combo is
TBbac4ST1orST2 <- function(S1,T1,S2,T2){
  1 - (1-TBbacsampletest(S1,T1)) * (1-TBbacsampletest(S2,T2))
}

## function to add in the Sample/Test combined probabilities of TB dx
## works by side effect
AddSampleTests <- function(D){
  
  ## Bacteriological test availability
  # use Odume et al., 2023 to approximate this
  # Plot in Odume et al., 2023 suggests 
  # 0-10% availability of Xpert Ultra in PHC
  # ~20-100% availability of Xpert Ultra in DH
  
  # DH
  # (p <- getAB(0.60, ((20-100)^2)/392^2))
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  D[,soc.dh.test:=rbeta(nrow(D),2.85744,1.90496)]
  D[,int.dh.test:=rbeta(nrow(D),2.85744,1.90496)] # to increase this with Truenat intorduction
  
  # PHC
  # (p <- getAB(0.05, ((0-10)^2)/392^2))
  # curve(dbeta(x, p$a, p$b), from=0, to=1, n=200)
  # summary(rbeta(1000, p$a, p$b))
  D[,soc.phc.test:=rbeta(nrow(D),3.59952,68.39088)]
  D[,int.phc.test:=rbeta(nrow(D),3.59952,68.39088)]
  
  ## ------- X on sputum/GA -------
  ## People receiving Xpert Ultra testing [either sputum or GA], in those identified as having presumptive TB
  D[,soc.dh.fracUltra:=ifelse(age=='5-14',dh.fracUltra,dh.fracUltra)]
  D[,int.dh.fracUltra:=ifelse(age=='5-14',dh.fracUltra,dh.fracUltra)]
  
  D[,soc.phc.fracUltra:=ifelse(age=='5-14',phc.fracUltra,phc.fracUltra)]
  D[,int.phc.fracUltra:=ifelse(age=='5-14',phc.fracUltra,phc.fracUltra)]
  
  # TB dx bac+ on Xpert Ultra on sputum
  D[,soc.dh.bact.tbdx:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  D[,int.dh.bact.tbdx:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  
  D[,soc.phc.bact.tbdx:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  D[,int.phc.bact.tbdx:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  
  D[,soc.dh.ptbxUltra:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  D[,int.dh.ptbxUltra:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  
  D[,soc.phc.ptbxUltra:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  D[,int.phc.ptbxUltra:=ifelse(tb=='TB',sens.xsputum,1-spec.xsputum)]
  
  # TB dx bac+ on TrueNat on stool
  D[,soc.dh.ptbxTrueNat:=ifelse(tb=='TB',sens.xstool,1-spec.xstool)]
  D[,int.dh.ptbxTrueNat:=ifelse(tb=='TB',sens.xstool,1-spec.xstool)]
  
  D[,soc.phc.ptbxTrueNat:=ifelse(tb=='TB',sens.xstool,1-spec.xstool)]
  D[,int.phc.ptbxTrueNat:=ifelse(tb=='TB',sens.xstool,1-spec.xstool)]
  
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
  
}

## ======= COMBINED LABELLER ===========


## additional labels from data (may overwrite some initial version currently)
AddDataDrivenLabels <- function(D){

  # generate some parameters
  
  # pre-treatment loss to follow-up
  D[,dh.ptltfu:=1-dh.att]
  D[,phc.ptltfu:=1-phc.att]
  
  # pre-treatment loss to follow-up
  D[,int.phc.rltfu:=soc.phc.rltfu]
  
  # RIF resistance
  D[,prr:=dh.prr]
  
  # Bacteriological confirmation
  # D[,soc.dh.bact.tbdx:=1-dh.att]
  # D[,int.dh.bact.tbdx:=1-phc.att]
  
  # 
  D[,hivprev.u5:=0]
  D[,hivprev.o5:=0]
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
    fx <- list(age=agelevels,tb=tblevels,hiv=hivlevels,art=artlevels)
    cofx <- expand.grid(fx)
    cat('Attribute combinations used:\n')
    print(cofx)
    D <- D[rep(1:nrow(D),each=nrow(cofx))] #expand out data
    D[,names(cofx):=cofx[rep(1:nrow(cofx),nrep),]]
    
    ## --- age
    # D[,F.u5:=ifelse(isoz=='CMR', cmr.frac.u5, uga.frac.u5)]
    # D[,artcov:=ifelse(isoz=='CMR', cmr.artcov, uga.artcov)]
    # D[,artcov:=0.99]
 
    # D[,hivprev.u5:=ifelse(isoz=='CMR', cmr.frac.hiv, uga.frac.hiv)]
    # D[,hivprevO:=hivprev.u5]

    D[,value:=ifelse(age=='5-14',1-dh.fracO,dh.fracO)] #NOTE value first set here

    ## --- HIV/ART
    D[,hivprev.u5:=0.01]
    D[,hivprev.o5:=0.01]
    # D[,artcov:=0]
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
    ## D[,tbi:= d.soc.pphc * d.phc.tbinprsmptv + (1-d.soc.pphc) * d.dh.tbinprsmptv]
    ## D[tb=='noTB',value:=value*(1-tbi)]
    ## D[tb=='TB-',value:=value*tbi*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    ## D[tb=='TB',value:=value*tbi*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
    ## D[,tbi:=NULL]                            #remove temporary variable
    # ## (new version) based on data:
    # D[,tbi:=ifelse(age=='5-14',d.TBprev.ICS.o5,d.TBprev.ICS.u5)]
    # D[tb=='noTB',value:=value*(1-tbi)]
    D[tb=='noTB',value:=value*ifelse(age=='5-14',1-Fbc.o5,1-Fbc.u5)] #NOTE assuming no TB outside of presumptive?
    D[tb=='TB',value:=value*ifelse(age=='5-14',Fbc.o5,Fbc.u5)]
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
## NOTE needs a folder called graphs/ creating (which is currently excluded from the repo)
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
  fn <- paste0(here('graphs/CEAC_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  ceplane.plot(M,graph='ggplot2',wtp=wtp)+
    scale_x_continuous(label=comma) +
    theme_classic() +
    theme(legend.position = 'top') + ggpubr::grids()
  fn <- paste0(here('graphs/CE_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  eib.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EIB_'),file.id,'.png')
  ggsave(file=fn,w=7,h=7)

  evi.plot(M,graph='ggplot2',wtp=wtp) +
    scale_x_continuous(label=comma) +
    theme_classic() + ggpubr::grids()
  fn <- paste0(here('graphs/EVI_'),file.id,'.png')
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
## MLH(out[,.(DcostperATT.int,DcostperATT.soc)]) #test


## =========== output formatters
outsummary <- function(out){

  ## mid/lo/hi
  outa <- MLH(out[,.(costperATT.soc,costperATT.int,
                     DcostperATT.int,
                     costperTPT.soc,costperTPT.int,
                     DcostperTPT.int,
                     Ddeaths.int,
                     DLYL.int,
                     DLYL0.int,
                     DcostperLYS0.int,
                     DcostperLYS.int,
                     Dcostperdeaths.int,
                     Dcost.int)])

  ## more bespoke statistics
  outi <- out[,.(ICER.int= -mean(Dcost.int) / mean(DLYL.int))]

  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine

  ## pretty version
  pouts <- outs[,.(costperATT.soc = brkt(costperATT.soc.mid,costperATT.soc.lo,costperATT.soc.hi),
                   costperATT.int = brkt(costperATT.int.mid,costperATT.int.lo,costperATT.int.hi),
                   DcostperATT.int = brkt(DcostperATT.int.mid,DcostperATT.int.lo,DcostperATT.int.hi), # TODO::quick gap measure check!!
                   costperTPT.soc = brkt(costperTPT.soc.mid,costperTPT.soc.lo,costperTPT.soc.hi),
                   costperTPT.int = brkt(costperTPT.int.mid,costperTPT.int.lo,costperTPT.int.hi),
                   DcostperTPT.int = brkt(DcostperTPT.int.mid,DcostperTPT.int.lo,DcostperTPT.int.hi),
                   DcostperLYS0.int = brkt(DcostperLYS0.int.mid,
                                           DcostperLYS0.int.lo,DcostperLYS0.int.hi),
                  DcostperLYS.int = brkt(DcostperLYS.int.mid,DcostperLYS.int.lo,DcostperLYS.int.hi),
                  Dcostperdeaths.int = brkt(Dcostperdeaths.int.mid,
                                            Dcostperdeaths.int.lo,Dcostperdeaths.int.hi),
                  DcostPerChildContact.int = brkt(Dcost.int.mid,Dcost.int.lo,Dcost.int.hi),
                  DdeathsPer100kChildContacts.int = brkt(-1e5*Ddeaths.int.mid,
                                               -1e5*Ddeaths.int.hi,-1e5*Ddeaths.int.lo),
                  DLYS0Per100kChildContacts.int = brkt(-1e5*DLYL0.int.mid,
                                               -1e5*DLYL0.int.hi,-1e5*DLYL0.int.lo),
                  DLYSPer100kChildContacts.int = brkt(-1e5*DLYL.int.mid,
                                            -1e5*DLYL.int.hi,-1e5*DLYL.int.lo),
                  ICER.int=round(ICER.int,0))]

  ## return value
  list(outs=outs,pouts=pouts)
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
    Dcontacts,contacts.int,contacts.soc,
    Dtpt,tpt.int,tpt.soc,
    Datt,att.int,att.soc,
    DLYL,LYL.int,LYL.soc,
    Dprevtb,prevtb.int,prevtb.soc,
    Dinctb,inctb.int,inctb.soc,
    Dincdeaths,incdeaths.int,incdeaths.soc,
    Dprevdeaths,prevdeaths.int,prevdeaths.soc,
    Ddeaths,deaths.int,deaths.soc,
    Dcost.screen,cost.screen.int,cost.screen.soc,
    Dcost.tpt,cost.tpt.int,cost.tpt.soc,
    Dcost.prev.att,cost.prev.att.int,cost.prev.att.soc,
    Dcost.inc.att,cost.inc.att.int,cost.inc.att.soc,
    Dcost,cost.int,cost.soc
  )])
  
  ## more bespoke statistics
  outi <- dat[,.(ICER.int= -mean(Dcost) / mean(DLYL))]
  
  ## join
  outs <- do.call(cbind,list(outa$M,outa$L,outa$H,outi)) #combine
  
  ## pretty version
  fac <- 1e3 #per fac index cases
  pouts <- outs[,.(
    Dcontacts = brkt(fac*Dcontacts.mid,fac*Dcontacts.lo,fac*Dcontacts.hi),
    contacts.int = brkt(fac*contacts.int.mid,fac*contacts.int.lo,fac*contacts.int.hi),
    contacts.soc = brkt(fac*contacts.soc.mid,fac*contacts.soc.lo,fac*contacts.soc.hi),
    Dtpt = brkt(fac*Dtpt.mid,fac*Dtpt.lo,fac*Dtpt.hi),
    tpt.int = brkt(fac*tpt.int.mid,fac*tpt.int.lo,fac*tpt.int.hi),
    tpt.soc = brkt(fac*tpt.soc.mid,fac*tpt.soc.lo,fac*tpt.soc.hi),
    Datt = brkt(fac*Datt.mid,fac*Datt.lo,fac*Datt.hi),
    att.int = brkt(fac*att.int.mid,fac*att.int.lo,fac*att.int.hi),
    att.soc = brkt(fac*att.soc.mid,fac*att.soc.lo,fac*att.soc.hi),
    DLYL = brkt(fac*DLYL.mid,fac*DLYL.lo,fac*DLYL.hi),
    LYL.int = brkt(fac*LYL.int.mid,fac*LYL.int.lo,fac*LYL.int.hi),
    LYL.soc = brkt(fac*LYL.soc.mid,fac*LYL.soc.lo,fac*LYL.soc.hi),
    Dprevtb = brkt(fac*Dprevtb.mid,fac*Dprevtb.lo,fac*Dprevtb.hi),
    prevtb.int = brkt(fac*prevtb.int.mid,fac*prevtb.int.lo,fac*prevtb.int.hi),
    prevtb.soc = brkt(fac*prevtb.soc.mid,fac*prevtb.soc.lo,fac*prevtb.soc.hi),
    Dinctb = brkt(fac*Dinctb.mid,fac*Dinctb.lo,fac*Dinctb.hi),
    inctb.int = brkt(fac*inctb.int.mid,fac*inctb.int.lo,fac*inctb.int.hi),
    inctb.soc = brkt(fac*inctb.soc.mid,fac*inctb.soc.lo,fac*inctb.soc.hi),
    Dincdeaths = brkt(fac*Dincdeaths.mid,fac*Dincdeaths.lo,fac*Dincdeaths.hi),
    incdeaths.int = brkt(fac*incdeaths.int.mid,fac*incdeaths.int.lo,fac*incdeaths.int.hi),
    incdeaths.soc = brkt(fac*incdeaths.soc.mid,fac*incdeaths.soc.lo,fac*incdeaths.soc.hi),
    Dprevdeaths = brkt(fac*Dprevdeaths.mid,fac*Dprevdeaths.lo,fac*Dprevdeaths.hi),
    prevdeaths.int = brkt(fac*prevdeaths.int.mid,fac*prevdeaths.int.lo,fac*prevdeaths.int.hi),
    prevdeaths.soc = brkt(fac*prevdeaths.soc.mid,fac*prevdeaths.soc.lo,fac*prevdeaths.soc.hi),
    Ddeaths = brkt(fac*Ddeaths.mid,fac*Ddeaths.lo,fac*Ddeaths.hi),
    deaths.int = brkt(fac*deaths.int.mid,fac*deaths.int.lo,fac*deaths.int.hi),
    deaths.soc = brkt(fac*deaths.soc.mid,fac*deaths.soc.lo,fac*deaths.soc.hi),
    Dcost.screen = brkt(fac*Dcost.screen.mid,fac*Dcost.screen.lo,fac*Dcost.screen.hi),
    cost.screen.int = brkt(fac*cost.screen.int.mid,fac*cost.screen.int.lo,fac*cost.screen.int.hi),
    cost.screen.soc = brkt(fac*cost.screen.soc.mid,fac*cost.screen.soc.lo,fac*cost.screen.soc.hi),
    Dcost.tpt = brkt(fac*Dcost.tpt.mid,fac*Dcost.tpt.lo,fac*Dcost.tpt.hi),
    cost.tpt.int = brkt(fac*cost.tpt.int.mid,fac*cost.tpt.int.lo,fac*cost.tpt.int.hi),
    cost.tpt.soc = brkt(fac*cost.tpt.soc.mid,fac*cost.tpt.soc.lo,fac*cost.tpt.soc.hi),
    Dcost.prev.att = brkt(fac*Dcost.prev.att.mid,fac*Dcost.prev.att.lo,fac*Dcost.prev.att.hi),
    cost.prev.att.int = brkt(fac*cost.prev.att.int.mid,fac*cost.prev.att.int.lo,fac*cost.prev.att.int.hi),
    cost.prev.att.soc = brkt(fac*cost.prev.att.soc.mid,fac*cost.prev.att.soc.lo,fac*cost.prev.att.soc.hi),
    Dcost.inc.att = brkt(fac*Dcost.inc.att.mid,fac*Dcost.inc.att.lo,fac*Dcost.inc.att.hi),
    cost.inc.att.int = brkt(fac*cost.inc.att.int.mid,fac*cost.inc.att.int.lo,fac*cost.inc.att.int.hi),
    cost.inc.att.soc = brkt(fac*cost.inc.att.soc.mid,fac*cost.inc.att.soc.lo,fac*cost.inc.att.soc.hi),
    Dcost = brkt(fac*Dcost.mid,fac*Dcost.lo,fac*Dcost.hi),
    cost.int = brkt(fac*cost.int.mid,fac*cost.int.lo,fac*cost.int.hi),
    cost.soc = brkt(fac*cost.soc.mid,fac*cost.soc.lo,fac*cost.soc.hi),
    ICER.int=round(ICER.int,0)
  )]
  
  ## return value
  list(outs=outs,pouts=pouts)
}
