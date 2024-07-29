toget2 <- c('id',
'PHC.presumptive.soc','PHC.presumptive.int',
'PHC.evaluated.soc','PHC.evaluated.int',
'PHC.bacassess.soc', 'PHC.bacassess.int',
'PHC.diagnosed.soc', 'PHC.diagnosed.int',
'DH.presumptive.soc','DH.presumptive.int',
'DH.evaluated.soc','DH.evaluated.int',
'DH.bacassess.soc', 'DH.bacassess.int',
'DH.diagnosed.soc', 'DH.diagnosed.int',
'refers.soc','refers.int',
'bacassess.soc', 'bacassess.int',
'dxc.soc', 'dxc.int',
'dxb.soc', 'dxb.int',
'bactbtx.soc', 'bactbtx.int',
'PHC.treated.soc','PHC.treated.int',
'DH.treated.soc','DH.treated.int',
'att.soc','att.int',
'cost.soc','cost.int',
costsbystg,
'LYS','LYS0','value',
'deaths.soc','deaths.int')

toget3 <- c('id', 'value', 'tb',
'bacassess.soc', 'bacassess.int',
'dxb.soc', 'dxb.int',
'att.soc', 'att.int',
'PHC.treated.soc', 'PHC.treated.int')

tosum2 <- c(setdiff(toget2,notwt),lyarm)

## containers & loop
allout <- allpout <- allscout <- psapout <- list() #tabular outputs
allout2 <- allpout2 <- allscout2 <- psapout2 <- list() #tabular outputs
ceacl <- NMB <- list()             #CEAC outputs etc

for(cn in isoz){
cat('running model for:',cn, '0-14 years\n')
## --- costs
cols_to_remove <- intersect(cnmz, names(D))
D[, (cols_to_remove) := NULL]
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
out <- out[,lapply(.SD,function(x) sum(x*value, na.rm = T)),.SDcols=tosum,by=id] #sum against popn
## non-incremental cost per ATT
out[,costperATT.soc:=cost.soc/att.soc];
out[,costperATT.int:=cost.int/att.int];
out[,(psoc.sc):=lapply(.SD,function(x) x/att.soc),.SDcols=soc.sc]
out[,costperATT.int:=cost.int/att.int];
out[,(pint.sc):=lapply(.SD,function(x) x/att.int),.SDcols=int.sc]
## increments wrt SOC (per child presenting at either DH/PHC)
out[,Dcost:=cost.int-cost.soc] #inc costs
out[,Datt:=att.int-att.soc] #inc atts
out[,attPC:=1e2*att.int/att.soc] #rel inc atts
out[,Ddeaths:=deaths.int-deaths.soc] #inc deaths
out[,DLYL0:=LYL0.int-LYL0.soc] #inc LYLs w/o discount
out[,DLYL:=LYL.int-LYL.soc] #inc LYLs
## per whatever
out[,DcostperATT:=cost.int/att.int-cost.soc/att.soc];
out[,Dcostperdeaths:=-cost.int/deaths.int+cost.soc/deaths.soc]
out[,DcostperLYS0:=-cost.int/LYL0.int+cost.soc/LYL0.soc]
out[,DcostperLYS:=-cost.int/LYL.int+cost.soc/LYL.soc]
## D/D
out[,DcostperDATT:=Dcost/Datt];
out[,DcostperDdeaths:=-Dcost/Ddeaths]
out[,DcostperDLYS0:=-Dcost/DLYL0]
out[,DcostperDLYS:=-Dcost/DLYL]
## summarize
smy <- outsummary(out)
outs <- smy$outs; pouts <- smy$pouts;
outs[,iso3:=cn]; pouts[,iso3:=cn]
## capture tabular
allout[[cn]] <- outs; allpout[[cn]] <- pouts
## capture data for NMB
NMB[[cn]] <- out[,.(iso3=cn,DLYL,Dcost)]
## ceac data
ceacl[[cn]] <- data.table(iso3=cn,
                          int=make.ceac(out[,.(Q=-DLYL,P=Dcost)],lz),
                          threshold=lz)

## --- grather outcomes for Table 2
out2 <- D[,..toget2]
out2[,c(lyarm):=.(LYS*deaths.soc,LYS*deaths.int,
                  LYS0*deaths.soc,LYS0*deaths.int)] #LYL per pop by arm
## out[,sum(value),by=id]                                       #CHECK
## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
# cntcts <- out2[,sum(value),by=id]
out2 <- out2[,lapply(.SD,function(x) sum(x*value)),.SDcols=tosum2,by=id] #sum against popn
# out2[,contacts.int:=cntcts$V1]
# # out2[,contacts.soc:=2.0] #TODO needs changing
# out2[,contacts.soc:=ifelse(cn=='CMR', 1.202703, 1.941176)] # by country TODO needs checking
# out2[,c('LYL0.soc','LYL0.int'):=NULL] #drop
# out2[,c('prevdeaths.soc','prevdeaths.int'):=.(deaths.soc-incdeaths.soc,deaths.int-incdeaths.int)]
#
out3 <- D[,..toget3]
ttbs <- out3[tb!='noTB',.(ttb=sum(value),tx=sum(att.soc*value)),by=id]
ttbi <- out3[tb!='noTB',.(ttb=sum(value),tx=sum(att.int*value)),by=id]
ftbs <- out3[tb=='noTB',.(ftb=sum(value),tx=sum(att.soc*value)),by=id]
ftbi <- out3[tb=='noTB',.(ftb=sum(value),tx=sum(att.int*value)),by=id]
ftbs.phc <- out3[tb=='noTB',.(fptx=sum(att.soc*value),fptxphc=sum(PHC.treated.soc*value)),by=id]
ftbi.phc <- out3[tb=='noTB',.(fptx=sum(att.int*value),fptxphc=sum(PHC.treated.int*value)),by=id]
## % FP tx at PHC
pfpphc <- data.table(id=ftbs.phc[,(id)], pfp.soc=ftbs.phc[,(fptxphc/fptx)],
                     pfp.int=ftbi.phc[,(fptxphc/fptx)],
                     Dpftb=(ftbi.phc[,fptxphc/fptx]-ftbs.phc[,fptxphc/fptx]))
## % True TB
ptt <- data.table(pttb.soc=ttbs[,(ttb)],pttb.int=ttbi[,(ttb)], Dpttb=(ttbi$ttb-ttbs$ttb))
## % true TB on ATT
pttatt <- data.table(pttbtx.soc=ttbs[,(tx/ttb)],pttbtx.int=ttbi[,(tx/ttb)], Dpttbtx=(ttbi[,tx/ttb]-ttbs[,tx/ttb]))
## % ATT FP
pfpatt <- data.table(pftbtx.soc=ftbs[,(tx/out2$att.soc)],
                     pftbtx.int=(ftbi$tx/out2$att.int),Dpftbtx=(ftbi$tx/out2$att.int-ftbs$tx/out2$att.soc))
# join
ttb <- cbind(pfpphc, ptt,pttatt,pfpatt)
## add to out2
out2 <- merge(out2,ttb,by='id',all.x = TRUE)
#
out2[,assessments.soc:=DH.evaluated.soc+PHC.evaluated.soc];
out2[,assessments.int:=DH.evaluated.int+PHC.evaluated.int];
out2[,passessphc.soc:=PHC.evaluated.soc/assessments.soc];
out2[,passessphc.int:=PHC.evaluated.int/assessments.int];
out2[,pbacassessphc.soc:=(PHC.bacassess.soc)/bacassess.soc]; #TODO change/check
out2[,pbacassessphc.int:=(PHC.bacassess.int)/bacassess.int]; #TODO change/check
out2[,presumptive.soc:=DH.presumptive.soc+PHC.presumptive.soc];
out2[,presumptive.int:=DH.presumptive.int+PHC.presumptive.int];
out2[,dx.soc:=(dxc.soc+dxb.soc)];
out2[,dx.int:=(dxc.int+dxb.int)];
out2[,pdxphc.soc:=(PHC.diagnosed.soc)/dx.soc]; #TODO change/check
out2[,pdxphc.int:=(PHC.diagnosed.int)/dx.int]; #TODO change/check

out2[,pdxb.soc:=dxb.soc/(dxc.soc+dxb.soc)];
out2[,pdxb.int:=dxb.int/(dxc.int+dxb.int)];
# out2[,pdxb.soc:=dxb.soc];
# out2[,pdxb.int:=dxb.int];
out2[,cost.assessments.soc:=DH.evaluated.cost.soc+PHC.evaluated.cost.soc];
out2[,cost.assessments.int:=DH.evaluated.cost.int+PHC.evaluated.cost.int];
out2[,patt.phc.soc:=PHC.treated.soc/att.soc];
out2[,patt.phc.int:=DH.treated.int/att.int];
## increments
out2[,DPHC.presumptive:=PHC.presumptive.int-PHC.presumptive.soc]
out2[,DDH.presumptive:=DH.presumptive.int-DH.presumptive.soc]
out2[,DPHC.evaluated:=PHC.evaluated.int-PHC.evaluated.soc]
out2[,DDH.evaluated:=DH.evaluated.int-DH.evaluated.soc]
out2[,Dassessments:=assessments.int-assessments.soc];
out2[,Dpassessphc:=passessphc.int-passessphc.soc];
out2[,Dbacassess:=bacassess.int-bacassess.soc];
out2[,Dpbacassessphc:=pbacassessphc.int-pbacassessphc.soc];
out2[,Dpresumptive:=presumptive.int-presumptive.soc];
out2[,Drefers:=refers.int-refers.soc];
out2[,Ddx:=dx.int-dx.soc];
out2[,Dpdxphc:=pdxphc.int-pdxphc.soc];
out2[,Dpdxb:=pdxb.int-pdxb.soc];
out2[,DPHC.treated:=PHC.treated.int-PHC.treated.soc]
out2[,DDH.treated:=DH.treated.int-DH.treated.soc]
out2[,Datt:=att.int-att.soc]
out2[,Dpatt.phc:=patt.phc.int-patt.phc.soc]
out2[,Dbactbtx:=bactbtx.int-bactbtx.soc]
out2[,DLYL:=LYL.int-LYL.soc]
out2[,Ddeaths:=deaths.int-deaths.soc]
out2[,DPHC.evaluated.cost:=PHC.evaluated.cost.int-PHC.evaluated.cost.soc]
out2[,DDH.evaluated.cost:=DH.evaluated.cost.int-DH.evaluated.cost.soc]
out2[,DPHC.treated.cost:=PHC.treated.cost.int-PHC.treated.cost.soc]
out2[,DDH.treated.cost:=DH.treated.cost.int-DH.treated.cost.soc]
out2[,Dcost.assessments:=cost.assessments.int-cost.assessments.soc];
out2[,Dcost:=cost.int-cost.soc]
## summarize
smy2 <- Table2(out2) #NOTE set per 1000 index cases or HHs - adjust fac in contact_functions.R
outs2 <- smy2$outs; pouts2 <- smy2$pouts;
outs2[,iso3:=cn]; pouts2[,iso3:=cn]
psapout2[[cn]] <- out2[,iso3:=cn]
## capture tabular
allout2[[cn]] <- outs2; allpout2[[cn]] <- pouts2
}

allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
allscout <- rbindlist(allscout)
ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)
allout2 <- rbindlist(allout2)
allpout2 <- rbindlist(allpout2)
