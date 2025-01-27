# rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(data.table)

gh <- function(x) glue(here(x))

## === making CDRs and getting some relevant HIV data ===
isoz <- 'NGA' #Nigeria
whodir <- glue(here('who_data/data2023//'))

## === load WHO notifications
## http://www.who.int/tb/country/data/download/en/
fn <- whodir+'TB_notifications_2024-07-23.csv'
N <- fread(fn)

nmz <- paste0('newrel_',c('m04','f04',
                          'm59','f59',
                          'm1014','f1014',
                          'm014','f014'))
nmz <- c('iso3','year',nmz)

## reduce to relevant data: 2018 has no NA
NP <- N[year==max(year),..nmz]
NP <- melt(NP,id.vars = c('iso3','year'))
NP[,sex:=ifelse(grepl("f",variable),'F','M')]
NP[,age:=gsub("[a-z]|_","",variable)]
NP <- NP[iso3 %in% isoz]
NP
NP <- NP[age %in% c('04','014')]
NP

## === load WHO age-specific incidence estimates
fn <- whodir+'TB_burden_age_sex_2024-07-23.csv'
A <- fread(fn)

## keep only relevant categories
A <- A[year==max(year),]
A <- A[sex!='a']
A <- A[age_group %in% c('0-4','0-14','5-14')]
A <- A[risk_factor=='all']
A[,age:=gsub('-','',age_group)]
## harmonize namings
A[sex=='f',sex:='F']
A[sex=='m',sex:='M']
unique(A[,.(sex,age)])                  #check
A[,best.sd:=(hi-lo)/3.92]
A <- A[iso3 %in% isoz]
A

## HIV
fn <- here('who_data/data2024/TB_burden_countries_2025-01-17.csv')
H <- fread(fn)
H <- H[year==unique(A$year),.(iso3,e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]
H <- H[iso3 %in% isoz]
H[,hiv:=e_tbhiv_prct/100]
H[,hiv.sd:=(e_tbhiv_prct_hi-e_tbhiv_prct_lo)/392] #adults of course
H <- H[,.(iso3,hiv,hiv.sd)]
H

hfn <- here('outdata/H.Rdata')
save(H,file=hfn)

## HIV/ART
fn <- whodir+'TB_notifications_2024-07-23.csv'
HA <- fread(fn)
HA <- HA[iso3 %in% isoz & year==max(year),.(iso3,newrel_hivtest,newrel_hivpos,newrel_art)]
HA[,.(iso3,newrel_hivpos/newrel_hivtest,newrel_art/newrel_hivpos)]
HA[,hivprop:=newrel_hivpos/newrel_hivtest]
HA[,artprop:=newrel_art/newrel_hivpos]
HA

hfn <- here('outdata/HA.Rdata')
save(HA,file=hfn)

## === merge data
AN <- merge(NP[,.(iso3,sex,age,notes=value)],
            A[,.(iso3,sex,age,inc=best,lo,hi)],
            by=c('iso3','sex','age'),all.x=TRUE,all.y=FALSE)


ANO <- AN[age=='014']
ANY <- AN[age=='04']
ANB <- merge(ANY[,.(iso3,sex,notes.04=notes,inc.04=inc,
                    inc.sd.04=(hi-lo)/3.92)],
             ANO[,.(iso3,sex,notes.514=notes,inc.514=inc,
                    inc.sd.514=(hi-lo)/3.92)],
             by=c('iso3','sex')
)
CDRu <- ANB[,.(rel.sd=mean(inc.sd.514/inc.514)),by=iso3]
ANB[,c('notes.514','inc.514'):=.(notes.514-notes.04,inc.514-notes.04)]
CDR <- ANB[,.(notes.04=sum(notes.04),inc.04=sum(inc.04),
              notes.514=sum(notes.514),inc.514=sum(inc.514)),by=iso3]
CDR[,c('cdr04','cdr514'):=.(notes.04/inc.04,notes.514/inc.514)]
CDR[,totnotes:=notes.04+notes.514]
CDR[,c('frac04','frac514'):=.(notes.04/totnotes,notes.514/totnotes)]
CDR <- melt(CDR[,.(iso3,frac04,frac514,cdr04,cdr514)],id='iso3')
CDR[,qty:='cdr']
CDR[grepl('frac',variable),qty:='frac']
CDR[,age:='0-4']
CDR[grepl('14',variable),age:='5-14']
CDR <- CDR[,.(iso3,qty,age,value)]
CDR <- merge(CDR,CDRu,by='iso3',all.x=TRUE)
CDR[,cdr.v:=NA_real_]
CDR[,cdr.v:=(rel.sd*value)^2]
CDR[,rel.sd:=NULL]
CDR <- CDR[iso3 %in% isoz]

cdrfn <- here('outdata/CDR.Rdata')
save(CDR,file=cdrfn)

