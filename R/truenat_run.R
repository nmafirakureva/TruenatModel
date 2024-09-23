## flags for sensitivity analyses
shell <- FALSE # whether running from shell script or not
if (shell) {
  ## running from shell
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  SA <- args[1] # none,base/lo/hi,tptru,hicoprev
  if (SA == "none") {
    SA <- ""
  }
} else { # set by hand
  rm(list = ls()) # clear all
  shell <- FALSE # whether running from shell script or not
  ## sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'truenatbl' = no Truenat testing @ PHC under SOC
  ## 'truenatint' = Truenat testing @ DH under INT
  ## 'artcov' = higher ART coverage
  ## 'fracphc' = low PHC presented
  sacases <- c("", "lo", "hi", "truenatbl", "truenatint", "artcov", "fracphc")
  SA <- sacases[1]
}

# rm(list=ls())
library(here)
library(tidyverse)

## load other scripts
source(here("R/truenat_pathways_tree.R")) # tree structure and namings: also tree functions & libraries
source(here("R/truenat_functions.R")) # functions for tree parameters

## number of reps
nreps <- 1e3
set.seed(1234)

## attributes to use
hivlevels <- c(0, 1)
artlevels <- c(0, 1)
# tblevels <- c('TB', 'noTB') # Active TB disease, no TB
tblevels <- c("TB+", "TB-", "noTB") # bac confirmable TB, bac unconfirmable TB, not TB
agelevels <- c("0-4", "5-14")

isoz <- c("NGA") # relevant countries

## --- life years and other outputs NOTE needs to be set FALSE on first run thru
LYSdone <- FALSE
if (!LYSdone) {
  ## make discounted life-years if they haven't been done
  LYKc <- GetLifeYears(isolist = isoz, discount.rate = 0.03, yearfrom = 2021)
  LYKc0 <- GetLifeYears(isolist = isoz, discount.rate = 0.00, yearfrom = 2021)
  LYKc5 <- GetLifeYears(isolist = isoz, discount.rate = 0.05, yearfrom = 2021)
  LYKc <- merge(LYKc, LYKc0[, .(iso3, age, LYS0 = LYS)], by = c("iso3", "age"))
  LYKc <- merge(LYKc, LYKc5[, .(iso3, age, LYS5 = LYS)], by = c("iso3", "age"))
  LYK <- LYKc[, .(LYS = mean(LYS), LYS0 = mean(LYS0), LYS5 = mean(LYS5)), by = .(age)] # averaged life-years 4 generic tests
  save(LYKc, file = here("indata/LYKc.Rdata"))
  save(LYK, file = here("indata/LYK.Rdata"))
} else {
  load(file = here("indata/LYKc.Rdata"))
  load(file = here("indata/LYK.Rdata"))
}

# # Sensitivity analysis: 0% & 5% discount rates
if (SA %in% c("hi", "lo")) {
  LYKc[, LYS := ifelse(SA == "lo", LYS0,
    ifelse(SA == "hi", LYS5, LYS)
  )]
}

## prior parameters
PD <- read.csv(here("indata/parms.csv")) # read in main parameters
PDI <- read.csv(here("indata/CASCI.csv")) # read in presentation parameters
PDA <- read.csv(here("indata/CASCA.csv")) # read in age splits
PDC <- read.csv(here("indata/CASCP.csv")) # read in cascade parameters
PDO <- read.csv(here("indata/CASCPO.csv")) # read in outcome parameters
# AD <- read.csv(here('indata/DiagnosticAccuracy.csv')) #read in accuracy parameters
# RD <- fread(gh('indata/RUParms.csv'))    #read resource use data
CDD <- fread(gh("indata/costs.csv")) # read dummy cost data
CD <- fread(gh("indata/model_unit_costs.csv")) # read actual cost data

CET <- fread(gh("indata/CET.csv")) # read cost-effectiveness thresholds

PDCO <- rbind(PDC, PDO, PDI, PDA) # combine cascade and outcome parameters
names(PD)
names(PDCO)

# Quick check & fix
unique(PDCO$parm)
PDCO$parm <- gsub("prusumed", "presumed", PDCO$parm)
PDCO <- PDCO |>
  # filter(grepl('att',parm)) |>
  mutate(parm = case_when(
    age == "0-4" & grepl("att|xray", parm) ~ paste0("u5.", parm),
    age == "5-14" & grepl("att|xray", parm) ~ paste0("o5.", parm),
    .default = parm
  ))

unique(PDCO$parm)
PDCO |>
  filter(grepl("att", parm))

## parameters to be determined from cascade data
PD0 <- PD |>
  mutate(DISTRIBUTION = as.numeric(DISTRIBUTION)) |>
  filter(!is.na(DISTRIBUTION))

## the rest
PD1 <- PD |>
  mutate(dist = as.numeric(DISTRIBUTION)) |>
  filter(is.na(dist)) |>
  select(-dist)

# Adding in HIV/ART data - basically just ART data for now
if (!file.exists(here("outdata/HA.Rdata"))) {
  source(here("R/tb_cdr.R"))
} else {
  load(file = here("outdata/HA.Rdata")) # CDR
}

# merge in WHO HIV/ART data
tmp <- HA[, .(isoz, frac.hiv = hivprop, art.cov = artprop)]
tmp <- melt(tmp, id.vars = "isoz")

tmp <- tmp |>
  rename(NAME = variable, DISTRIBUTION = value)

names(PDCO)
names(PD0)

costparms <- CD |>
  filter(!unit_cost %in% c("cost.art", "cost.phc.referral", "cost.antibiotics", "cost.cxr.exam")) |>
  mutate(
    cost.m = round(cost.m, 2),
    cost.sd = round(cost.sd, 2),
    "Cost (SD)" = paste0(cost.m, " (", cost.sd, ")"),
    facility = ifelse(grepl("phc", unit_cost), "PHC", "Hospital"),
    resource = gsub(" at.*|@PHC |@DH ", "", resource),
    unit_cost = gsub(".phc|.dh", "", unit_cost)
  ) |>
  select(unit_cost, facility, resource, "Cost (SD)") |>
  pivot_wider(names_from = facility, values_from = "Cost (SD)") |>
  mutate(
    PHC = coalesce(PHC, Hospital) # Fill missing PHC values with Hospital values
  )


write.csv(PDCO, file = here("outdata/allparmsFixed.csv"), row.names = FALSE)
write.csv(PD1, file = here("outdata/allparmsDistributions.csv"), row.names = FALSE)
write.csv(costparms, file = here("outdata/costparms.csv"), row.names = FALSE)

unique(PDCO$parm)
PD2 <- rbind(
  PDCO |>
    rename(NAME = parm, DISTRIBUTION = frac) |>
    filter(!grepl("clinbac.assess|TrueNat|TBLamp", NAME)) |> # dropping things not needed for now
    select(NAME, DISTRIBUTION),
  PD0 |>
    select(NAME, DISTRIBUTION),
  tmp |>
    select(NAME, DISTRIBUTION)
) |>
  distinct(NAME, .keep_all = TRUE) |>
  pivot_wider(names_from = NAME, values_from = DISTRIBUTION)

names(PD2)

write.csv(PD1, file = here("indata/allparms.csv"))

# convert into parameter object
P <- parse.parmtable(PD1)

## make base PSA dataset
set.seed(1234) # random number seed

D0 <- makePSA(nreps, P, dbls = list(c("cfrhivor", "cfrartor")))

# merge in fixed parameter
names(PD2)

# drop phc.presented
PD2 <- PD2 |>
  select(-contains(".presented"))
D0 <- cbind(D0, PD2)


if (SA == "artcov") {
  D0[, artcov := art.cov] # higher ART coverage
}

## use these parameters to construct input data by attribute
D0 <- makeAttributes(D0)
D0[, sum(value), by = .(id, tb)] # CHECK
D0[, sum(value), by = .(id, tb)][, mean(V1), by = tb] # CHECK

## read and make cost data
rcsts <- CD
(keep <- vrz[grepl("cost", vrz)])
setdiff(keep, rcsts$unit_cost)
setdiff(rcsts$unit_cost, keep)
(keep2 <- unique(rcsts$unit_cost)[grepl("screening|evaluation|cxr|sample|ATT|antibiotics", unique(rcsts$unit_cost))])
keep <- c(keep, keep2)
rcsts <- rcsts[rcsts$unit_cost %in% keep, ]

names(CD)

rcsts <- setDT(rcsts)

## turn cost data into PSA
rcsts[is.na(rcsts)] <- 0 # some quick fix >> setting NA to 0
rcsts[cost.sd == 0, cost.sd := cost.m / 40] # SD such that 95% UI ~ 10% of mean

allcosts <- rcsts[, .(cost = unit_cost, cost.m, cost.sd)]

C <- MakeCostData(allcosts, nreps) # make cost PSA NOTE using CMR cost data

names(D0)
names(C)

## NOTE can re-run from here to implement changes to MakeTreeParms
## add cost data
D <- merge(D0, C, by = c("id"), all.x = TRUE) # merge into PSA (differentiated D and D0 to facilitate rerunning)

if (SA == "truenatint") {
  # D[,soc.dh.fracUltra:=rbeta(nrow(D),8.251046,5.500698)] # Allowing for Truenat testing @ DH under SOC
  D[, int.dh.fracUltra := rbeta(nrow(D), 8.251046, 5.500698)] # Allowing for Truenat testing @ DH under INT
  D[, int.phc.test := 1] # universal Truenat testing @ PHC under INT
}

if (SA == "truenatbl") {
  D[, soc.phc.test := 0.0] # no PHC test
}

if (SA == "fracphc") {
  D[, phc.presented := 0.40] # low PHC presented
  D[, dh.presented := 1 - phc.presented]
}

## compute other parameters (adds by side-effect)
MakeTreeParms(D, P)

# names(D)[grepl('cost', names(D))]

## checks
D[, sum(value), by = .(id)] # CHECK
D[, sum(value), by = .(id, age)] # CHECK
D[, sum(value), by = .(id, age, tb)] # CHECK

## check for leaks
head(SOC.F$checkfun(D)) # SOC arm
head(INT.F$checkfun(D)) # INT arm

names(SOC.F)

# TODO: check if this is OKAY
# Approach to normalizing everything to presumptive TB
D[, phc.presented := phc.presented / phc.presumed]
D[, dh.presented := dh.presented / dh.presumed]
D[, phc.presumed := phc.presumed / phc.presumed]
D[, dh.presumed := dh.presumed / dh.presumed]

D |>
  select(tb, phc.presented, dh.presented, phc.presumed, dh.presumed) |>
  group_by(tb) |>
  summarise_all(mean)

## === RUN MODEL
arms <- c("SOC", "INT")
D <- runallfuns(D, arm = arms) # appends anwers

## restricted trees:
D[["soc_att_check"]] <- SOC.att.F$checkfun(D)
D[["soc_att_cost"]] <- SOC.att.F$costfun(D)
D[["soc_notx_check"]] <- SOC.notx.F$checkfun(D)
D[["soc_notx_cost"]] <- SOC.notx.F$costfun(D)

D[["int_att_check"]] <- INT.att.F$checkfun(D)
D[["int_att_cost"]] <- INT.att.F$costfun(D)
D[["int_notx_check"]] <- INT.notx.F$checkfun(D)
D[["int_notx_cost"]] <- INT.notx.F$costfun(D)

## cross checks: compute probability of endpoints from whole tree vs pruned trees
head(D[, int_att_check])
head(D[, attend.int])

## NOTE OK
all(D[, attend.int] == D[, int_att_check])
all(D[, attend.soc] == D[, soc_att_check])

## NOTE confusingly there is also a notxend variable, which is different
all(D[, notx.int] == D[, int_notx_check])
all(D[, notx.soc] == D[, soc_notx_check])

D[, table(tb)]

## create restricted PSA
DR <- D[, .(
  id, tb,
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
DR[, c(
  "soc_att_cost",
  "soc_notx_cost",
  "int_att_cost",
  "int_notx_cost"
) := .(
  soc_att_cost / soc_att_check,
  soc_notx_cost / soc_notx_check,
  int_att_cost / int_att_check,
  int_notx_cost / int_notx_check
)]

save(DR, file = here("outdata/DR.Rdata"))

## summary
DRS <- DR[, lapply(.SD, mean), .SDcols = names(DR)[-c(1, 2)], by = tb]
DRS <- melt(DRS, id = "tb")
DRS[, c("arm", "outcome", "quantity") := tstrsplit(variable, split = "_")]
options(scipen = 999)
(DRS <- dcast(data = DRS, formula = arm + quantity + outcome ~ tb, value.var = "value"))
fwrite(DRS, file = here("outdata/DRS.csv"))

(dxb <- names(D)[grepl("dxb|dxc|att.soc|att.int", names(D))])
summary(D[tb == "TB", ..dxb])
summary(D[tb == "noTB", ..dxb])

## Cost-effectiveness analysis
## --- run over different countries
cnmz <- names(C)
cnmz <- cnmz[cnmz != "id"]

names(labz)[grepl("cost$", names(labz))]
costsbystg <- c(
  "DH.presumptive.cost.soc", "DH.evaluated.cost.soc", "DH.treated.cost.soc",
  "DH.presumptive.cost.int", "DH.evaluated.cost.int", "DH.treated.cost.int",
  "PHC.presumptive.cost.soc", "PHC.evaluated.cost.soc", "PHC.treated.cost.soc",
  "PHC.presumptive.cost.int", "PHC.evaluated.cost.int", "PHC.treated.cost.int"
)

# outcomesbystg <- gsub('.cost','',costsbystg)
toget <- c(
  "id",
  "cost.soc", "cost.int",
  costsbystg,
  "att.soc", "att.int",
  # 'evaluated.soc','evaluated.int',
  "refers.soc", "refers.int",
  "dxc.soc", "dxc.int",
  "dxb.soc", "dxb.int",
  "deaths.soc", "deaths.int",
  "LYS", "LYS0", "value"
)
notwt <- c("id", "LYS", "LYS0", "value") # variables not to weight against value
lyarm <- c("LYL.soc", "LYL.int")
lyarm <- c(lyarm, gsub("\\.", "0\\.", lyarm)) # include undiscounted
tosum <- c(setdiff(toget, notwt), lyarm)

toget2 <- c(
  "id",
  "PHC.presumptive.soc", "PHC.presumptive.int",
  "PHC.evaluated.soc", "PHC.evaluated.int",
  "PHC.bacassess.soc", "PHC.bacassess.int",
  "PHC.diagnosed.soc", "PHC.diagnosed.int",
  "DH.presumptive.soc", "DH.presumptive.int",
  "DH.evaluated.soc", "DH.evaluated.int",
  "DH.bacassess.soc", "DH.bacassess.int",
  "DH.diagnosed.soc", "DH.diagnosed.int",
  "refers.soc", "refers.int",
  "bacassess.soc", "bacassess.int",
  "dxc.soc", "dxc.int",
  "dxb.soc", "dxb.int",
  "bactbtx.soc", "bactbtx.int",
  "PHC.treated.soc", "PHC.treated.int",
  "DH.treated.soc", "DH.treated.int",
  "att.soc", "att.int",
  "cost.soc", "cost.int",
  costsbystg,
  "LYS", "LYS0", "value",
  "deaths.soc", "deaths.int"
)

toget3 <- c(
  "id", "value", "tb",
  "bacassess.soc", "bacassess.int",
  "dxb.soc", "dxb.int",
  "att.soc", "att.int",
  "PHC.treated.soc", "PHC.treated.int"
)

tosum2 <- c(setdiff(toget2, notwt), lyarm)

## heuristic to scale top value for thresholds:
heur <- c("id", "value", "deaths.int", "deaths.soc")
out <- D[, ..heur]
out <- out[, lapply(.SD, function(x) sum(x * value)), .SDcols = c("deaths.int", "deaths.soc"), by = id] # sum against popn
## topl <- 0.25/out[,mean(deaths.soc-deaths.iph)]
topl <- 500 # 100
lz <- seq(from = 0, to = topl, length.out = 1000) # threshold vector for CEACs
## staged costs by arm
soc.sc <- grep("soc", costsbystg, value = TRUE)
psoc.sc <- paste0("perATT.", soc.sc)
int.sc <- grep("int", costsbystg, value = TRUE)
pint.sc <- paste0("perATT.", int.sc)


## containers & loop
allout <- allpout <- allscout <- psapout <- list() # tabular outputs
allout2 <- allpout2 <- allscout2 <- psapout2 <- list() # tabular outputs
ceacl <- NMB <- list() # CEAC outputs etc

for (cn in isoz) {
  dc <- D[isoz == cn]
  cat("running model for:", cn, "0-14 years\n")
  ## --- costs
  cols_to_remove <- intersect(cnmz, names(D))
  dc[, (cols_to_remove) := NULL]
  ## add cost data
  C <- MakeCostData(allcosts, nreps) # make cost PSA
  dc <- merge(dc, C, by = "id", all.x = TRUE) # merge into PSA
  ## --- DALYs
  ## drop any that are there
  if ("LYS" %in% names(dc)) dc[, c("LYS", "LYS0") := NULL]
  dc <- merge(dc, LYKc[, .(age, LYS, LYS0)], by = "age", all.x = TRUE) # merge into PSA
  ## --- run model (quietly)
  invisible(capture.output(dc <- runallfuns(dc, arm = arms)))
  ## --- grather outcomes
  out <- dc[, ..toget]
  out[, c(lyarm) := .(
    LYS * deaths.soc, LYS * deaths.int,
    LYS0 * deaths.soc, LYS0 * deaths.int
  )] # LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  out <- out[, lapply(.SD, function(x) sum(x * value, na.rm = T)), .SDcols = tosum, by = id] # sum against popn
  ## non-incremental cost per ATT
  out[, costperATT.soc := cost.soc / att.soc]
  out[, costperATT.int := cost.int / att.int]
  out[, (psoc.sc) := lapply(.SD, function(x) x / att.soc), .SDcols = soc.sc]
  out[, costperATT.int := cost.int / att.int]
  out[, (pint.sc) := lapply(.SD, function(x) x / att.int), .SDcols = int.sc]
  ## increments wrt SOC (per child presenting at either DH/PHC)
  out[, Dcost := cost.int - cost.soc] # inc costs
  out[, Datt := att.int - att.soc] # inc atts
  out[, attPC := 1e2 * att.int / att.soc] # rel inc atts
  out[, Ddeaths := deaths.int - deaths.soc] # inc deaths
  out[, DLYL0 := LYL0.int - LYL0.soc] # inc LYLs w/o discount
  out[, DLYL := LYL.int - LYL.soc] # inc LYLs
  ## per whatever
  out[, DcostperATT := cost.int / att.int - cost.soc / att.soc]
  out[, Dcostperdeaths := -cost.int / deaths.int + cost.soc / deaths.soc]
  out[, DcostperLYS0 := -cost.int / LYL0.int + cost.soc / LYL0.soc]
  out[, DcostperLYS := -cost.int / LYL.int + cost.soc / LYL.soc]
  ## D/D
  out[, DcostperDATT := Dcost / Datt]
  out[, DcostperDdeaths := -Dcost / Ddeaths]
  out[, DcostperDLYS0 := -Dcost / DLYL0]
  out[, DcostperDLYS := -Dcost / DLYL]
  ## summarize
  smy <- outsummary(out)
  outs <- smy$outs
  pouts <- smy$pouts
  outs[, iso3 := cn]
  pouts[, iso3 := cn]
  ## capture tabular
  allout[[cn]] <- outs
  allpout[[cn]] <- pouts
  ## capture data for NMB
  NMB[[cn]] <- out[, .(iso3 = cn, DLYL, Dcost)]
  ## ceac data
  ceacl[[cn]] <- data.table(
    iso3 = cn,
    int = make.ceac(out[, .(Q = -DLYL, P = Dcost)], lz),
    threshold = lz
  )

  ## --- grather outcomes for Table 2
  out2 <- dc[, ..toget2]
  out2[, c(lyarm) := .(
    LYS * deaths.soc, LYS * deaths.int,
    LYS0 * deaths.soc, LYS0 * deaths.int
  )] # LYL per pop by arm
  ## out[,sum(value),by=id]                                       #CHECK
  ## TODO: check why the above is generating some NAs. Quick fix is using na.rm=TRUE
  # cntcts <- out2[,sum(value),by=id]
  out2 <- out2[, lapply(.SD, function(x) sum(x * value)), .SDcols = tosum2, by = id] # sum against popn

  out3 <- dc[, ..toget3]
  ttbs <- out3[tb != "noTB", .(ttb = sum(value), tx = sum(att.soc * value)), by = id] # true TB SOC
  ttbi <- out3[tb != "noTB", .(ttb = sum(value), tx = sum(att.int * value)), by = id] # true TB INT
  ftbs <- out3[tb == "noTB", .(ftb = sum(value), tx = sum(att.soc * value)), by = id] # false TB SOC
  ftbi <- out3[tb == "noTB", .(ftb = sum(value), tx = sum(att.int * value)), by = id] # false TB INT
  ftbs.phc <- out3[tb == "noTB", .(fptx = sum(att.soc * value), fptxphc = sum(PHC.treated.soc * value)), by = id] # false TB @ PHC SOC
  ftbi.phc <- out3[tb == "noTB", .(fptx = sum(att.int * value), fptxphc = sum(PHC.treated.int * value)), by = id] # false TB @ PHC INT
  ## % FP tx at PHC
  pfpphc <- data.table(
    id = ftbs.phc[, (id)], pfp.soc = ftbs.phc[, (fptxphc / fptx)],
    pfp.int = ftbi.phc[, (fptxphc / fptx)],
    Dpftb = (ftbi.phc[, fptxphc / fptx] - ftbs.phc[, fptxphc / fptx])
  )
  ## % True TB
  ptt <- data.table(pttb.soc = ttbs[, (ttb)], pttb.int = ttbi[, (ttb)], Dpttb = (ttbi$ttb - ttbs$ttb))
  ## % true TB on ATT
  pttatt <- data.table(pttbtx.soc = ttbs[, (tx / ttb)], pttbtx.int = ttbi[, (tx / ttb)], Dpttbtx = (ttbi[, tx / ttb] - ttbs[, tx / ttb]))
  ## % ATT FP
  pfpatt <- data.table(
    pftbtx.soc = ftbs[, (tx / out2$att.soc)],
    pftbtx.int = (ftbi$tx / out2$att.int), Dpftbtx = (ftbi$tx / out2$att.int - ftbs$tx / out2$att.soc)
  )
  # join
  ttb <- cbind(pfpphc, ptt, pttatt, pfpatt)
  ## add to out2
  out2 <- merge(out2, ttb, by = "id", all.x = TRUE)
  #
  out2[, assessments.soc := DH.evaluated.soc + PHC.evaluated.soc]
  out2[, assessments.int := DH.evaluated.int + PHC.evaluated.int]
  out2[, passessphc.soc := PHC.evaluated.soc / assessments.soc]
  out2[, passessphc.int := PHC.evaluated.int / assessments.int]
  out2[, pbacassessphc.soc := (PHC.bacassess.soc) / bacassess.soc] # TODO change/check
  out2[, pbacassessphc.int := (PHC.bacassess.int) / bacassess.int] # TODO change/check
  out2[, presumptive.soc := DH.presumptive.soc + PHC.presumptive.soc]
  out2[, presumptive.int := DH.presumptive.int + PHC.presumptive.int]
  out2[, dx.soc := (dxc.soc + dxb.soc)]
  out2[, dx.int := (dxc.int + dxb.int)]
  out2[, pdxphc.soc := (PHC.diagnosed.soc) / dx.soc] # TODO change/check
  out2[, pdxphc.int := (PHC.diagnosed.int) / dx.int] # TODO change/check

  out2[, pdxb.soc := dxb.soc / (dxc.soc + dxb.soc)]
  out2[, pdxb.int := dxb.int / (dxc.int + dxb.int)]
  # out2[,pdxb.soc:=dxb.soc];
  # out2[,pdxb.int:=dxb.int];
  out2[, cost.assessments.soc := DH.evaluated.cost.soc + PHC.evaluated.cost.soc]
  out2[, cost.assessments.int := DH.evaluated.cost.int + PHC.evaluated.cost.int]
  out2[, patt.phc.soc := PHC.treated.soc / att.soc]
  out2[, patt.phc.int := DH.treated.int / att.int]
  out2[, patt.bac.soc := bactbtx.soc / att.soc]
  out2[, patt.bac.int := bactbtx.int / att.int]

  ## increments
  out2[, DPHC.presumptive := PHC.presumptive.int - PHC.presumptive.soc]
  out2[, DDH.presumptive := DH.presumptive.int - DH.presumptive.soc]
  out2[, DPHC.evaluated := PHC.evaluated.int - PHC.evaluated.soc]
  out2[, DDH.evaluated := DH.evaluated.int - DH.evaluated.soc]
  out2[, Dassessments := assessments.int - assessments.soc]
  out2[, Dpassessphc := passessphc.int - passessphc.soc]
  out2[, Dbacassess := bacassess.int - bacassess.soc]
  out2[, Dpbacassessphc := pbacassessphc.int - pbacassessphc.soc]
  out2[, Dpresumptive := presumptive.int - presumptive.soc]
  out2[, Drefers := refers.int - refers.soc]
  out2[, Ddx := dx.int - dx.soc]
  out2[, Dpdxphc := pdxphc.int - pdxphc.soc]
  out2[, Dpdxb := pdxb.int - pdxb.soc]
  out2[, DPHC.treated := PHC.treated.int - PHC.treated.soc]
  out2[, DDH.treated := DH.treated.int - DH.treated.soc]
  out2[, Datt := att.int - att.soc]
  out2[, Dpatt.phc := patt.phc.int - patt.phc.soc]
  out2[, Dpatt.bac := patt.bac.int - patt.bac.soc]
  out2[, DLYL := LYL.int - LYL.soc]
  out2[, Ddeaths := deaths.int - deaths.soc]
  out2[, DPHC.evaluated.cost := PHC.evaluated.cost.int - PHC.evaluated.cost.soc]
  out2[, DDH.evaluated.cost := DH.evaluated.cost.int - DH.evaluated.cost.soc]
  out2[, DPHC.treated.cost := PHC.treated.cost.int - PHC.treated.cost.soc]
  out2[, DDH.treated.cost := DH.treated.cost.int - DH.treated.cost.soc]
  out2[, Dcost.assessments := cost.assessments.int - cost.assessments.soc]
  out2[, Dcost := cost.int - cost.soc]
  ## summarize
  smy2 <- Table2(out2) # NOTE set per 1000 index cases or HHs - adjust fac in contact_functions.R
  outs2 <- smy2$outs
  pouts2 <- smy2$pouts
  outs2[, iso3 := cn]
  pouts2[, iso3 := cn]
  psapout2[[cn]] <- out2[, iso3 := cn]
  ## capture tabular
  allout2[[cn]] <- outs2
  allpout2[[cn]] <- pouts2
}

allout <- rbindlist(allout)
allpout <- rbindlist(allpout)
# allscout <- rbindlist(allscout)
ceacl <- rbindlist(ceacl)
NMB <- rbindlist(NMB)
allout2 <- rbindlist(allout2)
allpout2 <- rbindlist(allpout2)

fwrite(allout, file = gh("outdata/allout") + SA + ".csv")
fwrite(allpout, file = gh("outdata/allpout") + SA + ".csv")
fwrite(allout2, file = gh("outdata/allout2") + SA + ".csv")
fwrite(allpout2, file = gh("outdata/allpout2") + SA + ".csv")
# fwrite(allscout,file=gh('outdata/allscout') + SA + '.csv')
save(ceacl, file = gh("outdata/ceacl") + SA + ".Rdata")
save(NMB, file = gh("outdata/NMB") + SA + ".Rdata")

ICERS <- allpout2[, .(DLYL, Dcost, ICER = ICER)]
fwrite(ICERS, file = gh("outdata/ICERS") + SA + ".csv")
ICERS

# summary(out2[,.(Dpatt.bac*100,patt.bac.soc*100,patt.bac.int*100)])
# out2[,median(bactbtx.soc/att.soc)]*100
# out2[,median(bactbtx.int/att.int)]*100
#
# ## checks
# out[,.(att.int/att.soc)]
# out[,.(Datt/att.soc)]
# out[,.(att.soc,att.int)]
# out[,.(Ddeaths/Datt)] #OK
# out[,.(DLYL/Ddeaths)] #OK
# ## check
# allout[,.(costperATT.int.mid-costperATT.soc.mid,DcostperATT.mid)]

## CEAC plot
cbPalette <- c("#009E73")
ceaclm <- melt(ceacl, id = c("iso3", "threshold"))
ceaclm[, Intervention := ifelse(variable == "int", "Intervention", "Standard of care")]
## name key
ckey <- data.table(
  iso3 = c("NGA"),
  country = c("Nigeria")
)

ceaclm <- merge(ceaclm, ckey, by = "iso3", all.x = TRUE)

## plot: int only
GP <- ggplot(
  ceaclm[variable == "int" &
    iso3 %in% c("NGA")],
  aes(threshold, value,
    col = country, lty = Intervention
  )
) +
  # guides(color = "none") +
  geom_line(show.legend = F) +
  theme_classic(base_family = "Times-Roman") +
  theme(legend.position = "top", legend.title = element_blank()) +
  ggpubr::grids() +
  ylab("Probability cost-effective") +
  xlab("Cost-effectiveness threshold (USD/DALY averted)") +
  scale_colour_manual(values = cbPalette)
GP

ggsave(GP, file = gh("plots/CEAC") + SA + ".png", w = 7, h = 5)

# # CE plane
DS <- dc[, .(
  cost.SOC = sum(cost.soc * value),
  cost.INT = sum(cost.int * value),
  lyl.SOC = sum(deaths.soc * value * LYS),
  lyl.INT = sum(deaths.int * value * LYS)
),
by = id
] # PSA summary

DS[, dDALY := lyl.SOC - lyl.INT]
DS[, dCost := cost.INT - cost.SOC]

CET <- CET |>
  rename(CET = threshold) |>
  mutate(CET = factor(CET, levels = c("0.5x GDP", "1x GDP")))

icer <- DS %>%
  summarise(
    cst = mean(dCost),
    eff = mean(dDALY),
    icer = format(round(mean(dCost) / mean(dDALY)), big.mark = ",", scientific = FALSE)
  )

GP <- ggplot(DS, aes(dDALY, dCost)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.1, shape = 1) +
  geom_point(
    data = icer, aes(y = cst, x = eff),
    size = 2, color = "red", shape = 20
  ) +
  # geom_text(data = icer, aes(y = cst, x = 0.7, label=paste0('ICER=', 'US$',icer, '/DALY averted')), size = 3)+
  geom_abline(
    data = CET[CET %in% c("0.5x GDP", "1x GDP")],
    aes(intercept = 0, slope = value, col = CET)
  ) +
  scale_color_manual(
    name = "",
    values = c("#00BFC4", "#F8766D"),
    labels = c("0.5x GDP", "1x GDP")
  ) +
  # facet_wrap(~country, scales = 'free') +
  scale_y_continuous(label = comma) +
  xlab("Discounted disability adjusted life-years averted") +
  ylab("Incremental cost (USD)") +
  theme(legend.position = "top")
if (!shell) GP

fn1 <- glue(here("plots/CEplane")) + SA + ".png"
## fn2 <- glue(here('plots/CEplane')) + SAT + '.pdf'
ggsave(GP, file = fn1, w = 6, h = 5)
## ggsave(GP,file=fn2,w=10,h=10)
## PDF versions too big
#
#
ICERS

allout$ICER

ceaclm$threshold[ceaclm$value >= 0.5][1]
ceaclm$threshold[ceaclm$value == max(ceaclm$value)][1]
