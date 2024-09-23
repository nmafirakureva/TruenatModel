# Switch different parms on and off to test

# Set DH presented to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1

SOC.F$checkfun(A) # Not OK


# Set DH presented to 1 & SOC DH test to 0
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 0

SOC.F$checkfun(A) # NOTE OK

# Set DH presented to 1 & SOC DH test to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 1

SOC.F$checkfun(A) # NOT OK

# Set DH presented to 1 & SOC DH test to 1 & SOC DH fracUltra to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 1
A$soc.dh.fracUltra <- 1

SOC.F$checkfun(A) # NOTE OK

# Set DH presented to 1 & SOC DH test to 1 & SOC DH fracUltra to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 1
A$soc.dh.fracUltra <- 0

SOC.F$checkfun(A) # NOT OK

# Set DH presented to 1 & SOC DH test to 1 & SOC DH fracUltra to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 1
A$soc.dh.fracUltra <- 0
A$soc.dh.ptbxTrueNat <- 1

SOC.F$checkfun(A) # NOT OK

# Set DH presented to 1 & SOC DH test to 1 & SOC DH fracUltra to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$soc.dh.test <- 1
A$soc.dh.fracUltra <- 1
A$soc.dh.ptbxUltra <- 1
A$prr <- 0
A$dh.ptltfu <- 0

SOC.F$checkfun(A) # NOT OK

# Set DH presented to 1 & SOC DH test to 1 & SOC DH fracUltra to 1
A <- makeTestData(5e3,vrz)
A$dh.presented <- 1
A$dh.presumed <- 1
A$soc.dh.test <- 1
A$soc.dh.fracUltra <- 1
A$soc.dh.ptbxUltra <- 1
A$prr <- 0
A$dh.ptltfu <- 1
A$cfr.tx <- 1

SOC.F$checkfun(A) # NOT OK



names(A)[grep(".tpt$",names(A))]
A <- makeTestData(5e3,vrz)
A$soc.prop.igra.tested <- 0.5
any(round(SOC.F$checkfun(A))!=1) #NOT OKAY

# Filter variables in 'toget' with negative values
negative_cols <- sapply(D[, ..toget2], function(x) any(x < 0))
filtered_cols <- toget2[negative_cols]

# Summary of filtered columns
summary(D[, ..filtered_cols])

summary(D[,.(soc.prop.starting.att)])
names(D)[grepl('igra', names(D))]
checks <- names(D)[grepl('igra', names(D))]
summary(D[,..checks])

checks <- names(D)[grepl('cost.ppd|cost.nhs$', names(D))]
summary(D[,..checks])

# check if sum of cost.ppdsoc and cost.nhs.soc is equal to cost.soc
summary(D[,.(cost.soc, cost.ppd.soc, cost.nhs.soc)])
summary(D[,.(cost.int, cost.ppd.int, cost.nhs.int)])

# calculate proportion of cost.ppd and cost.nhs to cost
summary(D[,.(cost.ppd.soc/cost.soc)])
summary(D[,.(cost.nhs.soc/cost.soc)])
summary(D[,.(cost.ppd.int/cost.int)])
summary(D[,.(cost.nhs.int/cost.int)])

names(SOC.F)

check <- names(SOC.F)[grepl('cost',names(SOC.F))]
# check <- c('costperATT.soc', 'costperATT.int', 'costperTPT.soc', 'costperTPT.int')

check <- gsub('fun','',check)
check_soc <- paste0(check,'.soc')
check_int <- paste0(check,'.int')

summary(D[,..check_soc]) # TPT looks about right but ATT is off?
summary(D[,..check_int]) #should be exactly the same for now

summary(D[,.(screen.soc, tpt.soc, coprevtb.soc, att.soc)])
summary(D[,.(sens.symptom, spec.symptom, sens.any.abn.xray, spec.any.abn.xray)])
