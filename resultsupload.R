library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE new sheet
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1tn35_3kI6aRlZIwbyrsFtgRNJsVLjpbSBSVJdPfEBsY/edit?gid=0#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid
){
  filename <- gsub("\\.csv$","",filename) #safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

## upload relevant table data
# upload.to.sheets(here('outdata/'),'table2',shid) #first will need to re-authenticate

# ## rest can be run as block
shidneat <- as.character(as_sheets_id(yourl))

# Main results

## ---- Table 2 -------
Table2 <- fread(gh('outdata/table2.csv')) 

write_sheet(Table2,shidneat,sheet="Tab2RAW")

# Supplementary results

## 0-4 years 
Table2Y <- fread(gh('outdata/table2Y.csv')) 

write_sheet(Table2Y,shidneat,sheet="Tab2YRAW")

## 5-14 years 
Table2O <- fread(gh('outdata/table2O.csv')) 

write_sheet(Table2O,shidneat,sheet="Tab2ORAW")

## --- SA table ---
sa.base <- fread(here('outdata/ICERS.csv'))
sa.hi <- fread(here('outdata/ICERShi.csv'))
sa.lo <- fread(here('outdata/ICERSlo.csv'))
sa.artcov <- fread(here('outdata/ICERSartcov.csv'))
sa.fracphc <- fread(here('outdata/ICERSfracphc.csv'))
sa.truenatbl <- fread(here('outdata/ICERStruenatbl.csv'))
sa.truenatint <- fread(here('outdata/ICERStruenatint.csv'))

# Step 1: Find data frames with a specific pattern
all_objects <- ls()
sa_files <- all_objects[grep("^sa\\.", all_objects)]
print(sa_files)

# Step 2: Add an ID column and bind them together
add_id_column <- function(df_name) {
  df <- get(df_name)
  df <- df %>% mutate(ID = df_name)
  return(df)
}

SAll <- bind_rows(lapply(sa_files, add_id_column))

SAll <- SAll %>% select(ID, everything())

# Print the combined data frame
print(SAll)

id_levels <- c('sa.base',
               'sa.lo','sa.hi',
               'sa.fracphc',
               'sa.artcov',
               'sa.truenatbl',
               'sa.truenatint')

id_labels <- c('Basecase',
               '0% discount rate','5% discount rate',
               'Low PHC presentation',
               'Higher ART coverage',
               'No baseline testing at PHC',
               'Universal Truenat under intervention')

SAll <- SAll %>% 
  mutate(ID = factor(ID, levels = id_levels, labels = id_labels)) %>%
  arrange(ID)

SAll <- SAll %>% 
  rename("Assumption" = ID)

SAll
write_sheet(SAll,shidneat,sheet="SAll")

## --- Main parameters table ---
# ParmsTable1 <- fread(gh('outdata/parameters1.csv')) 
# 
# write_sheet(ParmsTable1,shidneat,sheet="ParmsTab1RAW")
# 
# # Other parameter tables
# flz1 <- c(
#   "tableS1.csv",   "tableS2.csv","tableS2a.csv",
#   "tableS3.csv", "tableS4.csv", "allpout2.csv")
# for( fn in flz1)
#   upload.to.sheets(here('outdata//'),fn,shidneat)

