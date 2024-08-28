#!/bin/bash
# NOTE (after changing shell flag in truenat_run.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, cdr (higher cdr for incidence), txd (completion of ATT/TPT included)

R --slave --vanilla --args <truenat_run.R hi & R --slave --vanilla --args <truenat_run.R lo &
R --slave --vanilla --args <truenat_run.R truenatbl & R --slave --vanilla --args <truenat_run.R truenatint &
R --slave --vanilla --args <truenat_run.R artcov & R --slave --vanilla --args <truenat_run.R fracphc &
R --slave --vanilla --args <truenat_run.R none


