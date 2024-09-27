#!/bin/bash
# NOTE (after changing shell flag in truenat_run.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, truenatbl (no Truenat testing @ PHC under SOC), truenatint (Truenat testing @ DH under INT), artcov (higher ART coverage), fracphc (low PHC presented)

R --slave --vanilla --args <truenat_run.R hi & R --slave --vanilla --args <truenat_run.R lo &
R --slave --vanilla --args <truenat_run.R truenatbl & R --slave --vanilla --args <truenat_run.R truenatint &
R --slave --vanilla --args <truenat_run.R artcov & R --slave --vanilla --args <truenat_run.R fracphc &
R --slave --vanilla --args <truenat_run.R none


