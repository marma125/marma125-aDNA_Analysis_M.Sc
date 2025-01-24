## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Timothy Farewell
##
## Date Created: 2025-01-24
##
## Copyright (c) Timothy Farewell, 2025
## Email: hello@timfarewell.co.uk
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("~/Google Drive/")      # Tim's working directory (mac)
setwd("C:/Users/tim/Google Drive/")    # Tim's working directory (PC)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on macs.

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
# source("functions/packages.R")       # loads up all the packages we need

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------