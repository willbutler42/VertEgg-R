####################################################################################
#################### Runner for VertEgg R package ##################################
####################################################################################

##############################################
## R scripts #################################
##############################################

## Global variables and egg distribution initialisation
source("R.initialise\\ve_grid.r")
source("R.initialise\\spawn.r") 
source("R.initialise\\ve_rand.r")

## Stationary coefficients
source("R.functions\\eggsact.r")
source("R.functions\\srcsact.r")
source("R.functions\\sstate.r")

## Terminal velocity
source("R.functions\\eggvelst.r")
source("R.functions\\eggvel.r")
source("R.functions\\molvisk.r")
source("R.functions\\dens0.r")

## Analysis tools
source("R.functions\\ve_int.r")
source("R.functions\\ve_mean.r")
source("R.functions\\ve_std.r")
source("R.functions\\eggmom.r")
source("R.functions\\ve_drint.r")
source("R.functions\\ve_rmsd.r")
source("R.functions\\ve_quantile.r")

## Transient problems
source("R.functions\\lwendrof.r")
source("R.functions\\ftcs.r")
source("R.functions\\minlim.r")
source("R.functions\\posmet.r")
source("R.functions\\upstream.r")

################################################
############## densities and diameters #########
################################################

# load(file="data\\csh.RData")

################################################
############## Variables #######################
################################################

## Gravitational acceleration
g <- 9.81

## Standard density of seawater
rho_s <- 1025

## StartDates & endDates of spawning gadoids in Iceland
cs <- as.Date("15/03/2006", format="%d/%m/%Y")
ce <- as.Date("20/05/2006", format="%d/%m/%Y")
hs <- as.Date("31/03/2006", format="%d/%m/%Y")
he <- as.Date("31/05/2006", format="%d/%m/%Y")
ss <- as.Date("01/02/2006", format="%d/%m/%Y")
se <- as.Date("15/04/2006", format="%d/%m/%Y")





