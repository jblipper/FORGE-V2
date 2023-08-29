#!/usr/bin/env python

import main
import sys

# this module calls the main module to run the FORGE model. Certain model parameters are varied, and are passed as one positional argument to this script. The parameters are separated by commas as follows: human_birthrate,maximum_human_fat_to_body_weight_ratio,human_bodyweight,metabolic_rate_of_anabolism,metabolic_rate_of_catabolism,animal_background_mortality_rate,human_background_mortality_rate.The variables pathin and pathout are the paths to the location of the input and output netCDF files, and can be modified as desired. The variables resand run_year are time resolution and number of years of run and can also be modified as desired.

#pathin ='/archive/jlipper/GPP/INPUT/gpp/' 
pathin ='/archive/jlipper/GPP/INPUT/gpp/' # the path where the input files are stored
pathout='/archive/jlipper/GPP/OUTPUT/' # the path where the output file will be generated
res='year'  # temporal resolution in the output file: 'year','month','daily'
run_year=75  # run how many years

pars=[0.015,0.1,1.0,200,2.5] # default parameter values,corresponding to ΦV,ΦA,emax,c,q in Table 1. 
phi_v =pars[0] # edible fractions of vegetable biomass to humans (ΦV)
phi_a =pars[1] # edible fracvtion of animal biomass to humans (ΦA)
eG_cst=pars[2] # maximum efficiency of gathering (emax)
eH_cst=pars[2] # maximum efficiency of hunting (emax)
c_dep =pars[3] # patch depletion rate (c)
q     =pars[4] # parameter used to tune meat craving (q)

#filein_pre='globe_S0_'  # for a single point (42N, 123W) as a test. As for a global 2-degree simulation with run_year=300, the cpu-time is about 6 hours with 8 cores in parallel. 
filein_pre='dailydata_' # common beginning of all input file names (without year) 
fileout=f'{sys.argv[1]}'  # name of the output file
parameters=sys.argv[2].split(',') # splits positional argument containing varied parameters into list

n5,d5=map(float, parameters[5].split('/')) # converts background animal mortality rate from fraction to float
parameters[5]=n5/d5
n6,d6=map(float, parameters[6].split('/')) # converts background human mortality rate from fraction to float
parameters[6]=n6/d6

main.run(pathin,pathout,res,run_year,phi_v,phi_a,eG_cst,eH_cst,c_dep,q,filein_pre,fileout,parameters) # runs FORGE model by calling main.py module
