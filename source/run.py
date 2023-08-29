#!/usr/bin/env python

import main
import sys

# this module calls the main module to run the FORGE model. Certain model parameters are varied, and are passed as one positional argument to this script. The parameters are separated by commas as follows: human_birthrate,maximum_human_fat_to_body_weight_ratio,human_bodyweight,metabolic_rate_of_anabolism,metabolic_rate_of_catabolism,animal_background_mortality_rate,human_background_mortality_rate.The variables pathin and pathout are the paths to the location of the input and output netCDF files, and can be modified as desired. The variables resand run_year are time resolution and number of years of run and can also be modified as desired.

#********************************************************
#******************** RUN PARAMETERS ********************

#pathin ='/archive/jlipper/GPP/INPUT/gpp/' 
pathin ='/archive/jlipper/GPP/INPUT/gpp/' # the path where the input files are stored
pathout='/archive/jlipper/GPP/OUTPUT/' # the path where the output file will be generated
res='year'  # temporal resolution in the output file: 'year','month','daily'
run_year=11  # run how many years

ninput_years=3 # number of files from the INPUT folder used. Each file corresponds to a year of GPP data starting in 2001. In this version, there are 19 files available and the same prec and t2m data is used for all files.
#filein_pre='globe_S0_'  # for a single point (42N, 123W) as a test. As for a global 2-degree simulation with run_year=300, the cpu-time is about 6 hours with 8 cores in parallel. 
filein_pre='dailydata_' # common beginning of all input file names (without year) 
fileout=f'{sys.argv[1]}'  # name of the output file

run_pars = {'pathin':pathin, 'pathout':pathout, 'res':res, 'run_year':run_year, 'ninput_years':ninput_years, 'filein_pre':filein_pre, 'fileout':fileout}

#********************************************************
if len(sys.argv)>2: # checks if positional argument containing varied parameters is given
     parameters=sys.argv[2].split(',') # splits positional argument containing varied parameters into list

     for i in range(len(parameters)): # converts all varied parameters given as fractions to float form
          if '/' in parameters[i]:
               n,d=map(float, parameters[i].split('/'))
               parameters[i]=n/d
else:
     parameters='none' # for if no model parameters are given as positional parameters

main.run(pathin,pathout,res,run_year,ninput_years,filein_pre,fileout,parameters,run_pars) # runs FORGE model by calling main.py module
