#!/usr/bin/env python

import main
import sys

# used to test the forge model. The input file consist of four pixels to run the model on

#********************************************************
#******************** RUN PARAMETERS ********************

pathin ='/archive/jlipper/GPP/INPUT/gpptest/' # the path where the input files are stored
pathout='/archive/jlipper/GPP/OUTPUT/' # the path where the output file will be generated
res='year'  # temporal resolution in the output file: 'year','month','daily'
run_year=30  # run how many years

ninput_years=19 # number of files from the INPUT folder used. Each file corresponds to a year of GPP data starting in 2001. In this version, there are 19 files available and the same prec and t2m data is used for all files.
#filein_pre='globe_S0_'  # for a single point (42N, 123W) as a test. As for a global 2-degree simulation with run_year=300, the cpu-time is about 6 hours with 8 cores in parallel. 
filein_pre='dailydatatest_' # common beginning of all input file names (without year) 
#filein_pre='dailydata_' # common beginning of all input file names (without year) 
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

