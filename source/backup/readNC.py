import numpy as np
from netCDF4 import Dataset


# this module reads the input file, obtaining the dimensions (lattitude, longidude, and time) of the input data as well as values of the data variables, gross primary production (GPP), annual precipitation (prec), and two-metre temperature (t2m) at each set of coordinates

def read_var(path,filename,varlist,varuse):

    adict=locals() # creates dictionary of path, filename, varlist, and varuse to their values 
    nc=Dataset(path+filename,'r')  # opens input file
    nlons=nc.dimensions[varlist[0]].size # obtains number of longitudes in input file
    nlats=nc.dimensions[varlist[1]].size # obtains number of lattitudes in input file
    #nlevs=nc.dimensions[varlist[2]].size
    nlevs=1 # number of functional types
    ntime=nc.dimensions[varlist[2]].size # obtains number of times in input file; Note: distinct from time in output file 

    for i in range(3,len(varlist)): # loops through data variables,   
        adict[varuse[i]]=nc.variables[varlist[i]][:] # creates a dictionary that assigns the name of each data variable (GPP, prec, and t2m) to an array that contains its values at all coordinates
        adict[varuse[i]]=adict[varuse[i]].filled(np.nan) # for each data variable, all missing values in the array are filled with nan (not a numbe)
        adict[varuse[i]][adict[varuse[i]]<-8e+20]=np.nan # NOTE erronious very large negative values are replaced with nan 
    nc.close()

    varin={varuse[3]: adict[varuse[3]]} # creates a dictionary to store all the data variables, assigning all variable names to arrays that contains their values at all coordinates
    for i in range(4,len(varlist)): # adds all variables to the dictionary
        varin.update({varuse[i]: adict[varuse[i]]}) #update the dictionary
    
    return nlons,nlats,nlevs,ntime,varin # returns dimensions of the input data and the dictionary storing all data variables
