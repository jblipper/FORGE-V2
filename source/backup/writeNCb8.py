import numpy as np
from netCDF4 import Dataset

# This module saves the output of the FORGE model as a netCDF file in the location specified in the 'run.py' module.

def write_var(res,pathin,pathout,filein,fileout,nlat,nlon,ntime,Nveg,Nani,Nhum,outlist,varout,shortname):
  
  ncref=Dataset(pathin+filein,'r') # opens input netCDF file
  nc=Dataset(pathout+fileout,'w',format='NETCDF3_CLASSIC') # createss output netCDF file
  
  nc.createDimension('lat',nlat) # creates lattitude dimension in output file with same size as input file
  nc.createDimension('lon',nlon) # creates longitude dimension in output file with same size as input file
  nc.createDimension('veg',Nveg) # creates plant functional type dimension in output file with same size as input file
  nc.createDimension('ani',Nani) # creates animal functional type dimension in output file with same size as input file
  nc.createDimension('hum',Nhum) # creates human functional type dimension in output file with same size as input file
  nc.createDimension('time',None) # creates time dimension with size of zero
  
  var=nc.createVariable('lat','f',('lat',),zlib=True) # creates new lattitude netCDF variable in output file
  var[:]=ncref.variables['lat'][:] # copies lattitudes from the input file to the output file
  vattr=['standard_name','long_name','units','axis'] # list of names of attributes for the lattitude variable
  for a in vattr: # retrieves each attribute of the lattitude variable from input file and assigns it to the lattitude variable in the output file
    attr=ncref.variables['lat'].__getattribute__(a)
    nc.variables['lat'].__setattr__(a,attr)
  
  var=nc.createVariable('lon','f',('lon',),zlib=True) # creates new longitude netCDF variable in output file
  var[:]=ncref.variables['lon'][:] # copies longitudes from the input file to the output file
  vattr=['standard_name','long_name','units','axis'] # list of names of attributes for the longitude  variable
  for a in vattr: # retrieves each attribute of the longitude variable from input file and assigns it to the longitude variable in the output file
    attr=ncref.variables['lon'].__getattribute__(a)
    nc.variables['lon'].__setattr__(a,attr)
  
  var=nc.createVariable('time','f',('time',),zlib=True) # creates time lattitude netCDF variable in output file
  # time axis is created in units of seconds since 1901-01-01 00:00:00 with the specified resolution of daliy, month, or year
  if res=='daily': var[:]=np.arange(ntime)*86400 # converts number of days since 1901-01-01 00:00:00 to seconds and assigns to time variable
  if res=='month': # converts the number of months since 1901-01 to seconds, choosing approximately the middle of the month and assigns to time variable
      tmp=np.cumsum(calendar.mdays)
      var[:]=np.array([ (i/12)*86400*365+(tmp[np.mod(i,12)]+15)*86400 for i in range(ntime)])
  if res=='year' : var[:]=np.arange(ntime)*86400*365+86400*182 # converts years since 1901 to seconds, choosing the middle of the year and assigns to time variable
  nc.variables['time'].__setattr__('units',u'seconds since 1901-01-01 00:00:00') # sets unit attribute of time variable to seconds since 1901-01-01 00:00:00
  nc.variables['time'].__setattr__('calendar',u'365_day') # sets calendar attribute of time variable to 365_day
  nc.variables['time'].__setattr__('long_name',u'Time axis') # sets long_name attribute of time variable to Time axis
  nc.variables['time'].__setattr__('standard_name',u'time') # sets standard_name attribute of time variable to time
  nc.variables['time'].__setattr__('axis',u'T') # sets axis attribute of time variable to T
  
  #####################################################################
  for loop in range(len(outlist)): #loops through list of state variables output from simulation
      tmp=varout[outlist[loop][0]]*1 # saves state variable (data array) as temporary variable 
      tmp2=outlist[loop][1] # saves type (ani or hum) of state variable as temporary variable
      name=outlist[loop][0] # saves long name of state variable as name
      short=shortname[loop] # saves short name of variable as short
      # creates array of appropriate dimensions for each variable based on its type
      if tmp2=='2':    var=nc.createVariable(name,'f',('lat','lon'),zlib=True)
      if tmp2=='3':    var=nc.createVariable(name,'f',('time','lat','lon'),zlib=True)
      if tmp2=='veg':  var=nc.createVariable(name,'f',('time','veg','lat','lon'),zlib=True)
      if tmp2=='ani':  var=nc.createVariable(name,'f',('time','ani','lat','lon'),zlib=True)
      if tmp2=='hum':  var=nc.createVariable(name,'f',('time','hum','lat','lon'),zlib=True)

      tmp[tmp!=tmp]=-9.e20 # creates placeholder value for missing values
      var[:]=tmp # saves state variable (data array) to nedCDF file 

      nc.variables[name].__setattr__('missing_value',-9.e20) # sets 'missing value' attribute to coordinates with missing values for the state variable
      nc.variables[name].setncattr('short_name',short) # adds 'short name' attribute to the state variable 
  #####################################################################
  ncref.close() # closes input netCDF file
  nc.close() # closes output netCDF file

  return #exits function
