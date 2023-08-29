import numpy as np


def time_aggre(res,varin): # takes as arguments desired time resolution and a data variable to have its time resolution changed

    nds=varin.shape # obtains dimensions of variable to be reshaped
    nyr=int(nds[0]/365) # selects size of time dimension (number of days to be simulated) and converts to years

    if res=='daily': varout = varin * 1 # if the desired resolution is daily, time resolution is unchanged
    
    if res=='month': # changes time resolution from daily to month 
        a=varin * 1 # saves copy of input variable as temporary variable 'a'
        a.shape=(nyr,365,nds[1],nds[2],nds[3]) # reshapes data to allow for calculation of monthly mean based on number of days in month
        b=np.zeros((nyr,12,nds[1],nds[2],nds[3]))*np.nan # creates new variable to hold data with monthly time resolution
        ind=np.cumsum(calendar.mdays) # keeps track of number of days since start of year at each month to deterime bounds between which to take monthly average 
        for imon in range(12): # loops through each month of year  
            b[:,imon]=np.nanmean(a[:,ind[imon]:ind[imon+1]],axis=1) # calculates mean value of the data for the month
        varout=b.reshape(nyr*12,nds[1],nds[2],nds[3]) # reshapes the monthly averaged data to a monthly time step
    
    if res=='year': # changes time resolution from daily to year
        a=varin * 1 # saves copy of input variable as temporary variable 'a'
        a.shape=(nyr,365,nds[1],nds[2],nds[3]) # reshapes data to allow for calculation of yearly mean
        varout=np.nanmean(a,axis=1) # calculates mean value of data for the year

    return varout # returns data variable with updated time resolution
