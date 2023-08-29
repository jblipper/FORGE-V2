import numpy as np
from scipy import stats

# This module predicts the time evolution of human hunter-gathers. The full methods section and equations referenced in the comments on this module are outlined in the publication  'Global hunter-gatherer population densities constrained by influence of seasonality on diet composition'.

def pop(restart,nlons,nlats,ntime,eG_cst,eH_cst,c_dep,q,t2m,food_veg_ini,popu,fat,influx_fruit,decay_veg,food_veg,food_ani,fat_ani,tforage,eG,eH,tH,expend,\
estab_date,estab_elapse,estab_humsum,starv_humsum,fat_max_rec,fat_max_date,parameters,Nhum,estab_hum_a, estab_hum_b,estab_hum_x0,expend_cst,expend_veg,expend_ani,NE_ani,NE_veg,m_anabolism,m_catabolism,cdfnor_x_coe,cdfnor_sd_coe,hum_bw,popu_init,tforage_max,tforage_min,A,mort_hum,fat_max,cdfnor_x,cdfnor_sd):
  
  #==============================================================
  #====== initialization at beginning of run ====================
  valid=(t2m==t2m) # only land pixel
  
  estab_actual = np.zeros((nlats,nlons))*np.nan # creates array to keep track of births
  
  if restart=='N': # for the first loop, all state variables are set to their initial values
    popu[ valid] = popu_init
    fat[  valid] = 0.
    food_veg = food_veg_ini[0]*1
    tforage[valid]= tforage_max
    tH[valid] = tforage_max - tforage_min
    eG[valid] = eG_cst
    eH[valid] = eH_cst
    expend[valid]= expend_cst
    
    estab_date[valid]=364
    estab_elapse[valid]=0
    estab_humsum[valid] = 0
    starv_humsum[valid] = 0
    fat_max_rec[valid]  = -9999
    fat_max_date[valid] = 0
 
  # arrays to store state variables at all times of simulation are initialized 
  popu_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  fat_all     =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  expend_all      =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  estab_hum_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  starv_hum_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan 
  estab_actual_all=np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  intake_veg_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  intake_ani_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  mort_starv_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  food_veg_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  food_ani_all    =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tforage_all     =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tG_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  tH_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  eG_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  eH_all          =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  benifit_ani_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  benifit_veg_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan   
  
  estab_date_all   =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  estab_elapse_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  estab_humsum_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  starv_humsum_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  fat_max_rec_all  =np.zeros((ntime,Nhum,nlats,nlons))*np.nan
  fat_max_date_all =np.zeros((ntime,Nhum,nlats,nlons))*np.nan

  for itime in range(ntime): # loops through all days in simulation
  
      DoY= np.mod(itime,365) # determines day of year from day in simulation

      food_veg_all[    itime,0]=food_veg*1 # adds available vegetable food for day to list with all days of simulation
      fat_old=fat*1 # saves fat before the day for later comparisson with fat after the day

      tG=tforage-tH # gathering time is calculated by subtracting hunting time from foraging (hunting or gathering) time
      tforage0=tforage*1 # saves fat before the day for later comparisson with fat after the day
      if np.any((tforage<=0)+(tG<=0)): # neither foraging time nor gathering time can be less than zero
        print ('ERROR!!!!!!!!! tforage or tG<=0, itime=',itime,'MIN(tforage)=',np.nanmin(tforage),'loc=',np.nanargmin(tforage),'MIN(tG)=',np.nanmin(tG),'loc=',np.nanargmin(tG))
        sys.exit()
      if np.sum(valid)!=np.sum(tH==tH): #the number of land pixels must be the same as the number of pixels where hunting time is defined
        print ('ERROR!!!!!!!!! (1) np.sum(valid)!=np.sum(tH==tH)',np.sum(valid),np.sum(tG==tG),np.sum(tH==tH),np.sum(tforage==tforage))
        sys.exit()
  
      # update foraging time
      if itime>0: # below conditions are not calculatable for first day since they require comparisson between days
        cond_lazy = (fat>fat_max*0.5)*(delta_fat>0) # if humans have more than half of maximum fat and fat has increased from the previous day, lazy condition is true
        cond_work=np.logical_not(cond_lazy) # work condition is true when lazy condition is fales
        tforage[cond_lazy]=tforage[cond_lazy]-1 # if lazy condition is true, foraging time decreases by one hour or as much is allowed by lower bound
        tforage[cond_work]=tforage[cond_work]+1 # if work condition is true, foraging time increases by one hour or as much is allowed by upper bound 

      tforage[tforage<tforage_min]=tforage_min # foraging time cannot be below lower
      tforage[tforage>tforage_max]=tforage_max # foraging time cannot exceed upper bound

      tG=tG/tforage0*tforage #gathering time is changed to keep proportion of foraging time spent gathering correct
      tH=tH/tforage0*tforage #hunting time is changed to keep proportion of foraging time spent hunting correct

      # energy expenditure of foraging
      expend=expend_cst+ tG*expend_veg + tH*expend_ani # energy expenditure for a day is equal to the sum of energy spent on non-foraging activities (expend_cst), gathering (tG*expend_veg), and, hunting (tH*expend_ani)

      benifit_veg=A * food_veg *eG * NE_veg # energy benefit of vegetable food is calculated as the product of serching area per person hour (A), biomass density of vegetable food (food_veg), efficiency of gathering (eG) and energy density of vegetable food (NE_veg) 
      benifit_ani=A * food_ani *eH * NE_ani # energy benefit of animal food is calculated as the product of serching area per person hour (A), biomass density of animal food (food_ani), efficiency of hunting (eH) and energy density of animal food (NE_ani) 

      # calculate daily intake: note that whenever fat is about to exceed fat_max, daily intake is reduced so that fat stays at fat_max
      tmp1= tG* benifit_veg # energy intake from vegetable food
      tmp2= tH* benifit_ani # energy intake from animal food

      tmp3 = tmp1+tmp2 # total energy intake, unit: MJ/day/ind.

      cond1 =((fat+(tmp3-expend)/m_anabolism) >= fat_max) # approaching fat_max
      tmp3[cond1] =( (fat_max-fat[cond1])*m_anabolism+expend[cond1] ) # cannot exceed fat_max

      tmp4 = tmp3 * tmp1/(tmp1+tmp2) / NE_veg # dry mass of fresh vegetable food consumed
      tmp5 = tmp3 * tmp2/(tmp1+tmp2) / NE_ani # dry mass of animal food consumed
      tmp4[(tmp1+tmp2)==0]=0. # if animals did not intake food, no dry mass is consumed
      tmp5[(tmp1+tmp2)==0]=0. # if animals did not intake food, no dry mass is consumed

      intake_veg =np.where(valid,tmp4,np.nan) # daily intake of plant food
      intake_ani =np.where(valid,tmp5,np.nan) # daily intake of animal food

      # reduce tG and tH in case humans are very full (fat is about to exceed fat_max)
      tG=np.where((food_veg>0.)*(intake_veg>1.e-6),intake_veg/food_veg/A/eG,tG)
      tH=np.where((food_ani>0.)*(intake_ani>1.e-6),intake_ani/food_ani/A/eH,tH)
      tG[tG<tforage_min]=tforage_min/2. # gathering time set to half minimum foraging time
      tH[tH<tforage_min]=tforage_min/2. # hunting time set to half minimum foraging time


      # intake cannot make food density below zero
      cond2=(food_veg - intake_veg*popu < 0.) # condition true if intake would make vegetable food density below zero
      tmp= food_veg / popu # amount of vegetable food available per person if split evenly
      intake_veg[cond2]= tmp[cond2]*1 # if the calculated intake would put vegetable food density below zero, the available vegetable food is split evenly among the human population

      cond3=(food_ani - intake_ani*popu < 0.) # condition true if intake would make animal food density below zero
      tmp= food_ani / popu # amount of animal food available per person if split evenly
      intake_ani[cond3]= tmp[cond3]*1 # if the calculated intake would put animal food density below zero, the available animal food is split evenly among the human population

      if np.any(intake_veg<0): # intake of vegetable food cannot be below zero
        print ('ERROR!!!!!!!!! intake_veg<0, itime=',itime,'MIN(intake_veg)=',np.nanmin(intake_veg),'loc=',np.nanargmin(intake_veg))
        sys.exit()
      if np.any(intake_ani<0): # intake of animal food cannot be below zero
        print ('ERROR!!!!!!!!! intake_ani<0, itime=',itime,'MIN(intake_ani)=',np.nanmin(intake_ani),'loc=',np.nanargmin(intake_ani))
        sys.exit()

      intake_total=intake_veg+intake_ani  # total daily intake of animal and vegetable food in dry mass unit
      intake_total_en=intake_veg*NE_veg+intake_ani*NE_ani  # total daily intake of animal and vegetable food in energy unit

      # update tforage, and tG/tH
      tmp=np.where(intake_total>0,intake_ani/intake_total,np.nan) # calculates fraction of total dry mass intake from animal food
      meat_desire= np.exp((0.-q)*(tmp-1.)) # calclulates human craving for meat for next day based on fraction of total dry mass intake from animal food using equation 6n

      tforage= tG + tH # foraging time (sum of hunting and gathering time) must be between minimum and maximum
      if np.any((tforage<tforage_min/1.0001)+(tforage>tforage_max*1.0001)):
        print ('ERROR!!!!!!!!! tforage<tforage_min or >max, itime=',itime,'MIN(tforage)=',np.nanmin(tforage),'loc=',np.nanargmin(tforage),'MAX(tforage)=',np.nanmax(tforage),'loc=',np.nanargmax(tforage))
        sys.exit()

      tH_frac = tH/(tH+tG) # fraction of foraging time spent hunting

      tmp = benifit_ani + benifit_veg # sum of relative energy benefit of animal food and plant food (see equation 6)
      cond = (tmp>0)*(intake_total>0) # condition true if the total relative energy and total dry mass of plant and animal food are both nonzero
      tH_frac[cond] = benifit_ani[cond] / tmp[cond] * meat_desire[cond] # fraction of time spent hunting is determined based on relative benefit of hunting and meat craving (equation 6)

      tmp=np.where(intake_total>0,intake_ani/intake_total,np.nan) # calculates meat fraction of diet
      tH_frac=np.where(tmp<0.1,0.95,tH_frac)  # when meat fraction of diet is below 10%, tH_frac is set to 95%

      tH_frac[tH_frac<0.05]=0.05 # ensures hunting fraction does not go below 5%
      tH_frac[tH_frac>0.95]=0.95 # ensures hunting fraction does not exceed 95% 

      tH = tH_frac * tforage # calculates hunting time for next day from hunting fraction (of foraging time) and foraging time 
      tG = (1.-tH_frac) * tforage # calculated gathering time for next day by subtracting hunting time from foraging (sum of hunting and gathering) time

      if not np.allclose(tforage[valid],tG[valid]+tH[valid]): # foraging time must equal the sum of hunting time and gathering time
        tmp=tforage-tG-tH; i=np.nanargmin(tmp); j=np.nanargmax(tmp)
        print ('ERROR!!!!!!!!! tG+tH!=tforage, itime=',itime,np.sum(valid),np.sum(tG==tG),np.sum(tH==tH),np.sum(tforage==tforage))
        sys.exit()

      # update eG and eH (efficiencies of hunting and gathering) based on equation 2
      tmp = food_veg/(popu*50.) # calculates ratio of available biomass of vegetable food to human biomass
      eG[valid]= eG_cst # sets efficiency of gathering to maximum as default
      cond=(tmp==tmp)*(c_dep>0) # condition true if the grid cell is valid (has human population) and patch depletion rate (c_dep or c in equation 2) is nonzero
      eG[cond] = eG_cst*tmp[cond]/(tmp[cond]+c_dep) # calculates efficiency of gathering according to equation 2

      tmp = food_ani/(popu*50.) # calculates ratio of available biomass of animal food to human biomass
      eH[valid]= eH_cst # sets efficiency of hunting to maximum as default
      cond=(tmp==tmp)*(c_dep>0) # condition true if the grid cell is valid (has human population) and patch depletion rate (c_dep or c in equation 2) is nonzero
      eH[cond] = eH_cst*tmp[cond]/(tmp[cond]+c_dep) # calculates efficiency of hunting according to equation 2 

      # update plant food density
      food_veg = food_veg - intake_veg*popu + influx_fruit[itime] # vegetable food density is updated daily by subtracting consumed vegetable food and adding influx (growth) of new vegetable food
      food_veg = food_veg - decay_veg[itime]*food_veg # a proportion of vegetable food that decays each day is subtracted
      if np.any(food_veg<0): # density of vegetable food cannot be below zero
        print ('ERROR!!!!!!!!! food_veg<0, itime=',itime,'MIN(food_veg)=',np.nanmin(food_veg),'loc=',np.nanargmin(food_veg))
        sys.exit()

      # update animal food density
      food_ani = food_ani - intake_ani*popu # animal food density is updated by subtracting hunted animals; new animals born are only added on a yearly time step using the herbivore module
      if np.any(food_ani<0): # density of animal food cannot be below zero
        print ('ERROR!!!!!!!!! food_ani<0, itime=',itime,'MIN=',np.nanmin(food_ani),'loc=',np.nanargmin(food_ani))
        sys.exit()

      # update fat reserve
      m_use=np.where(intake_total_en >= expend,m_anabolism,m_catabolism) # if energy intake is less than expenditure, catabolism is used; otherwise anabolism is used
      fat = fat + (intake_total_en - expend)/m_use # intaken energy is converted to fat using metabolic rate and added to fat reserves
      if np.any(fat>1.0001*fat_max): # fat reserves cannot exceed maximum value
        print ('ERROR!!!!!!!! fat>fat_max, itime=',itime,'MAX(fat)=',np.nanmax(fat),'loc=',np.nanargmax(fat))
        sys.exit()

      # record the date with maximum body fat, on which the birth will take place
      tmp=(fat>=fat_max_rec) # fat on each day is compared to the previous maximum since last birth
      fat_max_date[tmp]=DoY # day maximum fat occurs is recorded
      fat_max_rec[tmp] = fat[tmp]*1 # amount of fat on maximum fat day is recorded

      # calculate birth and mortality rate  
      estab_hum =estab_hum_a/(1.+np.exp(-estab_hum_b*(fat/fat_max-estab_hum_x0))) # number of humans born is calculated based on a sigmoidal function, using several parameters and the fat reserves of the human population (equation 10)
  
      starv_hum = np.zeros(fat.shape)*np.nan # for each grid cell, starvation per day is taken to be the probability that fat is below zero with the mean being the fat reserves and standard deviation being 0.125 times the maximum fat
      for i in range(nlats):
          for j in range(nlons):
                  starv_hum[i,j]=stats.norm.cdf(cdfnor_x,fat[i,j],cdfnor_sd)
  
      # temporary variables to calcuate annual mean birth and mortality rate later
      estab_humsum=estab_humsum+estab_hum # keeps track of total number of births predicted by birth rate  
      starv_humsum=starv_humsum+starv_hum # keeps track of total number of starvations

      # START birth date
      e=(DoY==estab_date) # birth is assumed to occur on the day of the year when fat reserves are maximum. e is a boolean that is true when birth occurs
      if np.any(estab_elapse[e]==0): # at least one day must have passed since the last birth
        print ('ERROR!!!!!!!!!!! check estab_elapse')
        sys.exit()
      e=(DoY==estab_date)*(estab_elapse>330) # if it is the day when birth is assumed to happen and more than 330 days have elapsed since the last birth, birth occurs

      estab_actual[e]=estab_humsum[e]/estab_elapse[e] # calculates birth rate as number of births divided by days elapsed since last bith if birth occurs if birth condition is true
      expend_hum=313.*4.184/1000. # daily energy expenditure of humans
      fat[e]=fat[e]-estab_actual[e]*(expend_hum*estab_elapse[e])/m_catabolism #fat expended to support birth rate is subtracted

      tmp=(estab_actual<0) # birth rate cannot be negative
      if np.any(tmp[e]==1):
        print ('ERROR!!!!!!!!!!! estab_actual is negative')
        sys.exit()

      estab_actual_all[itime,0][e]=estab_actual[e]*1 # adds bithrate for day to the array that records birthrates on each day
      mort_starv_all[  itime,0][e]=starv_humsum[e]/estab_elapse[e] #adds starvation rate for day to the array that records starvation rates on each day
 
      # when birth occurs, population is updated using birth and death rates (equation 11)
      cond_low=(estab_actual-mort_hum-starv_humsum <= 0.)*e # true when starvation rate is higher than birth rate and birth occurs
      ncond=   (estab_actual-mort_hum-starv_humsum >  0.)*e # true when starvation rate is not higher than birth rate and birth occurs
      popu[ncond]=popu[ncond]*(1+estab_actual[ncond]-mort_hum-starv_humsum[ncond]/estab_elapse[ncond]) # if starvation rate is not higher than birth rate, there are people who die of old age, so background mortality is used to calculate new population
      popu[cond_low]=popu[cond_low]*(1+estab_actual[cond_low]-starv_humsum[cond_low]/estab_elapse[cond_low]) # if starvation rate is higher than birth rate, people do not live long enough to die of old age, so background mortality is not used to calculate new population
  
      estab_elapse[e]=0 # when birth occurs, number of days since last birth is reset to zero
      estab_humsum[e]=0 # when birth occurs, number of births predicted by birth rate since last birth is reset to zero
      starv_humsum[e]=0 # when birth occurs, number of starvations since last birth is reset to zero

      estab_date[e]=fat_max_date[e]*1 # when birth occurs, next birth date is set to date of maximum fat reserves
      fat_max_rec[e]=fat[e]*1 # when birth occurs, new maximum fat reserve since last birth is set to amount of fat at the date of the birth

      # END birth date

      # reset to initial if popu too low
      reset=(popu<popu_init) # checks if human population is too low (effectively zero) 
      fat[reset*(fat<0)]=0. # if human population is too low, fat reserves are reset to zero
      popu[reset]=popu_init # if human population is too low, human population is reset to initial value

      estab_elapse=estab_elapse+1 # time since last birth is increased by one day after each day
 
      delta_fat=fat - fat_old # calculates change in fat reserves since previous day
  
      # after each day, the new values of each state vairable is added to the array that records its value on each day of the simulation 
      popu_all[    itime,0]=popu    *1
      fat_all[     itime,0]=fat     *1  
      expend_all[      itime,0]=expend      *1 
      estab_hum_all[   itime,0]=estab_hum   *1    
      starv_hum_all[   itime,0]=starv_hum   *1 
      intake_veg_all[  itime,0]=intake_veg*1
      intake_ani_all[  itime,0]=intake_ani*1
      tforage_all[     itime,0]=tforage *1  
      tG_all[          itime,0]=tG *1  
      tH_all[          itime,0]=tH *1  
      eG_all[          itime,0]=eG *1  
      eH_all[          itime,0]=eH *1  
      food_ani_all[    itime,0]=food_ani *1
      benifit_veg_all[ itime,0]=benifit_veg *1  
      benifit_ani_all[ itime,0]=benifit_ani *1  

      estab_date_all[  itime,0]=estab_date  *1
      estab_elapse_all[itime,0]=estab_elapse*1
      estab_humsum_all[itime,0]=estab_humsum*1
      starv_humsum_all[itime,0]=starv_humsum*1
      fat_max_date_all[itime,0]=fat_max_date*1
      fat_max_rec_all[ itime,0]=fat_max_rec *1

      if itime == ntime-1:  food_veg_all[    itime,0]=food_veg*1  # this is to avoid drift in food_veg,if one-year forcing is looped

  ###################################################################
  # after all days are looped through, the arrays that contain the values of the state variables on each day are returned as output 
  return popu_all,fat_all,expend_all,intake_veg_all,estab_hum_all,starv_hum_all,estab_actual_all,mort_starv_all,food_veg_all,tforage_all,tG_all,eG_all,eH_all,intake_ani_all,food_ani_all,tH_all,benifit_veg_all,benifit_ani_all,estab_date_all,estab_elapse_all,estab_humsum_all,starv_humsum_all,fat_max_rec_all,fat_max_date_all
