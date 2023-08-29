import numpy as np
from scipy import stats

# this modules predicts the time evolution of herbivore populations given a set of initial state variables. Several model parameters and equations are used to predict these state variables after a desired time interval. For the FORGE model, this time interval is typically selected to be a year, at which point the human model is run and the results are coupled and this model is run agian. This process is repeated for the desired number of simulated years. The full methods can be found in The large mean body size of mammalian herbivores explains the productivity paradox during the last glacial maximum.

def animal(restart,nlons,nlats,ntime,Nani,t2m,bm_avail_ini,lit_avail_ini,influx_bm,decay_bm,influx_lit,decay_lit,bm_avail,lit_avail,\
estab_date,estab_elapse,estab_anisum,starv_anisum,fat_max_rec,fat_max_date,\
ani_popu,fat_ani,parameters,dt,estab_ani_a,estab_ani_b,estab_ani_x0,expend_a,expend_b,ne_AGB,ne_litter,m_anabolism,m_catabolism,delai_ugb_max,cdfnor_x_coe,cdfnor_sd_coe, t2m_lt_min, ani_bw, KK, ani_popu_init, Ldecay, intake_AGB_max, intake_litter_max, try_able_grazing, mort_ani, fat_ani_max, cdfnor_x, cdfnor_sd):
  
  #==============================================================
  #====== initialization at beginning of run ====================
  valid=(t2m==t2m) # only land pixel
  
  # variables associated with the herbivore module are created with appropriate dimensions
  delai_ugb    = np.zeros((Nani,nlats,nlons))*np.nan # depletion of fresh biomass
  able_grazing = np.zeros((Nani,nlats,nlons))*np.nan # maximum dry mass of fresh plant biomass that population of herbivores can collectively consume, unit: kgDM/day
  ugb          = np.zeros((Nani,nlats,nlons))*np.nan # determines what herbivores eat. 0 for nothing, 1 for dead biomass (litter), and 2 for fresh biomass
  estab_actual = np.zeros((Nani,nlats,nlons))*np.nan # animal birth rate
  delai_ugb[valid] = delai_ugb_max # fresh biomass is not depleted at all initially

  if restart=='N': # for the first loop, all state variables are set to their initial values
    ani_popu[ valid] = ani_popu_init
    fat_ani[  valid] = 0.
    bm_avail = bm_avail_ini*1 
    lit_avail =lit_avail_ini*1
    estab_date[valid]=364
    estab_elapse[valid]=0
    estab_anisum[valid] = 0
    starv_anisum[valid] = 0
    fat_max_rec[valid]  = -9999
    fat_max_date[valid] = 0
  
  # arrays to store state variables at all times of simulation are initialized
  ani_popu_all    =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_ani_all     =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  delai_ugb_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  ugb_all         =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  expend_all      =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_ani_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  starv_ani_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_actual_all=np.zeros((ntime,Nani,nlats,nlons))*np.nan
  intake_ind_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  intake_ind_litter_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan   
  bm_avail_all    =np.zeros((ntime,Nani,nlats,nlons))*np.nan
  lit_avail_all   =np.zeros((ntime,Nani,nlats,nlons))*np.nan
  estab_date_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_elapse_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  estab_anisum_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  starv_anisum_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_max_rec_all  =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  fat_max_date_all =np.zeros((ntime,Nani,nlats,nlons))*np.nan 
  
  for itime in range(ntime): # loops through all days in simulation

      DoY= np.mod(itime,365) # determines day of year from day in simulation

      bm_avail_all[    itime]=bm_avail*1 # adds available fresh biomass for day to list with all days of simulation
      lit_avail_all[   itime]=lit_avail*1 # adds available dead biomass for day to list with all days of simulation
  
      ### (not applicable to this version) daily energy expenditure calculated based on badyweight and temperature
      #expend=np.where(t2m<(t2m_lt_min+273.15),\
      #           expend_a/np.exp(expend_b*t2m_lt_min)*(ani_bw*1000.)**0.75/1000.,\
      #           expend_a/np.exp(expend_b*(t2m-273.15))*(ani_bw*1000.)**0.75/1000.)
      
      expend=1.3*ani_bw # daily energy expenditure of animals is assumed to be directly proportional to body weight
      #expend[np.logical_not(valid)]=np.nan

      able_grazing = ani_popu*try_able_grazing #calculates maximum dry mass of plant biomass that all herbivores can consume collectively, unit: kgDM/day 

      cond1=(bm_avail>=able_grazing) # enough fresh biomass (more than all herbivores can consume collectively)
      delai_ugb[valid*cond1]=delai_ugb[valid*cond1]+1 # if there is enough fresh biomass, it regenerates each day by 1 on a scale of -5 (completely depleted) to 0 (completely regenerated and edible)
      delai_ugb[valid*(np.logical_not(cond1))]=delai_ugb_max # if there isn't enough fresh biomass on any day, it becomes completely depleted
  
      cond2=(delai_ugb>=0) # after a sufficient time (5 days), biomass becomes completely regenerated and edible
      cond3=(lit_avail>=able_grazing) # enough dead biomass (more than all herbivores can consume collectively)
  
      # determine what to eat
      ugb[valid]=0 # by default, herbivores eat nothing (0)
      ugb[cond1*cond2]=2 # if there is enough fresh biomass and it has completely regenerated, herbivores eat fresh biomass
      ugb[(ugb!=2)*cond3]=1 # if fresh biomass has not completely regenerated but there is enough dead biomass, herbivores eat dead biomass
  
      # start the eating
      intake_ind       =np.where(valid,0.,np.nan) # animals can only intake fresh biomass on land pixels
      intake_ind_litter=np.where(valid,0.,np.nan) # animals can only intake dead biomass on land pixels
  
      cond1=(ugb==2)*(ani_popu>0) # condition true if there are animals and they consume fresh biomass
      cond2=(ugb==1)*(ani_popu>0) # condition true if there are animals and they consume dead biomass
      ncond=np.logical_not(cond1+cond2)*valid # condition true if there are neither animals that consume fresh biomass, nor animals that consume dead biomass
  
      intake_ind[cond1] = intake_AGB_max/ne_AGB
      cond3 =((fat_ani+(intake_ind*ne_AGB-expend)/m_anabolism) >= fat_ani_max) # approaching fat_ani_max
      intake_ind[cond3] =( (fat_ani_max-fat_ani[cond3])*m_anabolism+expend )/ne_AGB # cannot exceed fat_ani_max
      #intake_ind[cond3] =( (fat_ani_max-fat_ani[cond3])*m_anabolism+expend[cond3] )/ne_AGB # cannot exceed fat_ani_max
      if np.any(intake_ind<0): # animal intake of food cannot be less than zero
        print ('ERROR! intake_ind<0, itime=',itime)
        sys.exit()
  
      intake_ind_litter[cond2] = intake_litter_max/ne_litter # if the animals consume dead biomass, their energy intake is calculated as their maximum intake of biomass converted to energy using the energy density of dead biomass (litter)

      # update biomass and litter pool
      bm_avail = bm_avail - intake_ind*ani_popu + influx_bm[itime] # fresh biomass is updated daily by subtracting consumed biomass and adding influx (growth) of new biomass
      bm_avail = bm_avail - decay_bm[itime]*bm_avail # a proportion of fresh biomass that decays is subrtracted each day
      if np.any(bm_avail<0): # amount of fresh biomass cannot be below zero
        print ('ERROR!!!!!!!!! bm_avail<0, itime=',itime,'MIN(bm_avail)=',np.nanmin(bm_avail),'loc=',np.nanargmin(bm_avail))
        sys.exit()
  
      lit_avail = lit_avail - intake_ind_litter*ani_popu + influx_lit[itime] # dead biomass (litter) is updated daily by subtracting consumed biomass and adding influx (generation of litter) of new biomass
      lit_avail = lit_avail - decay_lit[itime]*lit_avail # a proportion of dead biomass (litter) that decays is subrtracted each day
      lit_avail = Ldecay*lit_avail # an additional proportion of remaining dead biomass (litter) decays each day
      if np.any(lit_avail<0): # amount of dead biomass (litter) cannot be below zero
        print ('ERROR!!!!!!!!! lit_avail<0, itime=',itime,'MIN(lit_avail)=',np.nanmin(lit_avail),'loc=',np.nanargmin(lit_avail))
        sys.exit()
  
      # update state variables
      # change in animal fat is updated by adding intake energy and subtracting expended energy, and then dividing by the metabolic rate to convert energy to fat
      fat_ani[cond1]=fat_ani[cond1] + (intake_ind[cond1]*ne_AGB-expend)/m_anabolism # change in animal fat is calculated if dead biomass is consumed, in which case anabolism is assumed to be used
      fat_ani[cond2]=fat_ani[cond2] + (intake_ind_litter[cond2]*ne_litter-expend)/m_catabolism # change in animal fat is calculated if fresh biomass is consumed, in which case catabolism is assumed to be used
      fat_ani[ncond]=fat_ani[ncond]-expend/m_catabolism # change in animal fat is calculated if no biomass is consumed, in which case catabolism is assumed to be used
      #fat_ani[cond1]=fat_ani[cond1] + (intake_ind[cond1]*ne_AGB-expend[cond1])/m_anabolism
      #fat_ani[cond2]=fat_ani[cond2] + (intake_ind_litter[cond2]*ne_litter-expend[cond2])/m_catabolism
      #fat_ani[ncond]=fat_ani[ncond]-expend[ncond]/m_catabolism

      if np.any(fat_ani>fat_ani_max): # animal fat cannot exceed maximum value
        print ('ERROR! fat_ani>fat_ani_max, itime=',itime)
        sys.exit()

      # record the maximum fat date
      tmp=(fat_ani>=fat_max_rec) # fat on each day is compared to the previous maximum since last birth
      fat_max_date[tmp]=DoY # day maximum fat occurs is recorded
      fat_max_rec[tmp] = fat_ani[tmp]*1 # amount of fat on maximum fat day is recorded
 
      # calculate birth rate and mortality rate 
      estab_ani =estab_ani_a/(1.+np.exp(-estab_ani_b*(fat_ani/fat_ani_max-estab_ani_x0))) # number animals born is calculated based on a sigmoidal function, using several parameters and the fat reserves of the animal population
      estab_ani[fat_ani<=0]=0. # if fat reserves are not greater than zero, no animals are born
  
      starv_ani =np.where(fat_ani>=0,0.,np.nan) # for each grid cell, starvation per day is taken to be the probability that fat is below zero with the mean being the fat reserves and standard deviation being 0.125 times the maximum fat
      for k in range(Nani):
        for i in range(nlats):
            for j in range(nlons):
                if fat_ani[k,i,j]<0:
                    starv_ani[k,i,j]=stats.norm.cdf(cdfnor_x,fat_ani[k,i,j],cdfnor_sd)
      
      # temporary variables to calcuate annual mean birth and mortality rate later
      estab_anisum=estab_anisum+estab_ani # keeps track of total number of births predicted by birth rate    
      starv_anisum=starv_anisum+starv_ani # keeps track of total number of starvations
  
      # START birth date
      e=(DoY==estab_date) # birth is assumed to occur on the day of the year when fat reserves are maximum. e is a boolean that is true when birth occurs
      if np.any(estab_elapse[e]==0): # at least one day must have passed since the last birth
        print ('ERROR!!!!!!!!!!! check estab_elapse')
        sys.exit()
      e=(DoY==estab_date)*(estab_elapse>330) # if it is the day when birth is assumed to happen and more than 330 days have elapsed since the last birth, birth occurs

      cond1=(fat_ani*m_catabolism/expend > estab_anisum )*e # condition true if there is enough energy from fat to support total number of births predicted by the birth rate and birth condition (e) is true 
      cond2=(fat_ani*m_catabolism/expend <=estab_anisum )*e # condition true if there is not enough energy from fat to support total number of births predicted by the birth rate and birth condition (e) is true

      estab_actual[cond1]=estab_anisum[cond1]/estab_elapse[cond1] # calculates birth rate as number of births possible divided by days elapsed since last birth if there is enough energy from fat to support all new animals and birth condition is true
      estab_actual[cond2]=fat_ani[cond2]*m_catabolism/(expend*estab_elapse[cond2]) # calculates birth rate as number of births possible with given fat reserves divided by days elapsed since last birth if there is not enough energy from fat to support all new animals and birth condition is true
      #estab_actual[cond2]=fat_ani[cond2]*m_catabolism/(expend[cond2]*estab_elapse[cond2])
      estab_actual[estab_actual<0]=0. # birth rate cannot be below zero
  
      #fat_ani[e]=fat_ani[e]-estab_actual[e]*(expend[e]*estab_elapse[e])/m_catabolism
      fat_ani[e]=fat_ani[e]-estab_actual[e]*(expend*estab_elapse[e])/m_catabolism # fat expended to support birth rate is subtracted

      ani_popu[e]=ani_popu[e]*(1+estab_actual[e]-mort_ani-starv_anisum[e]/estab_elapse[e])-(estab_ani_a-mort_ani)/KK*ani_popu[e]*ani_popu[e]*1.e6 # animal population is updated based of birth and death rates

      estab_elapse[e]=0 # when birth occurs, number of days since last birth is reset to zero
      estab_anisum[e]=0 # when birth occurs, number of births predicted by birth rate since last birth is reset to zero
      starv_anisum[e]=0 # when birth occurs, number of starvations since last birth is reset to zero

      estab_date[e]=fat_max_date[e]*1 # when birth occurs, next birth date is set to date of maximum fat reserves
      fat_max_rec[e]=fat_ani[e]*1 # when birth occurs, new maximum fat reserve since last birth is set to amount of fat at the date of the birth

      # END birth date

      # reset to initial if popu too low
      reset=(ani_popu<ani_popu_init) #checks if animal population is too low (effectively zero) 
      fat_ani[reset]=0. # if animal population is too low, fat reserves are reset to zero
      ani_popu[reset]=ani_popu_init # if animal population is too low, animal population is reset to initial value

      estab_elapse=estab_elapse+1 # time since last birth is increased by one day after each day 
  
      # after each day, the new values of each state vairable is added to the array that records its value on each day of the simulation 
      ani_popu_all[    itime]=ani_popu    *1
      fat_ani_all[     itime]=fat_ani     *1  
      delai_ugb_all[   itime]=delai_ugb   *1 
      ugb_all[         itime]=ugb         *1 
      expend_all[      itime]=expend      *1 
      estab_ani_all[   itime]=estab_ani   *1    
      starv_ani_all[   itime]=starv_ani   *1 
      estab_actual_all[itime]=estab_actual*1
      intake_ind_all[  itime]=intake_ind*1
      intake_ind_litter_all[itime]=intake_ind_litter*1
      estab_date_all[  itime]=estab_date  *1
      estab_elapse_all[itime]=estab_elapse*1
      estab_anisum_all[itime]=estab_anisum*1
      starv_anisum_all[itime]=starv_anisum*1
      fat_max_date_all[itime]=fat_max_date*1
      fat_max_rec_all[ itime]=fat_max_rec *1

  ###################################################################
  # after all days are looped through, the arrays that contain the values of the state variables on each day are returned as output
  return ani_popu_all,fat_ani_all,delai_ugb_all,ugb_all,expend_all,intake_ind_all,intake_ind_litter_all,estab_ani_all,starv_ani_all,estab_actual_all,bm_avail_all,lit_avail_all,estab_date_all,estab_elapse_all,estab_anisum_all,starv_anisum_all,fat_max_rec_all,fat_max_date_all
