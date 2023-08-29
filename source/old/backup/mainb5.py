import numpy as np
import readNC
import writeNC
import herbivore 
import human
import functions as Fun

#Overview: The FORGE model simulates the time evolution of interacting hunter-gathers and animals tracking several important quantities relating to the human and animal populations, the full list of which can be found on line 30. In this version, observation-based data in one degree of latitude by one degree of longitude resolution from the SESAME database is used as an input, including two metre temperature, yearly precipitation rates, and Gross Primary Production sampled on an eight day time step. The GPP data is converted to a daily timestep using linear interpolation. 

def run(pathin,pathout,res,run_year,phi_v,phi_a,eG_cst,eH_cst,c_dep,q,filein_pre,fileout,parameters):
  
  ### variable name in the input file
  varlist=['lon','lat','time',\
          't2m','prec','GPP'] #the model takes three data variables from the input file, two-metre temperature (t2m) K or C, mean annual precipitation (prec) in mm/year, and Gross Primary Production (GPP) in gC/m^2/y (grams of carbon per square metre per year)
### variable name used in the FORGE model,corresponding to the names from the input
  varuse =['', '',  '',           \
          't2m','prec','GPP']  #variables from input are renamed for use in the FORGE model

  
#Model Parameters
  Nani=1 # number of animal functional types
  Nhum=1 # number of human funtional types
  CtoDM=0.45  # conversion factor between dry mass and carbon mass
  edi_frac=0.3 # proportion of NPP accessible to animals
  ani_edi=0.28 # conversion between live weight and dry mass after excluding water and bones for the animal food
  decay_frac=0.5 # proportion of fresh biomass available to herbivores which decays each day
  lit_prop=0.5 # ratio of amount of dead biomass (litter) generated per day to NPP
  decay_lit_frac=0.3 # proportion of dead biomass (litter) that decays per day 
  gpp_to_npp=0.7 #proportion of GPP that is not used for plant respiration and contributes to NPP
  phinpp = 0.7  # proportion of npp available to animals by biome  
  decay_frac_veg_food_hum = 0.25 # proportion of of plant food edible to humans that decays each day 
 



  ### variable list to be outputted
  adict=locals()
  
  

  outlist=[\
['Animal Population (ind.)','ani'],\
['Per Capita Energy Reserve from Fat in Animals (kg*ind.^-1)','ani'],\
['q1','ani'],\
['q2','ani'],\
['Per Capita Daily Energy Expenditure of Animals (MJ*day^-1*ind.^-1)','ani'],\
['Dry Mass Intake Rate for Fresh Biomass(kgDM*day^-1*ind.^-1)','ani'],\
['Dry Mass Intake Rate for Dead Biomass(kgDM*day^-1*ind.^-1)','ani'],\
['Animal Birth Rate (births*year^-1)','ani'],\
['Number of Animals Starved in Day (ind.)', 'ani'],\
['Animal Birh Rate Averaged Over Year (births*year^-1)','ani'],\
['Available Fresh Biomass(kgDM*m^-2)','ani'],\
['Available Dead Biomass (kgDM*m^-2)','ani'],\
['q3',  'ani'],\
['q4', 'ani'],\
['Influx of Available Fresh Biomass (kgDM*m^-2*day^-1)', 'ani'],\
['Influx of Available Dead Biomass (kgDM^*m^-2*day^-1)','ani'],\
['Human Population (ind.)','hum'],\
['Per Capita Energy Reserve from Fat in Humans (kg*ind.^-1)','hum'],\
['Per Capita Daily Energy Expenditure of Humans (MJ*day^-1*ind.^-1)','hum'],\
['Daily Intake of Plant Food by Humans (kgDM*ind.^-1*day^-1)','hum'],\
['Human Birth Rate (births*year^-1)','hum'],\
['Number of Humans Starved in Day (ind.)','hum'],\
['Human Birh Rate Averaged Over Year (births*year^-1)','hum'],\
['q5','hum'],\
['Biomass Density of Edible Vegitation (kgDM*m^-2)','hum'],\
['Time Spent Foraging (hours*day^-1)','hum'],\
['Time allocated to Gathering (hours*day^-1)','hum'],\
['Efficiency of Gathering (fraction of biomass aquired)','hum'],\
['Efficiency of Hunting (fraction of biomass aquired)','hum'],\
['Daily Dry Matter Intake of Meat (kgDM*day^-1*ind^-1)','hum'],\
['Biomass Density of Edible Animals (kgDM*m^-2)','hum'],\
['Time Allocated to Hunting (hours*day^-1)','hum'],\
['Relative Energy Benefit of Gathering (kgDM*m^-2)','hum'],\
['Relative Energy Benefit of Hunting (kgDM*m^-2)','hum'],\
['q6','hum'],\
['q7','hum'],\
['Influx of Plant Food for Humans (kgDm*m^-2*day^-1)','hum']] #list of variables to be tracked

  shortnames=['ani_popu',\
'ani_energy',\
'delai_ugb',\
'ugb_ani',\
'expend_ani',\
'intake_ind_ani',\
'intake_ind_litter_ani',\
'estab_ani',\
'starv_ani',\
'estab_actual_ani',\
'bm_avail',\
'lit_avail',\
'est_date',\
'est_elapse',\
'in_bm_avail',\
'in_lit_avail',\
'hum_popu',\
'fat_hum',\
'expend_hum',\
'intake_veg',\
'estab_hum',\
'starv_hum',\
'estab_actual_hum',\
'mort_starv_hum',\
'food_veg',\
'tforage',\
'tG',\
'eG',\
'eH',\
'intake_ani',\
'food_ani',\
'tH',\
'benifit_veg',\
'benifit_ani',\
'est_date_hum',\
'est_elapse_hum',\
'influx_v'] #list of variables to be tracked in short form
  
  DayInYr = 365 # days in a year
  ntime_run = 365*run_year # run_year is number of years the simulation is run and ntime run is the number of days the simulation is run
  
  for loop in range(run_year):
  
    print (loop)
    #==============================================================
    #=== read inputs from netcdf files ============================
    #==============================================================
    ninput_years=3 # number of files from the INPUT folder used. Each file corresponds to a year of GPP data starting in 2001. In this version, there are 19 files available and the same prec and t2m data is used for all files.
    input_year = 2001+np.mod(loop,ninput_years)  # the input files are cycled
    filein = filein_pre+str(input_year)+'.nc' #names of files in input folder
  
    nlons,nlats,nlevs,ntime,varin = readNC.read_var(pathin,filein,varlist,varuse) #input dimensions and data is obtained from the data files. nlons is the number of longitudes (360 for one degree resolution), nlats is the number of latitudes (180 for one degree resolution), nlevs is the number of plant or animal functional types (assumed to be one for this version), and ntime is the number of days in the year (365), and varin being the input variables(two meter temperature, precipitation, and gross primary production.
    if loop==0: print ('nlats=',nlats,'nlons=',nlons,'ntime=',ntime,'ntime_run=',ntime_run) #displays dimensions at beginning of running
    
    if np.mod(ntime_run,ntime)!=0 or DayInYr!=ntime: # years must be assumed to have 365 days so leap years must be trimmed
        print ('WRONG in ntime_run!')
        sys.exit()
    
    #==============================================================
    #=== prepare inputs for the herbiovore module =================
    #==============================================================
    
    #(not applicable for this version) # calculates grass and tree fractional covers
    #fgrass= varin['land_cov'][:,0]*1  # grass fractional cover
    #ftree = varin['land_cov'][:,1]*1  # tree fractional cover
    #fgrass_yr= np.mean(fgrass,axis=0) # annual mean
    #ftree_yr = np.mean(ftree, axis=0) # annual mean
  
    t2m_use=np.repeat(varin['t2m'].reshape(1,nlats*nlons),[Nani]*1,axis=0)  # to add the HFT dimension 
    t2m_use.shape=(Nani,nlats,nlons) 

    bm_avail_ini_Nani = np.zeros((Nani,nlats,nlons))  # initial food (plant living biomass) density for herbivores (set to zero)
    lit_avail_ini_Nani = np.zeros((Nani,nlats,nlons)) # initial food (plant litter) density for herbivores (set to zero)
    
    


    # The inputs for the herbivore model and some of the inputs for the human model, including Net primary production (NPP), daily influx of fresh and dead biomass for herbivores, and daily decay rate of fresh and dead biomass for herbivores are caculated using the input data for Gross Primary Production (GPP), precipitation (prec), and two meter temperature (t2m), and the model parameters.  

    
    phinpp= phinpp*np.ones((1,nlats,nlons)) #replace this array with biome-specific fraction of NPP available if land cover data is added as an input

    npp= varin['GPP']*gpp_to_npp/(1000. * CtoDM) #converts GPP in units of gC/m^2/day to NPP in units of kgDM/m^2/day
    
    #(not applicable in this version) npptotal=npp[:,0,:,:]+npp[:,1,:,:]+npp[:,2,:,:]  fractional covers of plant functional types PFTs are summed to obtain the total npp
    
    npptotal=npp #in this version, fractional covers are not considered so summation of fractional covers is not necessary

    influx_bm_avail_Nani=np.zeros((365,Nani,nlats,nlons)) #initializes array for daily influx of fresh biomass available to herbivoresi, unit: kgDM/m2/day

    influx_bm_avail_Nani[:,0,:,:]=npptotal[:,:,:]*phinpp #assumes that the daily influx of fresh biomnass available to herbivores is a constant fraction of npp 

    #(not applicable for this version) obtains daily influx of fresh biomass from variable in input files not present in this version
    #influx_bm_avail_Nani= varin['in_bm_ani']/(1000. * CtoDM) # influx to the edible litter for herbivores, unit: kgDM/m2/day    
    #influx_bm_avail_Nani[:,1,:,:] = influx_bm_avail_Nani[:,1,:,:] *edi_frac_tree   for browsers, only 10% is reachable.

    decay_bm_avail_Nani= np.ones(ntime_run) *decay_frac  # daily decay rate of the edible biomass for herbivores
  
    influx_lit_avail_Nani=np.zeros((365,Nani,nlats,nlons)) #initializes array for daily influx of available dead biomass (litter) available to herbivores

    influx_lit_avail_Nani[:,0,:,:]= npptotal[:,:,:]*lit_prop # influx to the edible litter for herbivores, unit: kgDM/m2/day. This value is assumed to be proportional to npp.

    #(not applicable for this version) obtains daily influx of dead biomass (litter) from variable in input files not present in this version
    #influx_lit_avail_Nani= varin['in_lit_ani']/(1000. * CtoDM) # influx to the edible litter for herbivores, unit: kgDM/m2/day
    # influx_lit_avail_Nani[:,1,:,:] = influx_lit_avail_Nani[:,1,:,:] *edi_frac_tree    

    decay_lit_avail_Nani= np.ones(ntime_run) *decay_lit_frac # daily decay rate of the edible litter for herbivores

    ### do some check ###
    # influx of fresh biomass cannot be less that zero
    if np.any(influx_bm_avail_Nani<0):
      loc=np.nanargmin(influx_bm_avail_Nani)
      print ('influx_bm_avail_Nani<0 MIN=',np.nanmin(influx_bm_avail_Nani),'loc=',loc)
      sys.exit()
    # influx of dead biomass cannot be less than zero
    if np.any(influx_lit_avail_Nani<0):
      loc=np.nanargmin(influx_lit_avail_Nani)
      print ('influx_lit_avail_Nani<0 MIN=',np.nanmin(influx_lit_avail_Nani),'loc=',loc)
      sys.exit()
    # proportion of fresh biomass that decays per day must be between zero and one
    if np.any((decay_bm_avail_Nani<0)+(decay_bm_avail_Nani>1)):
      print ('decay_bm_avail_Nani<0 or >1')
      sys.exit()
    # proportion of dead biomass that decays per day must be between zero and one
    if np.any((decay_lit_avail_Nani<0)+(decay_lit_avail_Nani>1)):
      print ('decay_lit_avail_Nani<0 or >1')
      sys.exit()

    #==============================================================
    #=== calculate inputs for the human module ====================
    #==============================================================
    #These inputs are relating to human consumption of plant food and do not include animal food, which is provided by the output of the herbivore module.

    food_veg_ini=np.zeros((nlats,nlons)) #initializes array for influx of plant food for humans in each grid cell
 
    ### calculate the influx of plant food for humans, unit: kgDM/m2/day. see Equation(7)
    tmp=varin['prec'][:,:]*1 # saves precipitation as temporary variable
    fBNPP=(88.3-0.0534*tmp)/100. ; fBNPP[fBNPP<0.2]=0.2 # calculates the fraction of grass NPP allocated to below ground (fBNPP) from annual precipitation    
    influx_veg = (npptotal*0.065+npptotal*0.5*fBNPP) * phi_v # calculates influx of plant food for humans from fBNPP, NPP, and the accesible fraction of vegetable biomass to humans (phi_v)
    #influx_veg = ( (npp[:,0]*fgrass + npp[:,1]*ftree)*0.065 + npp[:,0]*fgrass*fBNPP ) * phi_v

    decay_veg=np.ones(ntime_run)*decay_frac_veg_food_hum # assumes daily decay rate of the edible plants for humans is constant
    #decay_veg= varin['turn_bm_hum'][:,0,:,:]*1  # daily decay rate of the edible plants for humans
  
    ### do some check ###
    #influx of plant food for humans cannot be less than zero
    if np.any(influx_veg<0):
      loc=np.nanargmin(influx_veg)
      print ('influx_veg<0 MIN=',np.nanmin(influx_veg),'loc=',loc)
      sys.exit()
    #daily decay rate of edible plants must be between zero and one
    if np.any((decay_veg<0)+(decay_veg>1)):
      print ('decay_veg<0 or >1')
      sys.exit()

    #==============================================================
    #=== run the herbivore and human model ========================
    #==============================================================
    if loop==0: # at the beginning of simulation, do the initialization 
    # for the first year, simulated, all state variables are initialized to initial values specified in human and herbivore modules. For all other year, state variables have their values after the previous year as initial values
      #for the initialization, all state variables must be created as arrays with appropriate shapes
      ani_popu_0= np.zeros((Nani,nlats,nlons))*np.nan
      fat_ani_0 = np.zeros((Nani,nlats,nlons))*np.nan
      bm_avail_0= np.zeros((Nani,nlats,nlons))*np.nan
      lit_avail_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_date_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_elapse_0=np.zeros((Nani,nlats,nlons))*np.nan
      estab_anisum_0=np.zeros((Nani,nlats,nlons))*np.nan
      starv_anisum_0=np.zeros((Nani,nlats,nlons))*np.nan
      fat_max_rec_0=np.zeros((Nani,nlats,nlons))*np.nan
      fat_max_date_0=np.zeros((Nani,nlats,nlons))*np.nan
      
      hum_popu_0= np.zeros((nlats,nlons))*np.nan
      fat_hum_0 = np.zeros((nlats,nlons))*np.nan
      food_veg_0= np.zeros((nlats,nlons))*np.nan
      tforage_0 = np.zeros((nlats,nlons))*np.nan
      eG_0      = np.zeros((nlats,nlons))*np.nan
      eH_0      = np.zeros((nlats,nlons))*np.nan
      tH_0      = np.zeros((nlats,nlons))*np.nan
      expend_0  = np.zeros((nlats,nlons))*np.nan
      estab_date_hum_0  =np.zeros((nlats,nlons))*np.nan
      estab_elapse_hum_0=np.zeros((nlats,nlons))*np.nan
      estab_humsum_0    =np.zeros((nlats,nlons))*np.nan
      starv_humsum_0    =np.zeros((nlats,nlons))*np.nan
      fat_max_rec_hum_0 =np.zeros((nlats,nlons))*np.nan
      fat_max_date_hum_0=np.zeros((nlats,nlons))*np.nan
  
      restart = 'N' # causes state variables to be initialized to initial values in herbivore and human modules in for the first year
  
    else:
      restart = 'Y'  # causes state variables not to be initialized, meaning their initial value is their value from the previous year for years after the first year

    ### run the herbivore module
    Vs = herbivore.animal(restart,nlons,nlats,ntime,Nani,t2m_use,bm_avail_ini_Nani,lit_avail_ini_Nani,influx_bm_avail_Nani,decay_bm_avail_Nani,influx_lit_avail_Nani,decay_lit_avail_Nani,bm_avail_0,lit_avail_0,\
         estab_date_0,estab_elapse_0,estab_anisum_0,starv_anisum_0,fat_max_rec_0,fat_max_date_0,\
         ani_popu_0,fat_ani_0,parameters)
    ### record the state variables
    # the output of the herbivore module updates state variables, which are imputted as initial values for the human module as well as the herbivore module the following year
    ani_popu_0=Vs[0][-1]*1
    fat_ani_0 =Vs[1][-1]*1
    bm_avail_0=Vs[10][-1]*1
    lit_avail_0=Vs[11][-1]*1
    estab_date_0   = Vs[12][-1]*1
    estab_elapse_0 = Vs[13][-1]*1
    estab_anisum_0 = Vs[14][-1]*1
    starv_anisum_0 = Vs[15][-1]*1
    fat_max_rec_0  = Vs[16][-1]*1
    fat_max_date_0 = Vs[17][-1]*1
   
    ### prepare the output variables
    s=influx_bm_avail_Nani.shape ; s=(1,s[0],s[1],s[2],s[3])
    Vs_out=np.concatenate((Vs[:14],influx_bm_avail_Nani.reshape(s),influx_lit_avail_Nani.reshape(s)),axis=0) # write out the food influx for herbivores as well

    ### calculates density of animal food available for humans as input for human module
    ani_tmp1 = ani_popu_0[0] * 180.*ani_edi *phi_a ### convert from kg(live weight)/m2 to edible/accessible kgDM/m2
    # (not applicable to this version) ani_tmp1= ani_popu_0[0] * fgrass_yr*180.*ani_edi * phi_a
    # (not applicable to this version) ani_tmp2= ani_popu_0[1] * ftree_yr *180.*ani_edi * phi_a
    ani_DM_0 = ani_tmp1 # (not applicable to this version) sum up the density on tree and grass fractional covers to derive the density for the whole grid cell
    ani_old = ani_DM_0 *1 # saves animal density before human module is run to compare with animal density after
   
    ### run the human module
    Ps  = human.pop(restart,nlons,nlats,ntime,eG_cst,eH_cst,c_dep,q,varin['t2m'],food_veg_ini,hum_popu_0,fat_hum_0,influx_veg,decay_veg,food_veg_0,ani_DM_0,fat_ani_0,tforage_0,eG_0,eH_0,tH_0,expend_0,\
  estab_date_hum_0,estab_elapse_hum_0,estab_humsum_0,starv_humsum_0,fat_max_rec_hum_0,fat_max_date_hum_0,parameters)
    
    ### record the state variables
    # the output of the herbivore module updates state variables, which are imputted as initial values for the human module the following year
    hum_popu_0= Ps[0][-1,0]*1
    fat_hum_0 = Ps[1][-1,0]*1
    food_veg_0= Ps[8][-1,0]*1
    tforage_0 = Ps[9][-1,0]*1
    eG_0      =Ps[11][-1,0]*1
    eH_0      =Ps[12][-1,0]*1
    tH_0      =Ps[15][-1,0]*1
    ani_DM_0  =Ps[14][-1,0]*1
    estab_date_hum_0  = Ps[18][-1,0]*1
    estab_elapse_hum_0= Ps[19][-1,0]*1
    estab_humsum_0    = Ps[20][-1,0]*1
    starv_humsum_0    = Ps[21][-1,0]*1
    fat_max_rec_hum_0 = Ps[22][-1,0]*1
    fat_max_date_hum_0= Ps[23][-1,0]*1
  
    ### prepare the output variables
    s=influx_veg.shape ; s=(1,s[0],1,s[1],s[2])
    Ps_out=np.concatenate((Ps[:20],influx_veg.reshape(s)),axis=0) # write out the influx to plant food for humans as well
     
    ### (not applicable to this version) to subtract the hunted herbivores: assume a proportional reduction between grazers and browsers
    ani_new = ani_DM_0 *1 # new animal density is recorded as animal density after the human module is run (hunting)
    if np.any(ani_new>ani_old): # animal density cannot increase after hunting
      print ('ERROR!!!!!!!!! animal increases after hunting ')
      sys.exit()
    ani_dif1=np.where(ani_old==0,0,ani_tmp1*(ani_old-ani_new)/ani_old) # (not applicable in this version) splits up number of hunted animals between functional types
    #ani_dif2=np.where(ani_old==0,0,ani_tmp2*(ani_old-ani_new)/ani_old)
 
    ### hunted animals are subtracted from animal population before new animal population is used as input for the next year
    ani_popu_0[0]=ani_popu_0[0] - np.where(1>0, ani_dif1/(180.*ani_edi), 0.) # note: the unit of animal population density in the herbivore module is per unit PFT (grass or tree) area, instead of per ground area, so here it needs to divide by the fractional cover of grass or tree before being used in the herbivore module
   # ani_popu_0[1]=ani_popu_0[1] - np.where(ftree_yr >0, ani_dif2/(ftree_yr *180.*ani_edi), 0.)
    if np.any(ani_popu_0<0): # number of animals hunted cannot be greater than animal population
      print ('ERROR!!!!!!!!! animal less than 0')
      sys.exit()
    
    ### prepare the output list: average along the time dimension if res='month' or 'year'
    Numv=16
    if loop==0:
      for i in range(len(outlist)):
        if i<Numv:  
            adict[outlist[i][0] ]=Fun.time_aggre(res,Vs_out[i] )
        if i>=Numv: 
            adict[outlist[i][0] ]=Fun.time_aggre(res,Ps_out[i-Numv] )
    else:
      for i in range(len(outlist)):
        if i<Numv:  tmp=Fun.time_aggre(res,Vs_out[i] )
        if i>=Numv: tmp=Fun.time_aggre(res,Ps_out[i-Numv] )
        adict[outlist[i][0]] =np.concatenate((adict[outlist[i][0]],tmp),axis=0)
    
  
  #==============================================================
  #=== write outputs into netcdf files ==========================
  #==============================================================
  varout={outlist[0][0]: adict[outlist[0][0]]} # a dictionary to store all the variables  
  for i in range(1,len(outlist)): 
      varout.update({outlist[i][0]: adict[outlist[i][0]]}) #update the dictionary
  
  if res=='daily': ntime_out=ntime_run
  if res=='month': ntime_out=ntime_run/365*12
  if res=='year':  ntime_out=ntime_run/365
  
  writeNC.write_var(res,pathin,pathout,filein,fileout,nlats,nlons,ntime_out,nlevs,Nani,Nhum,outlist,varout,shortnames)
 
