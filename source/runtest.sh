#!/bin/bash

# runs FORGE model with only four pixels once with default parameters to test functionality of the model by calling the 'runtest.py' module once

python3 runtest.py 'control.nc' '0.1,0.3,50,54.6,39.3,1/25,1/80' & # calls 'runtest' module with default parameters as positional argument
#python3 run.py 'control2.nc' '0.1,0.3,50,54.6,39.3,1/25,1/80' &
#python3 run.py 'decreased_estab_hum_a(maximum_human_birthrate)_to_0.05.nc' '0.05,0.3,50,54.6,39.3,1/25,1/80' &
#python3 run.py 'decreased_fat_max_to_0.2_times_hum_bw.nc' '0.1,0.2,50,54.6,39.3,1/25,1/80' &
#python3 run.py 'decreased_hum_bw(human_body_weight)_to_40.nc' '0.1,0.3,40,54.6,39.3,1/25,1/80' &
#python3 run.py 'decreased_m_anabolism(anabolism)_to_34.6.nc' '0.1,0.3,50,34.6,39.3,1/25,1/80' &
#python3 run.py 'decreased_m_catabolism(catabolism)_to_29.3.nc' '0.1,0.3,50,54.6,29.3,1/25,1/80' &
#python3 run.py 'decreased_mort_ani(animal_death_rate)_to_30^-1.nc' '0.1,0.3,50,54.6,39.3,1/30,1/80' &
#python3 run.py 'decreased_mort_hum(background_human_death_rate)_to_100^-1.nc' '0.1,0.3,50,54.6,39.3,1/25,1/100' &
#python3 run.py 'increased_estab_hum_a(maximum_human_birthrate)_to_0.15.nc' '0.15,0.3,50,54.6,39.3,1/25,1/80' &
#wait
#python3 run.py 'increased_fat_max_to_0.4_times_hum_bw.nc' '0.1,0.4,50,54.6,39.3,1/25,1/80' &
#python3 run.py 'increased_hum_bw(human_body_weight)_to_60.nc' '0.1,0.3,60,54.6,39.3,1/25,1/80' &
#python3 run.py 'increased_m_anabolism(anabolism)_to_74.6.nc' '0.1,0.3,50,74.6,39.3,1/25,1/80' &
#python3 run.py 'increased_m_catabolism(catbolism)_to_49.3.nc' '0.1,0.3,50,54.6,49.3,1/25,1/80' &
#python3 run.py 'increased_mort_ani(background_animal_deathrate)_to_20^-1.nc' '0.1,0.3,50,54.6,39.3,1/20,1/80' &
#python3 run.py 'increased_mort_hum(background_human_death_rate)_to_60^-1.nc' '0.1,0.3,50,54.6,39.3,1/25,1/60' &


