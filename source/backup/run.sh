#!/bin/bash

# this script runs several trials of the forge model with varied input parameters, calling the 'run.py' script for each trial. The trials are run in parallel ten at a time. The parameters that are varied are maximum human birthrate, with a default value of 0.1, maximum human fat to bodyweight ratio with a default value of 0.3, human bodyweight in kilograms with a default value of 50, rate of metabolism by anabolism with a default value of 54.6, rate of metabolism by catabolism with a default value of 39.3, annual animal background mortality rate with a default value of 1/25, and annual human background mortality rate with a default value of 1/80.

python3 run.py 'control.nc' '0.1,0.3,50,54.6,39.3,1/25,1/80' & # control; no parameters changed from default
python3 run.py 'control2.nc' '0.1,0.3,50,54.6,39.3,1/25,1/80' & # another control trial to test stability with no changed parameters
python3 run.py 'decreased_estab_hum_a(maximum_human_birthrate)_to_0.05.nc' '0.05,0.3,50,54.6,39.3,1/25,1/80' & # decreased maximum annual human birthrate to 0.05
python3 run.py 'decreased_fat_max_to_0.2_times_hum_bw.nc' '0.1,0.2,50,54.6,39.3,1/25,1/80' & # decreased maximum human fat to bodyweight ratio to 0.2
python3 run.py 'decreased_hum_bw(human_body_weight)_to_40.nc' '0.1,0.3,40,54.6,39.3,1/25,1/80' & # decreased human bodyweight to 40 kg
python3 run.py 'decreased_m_anabolism(anabolism)_to_34.6.nc' '0.1,0.3,50,34.6,39.3,1/25,1/80' & # decreased rate of anabolism to 34.6
python3 run.py 'decreased_m_catabolism(catabolism)_to_29.3.nc' '0.1,0.3,50,54.6,29.3,1/25,1/80' & # decreased rate of catabolism to 29.3
python3 run.py 'decreased_mort_ani(animal_death_rate)_to_30^-1.nc' '0.1,0.3,50,54.6,39.3,1/30,1/80' & # decreased annual animal death rate to 1/30
python3 run.py 'decreased_mort_hum(background_human_death_rate)_to_100^-1.nc' '0.1,0.3,50,54.6,39.3,1/25,1/100' & # decreased annual human death rate to 1/100
python3 run.py 'increased_estab_hum_a(maximum_human_birthrate)_to_0.15.nc' '0.15,0.3,50,54.6,39.3,1/25,1/80' & # increased maximum annual human birthrate to 0.15
wait # wait for all previous processes to finish before continuing
python3 run.py 'increased_fat_max_to_0.4_times_hum_bw.nc' '0.1,0.4,50,54.6,39.3,1/25,1/80' & # increased maximum human fat to bodyweight ratio to 0.4
python3 run.py 'increased_hum_bw(human_body_weight)_to_60.nc' '0.1,0.3,60,54.6,39.3,1/25,1/80' & # increased human bodyweight to 60 kg
python3 run.py 'increased_m_anabolism(anabolism)_to_74.6.nc' '0.1,0.3,50,74.6,39.3,1/25,1/80' & # increased rate of anabolism to 74.6
python3 run.py 'increased_m_catabolism(catbolism)_to_49.3.nc' '0.1,0.3,50,54.6,49.3,1/25,1/80' & # increased rate of catabolism to 49.3
python3 run.py 'increased_mort_ani(background_animal_deathrate)_to_20^-1.nc' '0.1,0.3,50,54.6,39.3,1/20,1/80' & # increased annual animal death rate to 1/20
python3 run.py 'increased_mort_hum(background_human_death_rate)_to_60^-1.nc' '0.1,0.3,50,54.6,39.3,1/25,1/60' & # increased annual human death rate to 1/60


