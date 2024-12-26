## This repo contains the following:
- data: This folder contains the data for few basins to try the HYMOD model
- output: This folder contains two subfolders 'nse' and 'parameter' which stores the output nse and parameter values after model calibration
- hymod_model.py: This script contains the HYMOD model
- hymod_calibrate.py: This scripts contains a function to calibrate the HYMOD model using genetic algorithm, the script then reads the data from the `data` folder and calibrates the model for each basin and stores the output in the `output` folder
- station_id.csv: This file contains the station id for few basins to try the HYMOD model