#import packages
import numpy as np
import pandas as pd
import os
from hymod_model import hymod #imported from local python script
from geneticalgorithm import geneticalgorithm as ga # install package first


#function that reads station ID, takes input for that ID and outputs calibrated parameters and nse values
def calibNSE(station_id):
    #read input csv file
    df = pd.read_csv(f"data/hbv_input_{station_id}.csv")
    #only used data upto 2005 for calibration
    df = df[df["year"] < 2006]
    #hbv model input
    p = df["precip"]
    temp = df["tavg"]
    date = df["date"]
    latitude = df["latitude"]
    routing = 1 # 0: no routing, 1 allows running
    q_obs = df["qobs"] #validation data / observed flow

    ##genetic algorithm for hbv model calibration
    #reference: https://github.com/rmsolgi/geneticalgorithm
    #write a function you want to minimize
    def nse(pars):
        q_sim = hymod(pars, p, temp, date, latitude, routing)
        denominator = np.sum((q_obs - (np.mean(q_obs)))**2)
        numerator = np.sum((q_obs - q_sim)**2)
        nse_value = 1 - (numerator/denominator)
        return -nse_value #minimize this (use negative sign if you need to maximize)

    varbound = np.array([[0.001,0.999], #kpwp #lower and upper bounds of the parameters
                        [0.01,1.99], #etexp
                        [5,1000], #hmax
                        [0.01,1.99], #b 
                        [0.01,0.99], #alpha
                        [0.0005,0.99], #ks
                        [1, 2000], #lmax
                        [0.5, 2], #coeff_pet
                        [0.05,10], #ddf
                        [0.5,2], #scf
                        [-1,4], #ts
                        [-1,4], #tm
                        [-1,4], #tti
                        [0, 0.2], #whc
                        [0.1,1], #crf
                        [1,10]]) #maxbas

    algorithm_param = {
        'max_num_iteration': 10, #100              # Generations, higher is better, but requires more computational time
        'max_iteration_without_improv': None,   # Stopping criterion for lack of improvement
        'population_size': 20, #200                # Number of parameter-sets in a single iteration/generation(to start with population 10 times the number of parameters should be fine!)
        'parents_portion': 0.3,                 # Portion of new generation population filled by previous population
        'elit_ratio': 0.01,                     # Portion of the best individuals preserved unchanged
        'crossover_probability': 0.3,           # Chance of existing solution passing its characteristics to new trial solution
        'crossover_type': 'uniform',            # Create offspring by combining the parameters of selected parents
        'mutation_probability': 0.01            # Introduce random changes to the offspring’s parameters (0.1 is 10%)
    }

    model = ga(function = nse,
            dimension = 16, #number of parameters to be calibrated
            variable_type= 'real',
            variable_boundaries = varbound,
            algorithm_parameters = algorithm_param)

    model.run()
    #end of genetic algorithm

    #output of the genetic algorithm/best parameters
    best_parameters = model.output_dict
    param_value = best_parameters["variable"]
    nse_value = best_parameters["function"]
    nse_value = -nse_value #nse function gives -ve values, which is now reversed here to get true nse
    #convert into a dataframe
    df_param = pd.DataFrame(param_value).transpose()
    df_param = df_param.rename(columns={0:"kpwp", 1:"etexp", 2:"hmax", 3:"bexp", 4:"alpha", 5:"ks",
                             6:"lmax", 7:"coeff_pet", 8:"ddf", 9:"scf", 10:"ts",
                             11:"tm", 12:"tti", 13:"whc", 14:"crf", 15:"maxbas"})
    df_param["station_id"] = str(station_id)
    df_nse = pd.DataFrame([nse_value], columns=["nse"])
    df_nse["station_id"] = str(station_id)
    #save as a csv file
    df_param.to_csv(f"output/parameter/param_{station_id}.csv", index = False)
    df_nse.to_csv(f"output/nse/nse_{station_id}.csv", index = False)
    #End of function

#read stationID
stationid = pd.read_csv("station_id.csv", dtype={"station_id":str})
stationid = stationid["station_id"]
#run the function for each station
for i in stationid:
    calibNSE(i)
    print(f"Calibration done for station {i}")

print("Calibration done for all stations")
#end of script
