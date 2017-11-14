#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:58:47 2017

@author: Chaitanya Sampat
"""
"""
This is class uses the Data Reader and Data Interpretor classes to extract the data from d50 and
the particles files and compare them with the initial time step values. If the values differ by
more than 15%, the DEM is executed if the conditions are satisfied. If the value does not change
by 15% for over 10 seconds, it is assumed that steady state has been reached that the execution
can be terminated. If neither of the 2 conditions are met then the execution continues.

==================================================================================================
OUTPUT FLAG VALUE STATUS
==================================================================================================
0 - Continue with PBM
1 - Kill PBM and execute DEM with the new files
2 - The system has reached state with the PBM execution, kill all executions
==================================================================================================

"""

import numpy as np
import controllerPBMDataInterpretor as PBMinter
import controllerPBMDataReader as PBMreader
import os
import json
import time

class controllerPBMresourceMain(object):
# class variables
    def __init__(self, init_ts, compartments, bins1, bins2, pbm_out_path, mixingtime):
        self.initial_timestep = init_ts
        self.compartments = compartments
        self.bins1 = bins1
        self.bins2 = bins2
        self.pbm_output_path = pbm_out_path
        self.mixing_time = mixingtime

    def main(self):
    # This method executes all the necessary function calls to monitor the PBM execution and decides on the status of the simulation
        if os.path.isfile(self.pbm_output_path + '/d50_' + str(self.initial_timestep) + '.csv'):
            obj_reader = PBMreader.controllerPBMDataReader(self.initial_timestep, self.compartments, self.bins1, self.bins2, self.pbm_output_path)
        else:
            time.sleep(2)
        obj_reader = PBMreader.controllerPBMDataReader(self.initial_timestep, self.compartments, self.bins1, self.bins2, self.pbm_output_path)
        obj_inter = PBMinter.controllerPBMDataInterpretor(self.initial_timestep, self.compartments, self.bins1, self.bins2, self.pbm_output_path)
        new_timestep = obj_reader.nextfile_time_finder(self.initial_timestep)
        while (new_timestep < self.mixing_time):
            new_timestep = obj_reader.nextfile_time_finder(new_timestep)
        # opening 2 files, 1 for dumping the d50 and 2nd for dumping for particles collected over time
        dump_d50 = open("d50_overtime_starting_%d.txt"%self.initial_timestep, "w")
        dump_particles = open("particles_overtime_starting_%d.txt"%self.initial_timestep, "w")
        flag = 0
        d50_temp = np.zeros(self.compartments + 1)
#        count = 0
        particles_temp = 0
        flag1 = True
        while (flag1):
#            print(count)
            new_timestep = obj_reader.nextfile_time_finder(new_timestep)
            next_d50_filename = "/d50_%.2f.csv"%new_timestep
            next_particles_filename = "/particles_%.2f.csv"%new_timestep
            d50_filecheck = self.pbm_output_path + str(next_d50_filename)
            particles_filecheck = self.pbm_output_path + str(next_particles_filename)
#                if (os.path.isfile(d50_filecheck) and os.path.isfile(particles_filecheck)):
            t_1 = obj_inter.new_data_storage(new_timestep)
#                print(obj_inter.data_comparison(new_timestep))
            flag = obj_inter.data_comparison(new_timestep)
            d50_temp = obj_reader.data_d50_extractor(new_timestep)
            particles_temp = obj_reader.data_particles_extractor(new_timestep)
            for x in range(0,len(d50_temp)):
                dump_d50.writelines("%f "%np.float(d50_temp[0][x]))
            dump_particles.write("%f %f\n"%(new_timestep,particles_temp))
#            print(new_timestep)
            if (flag == 1 or flag == 2):
                flag1 = False
#            else:
#                continue
        dump_d50.close()
        dump_particles.close()
        if (flag == 1):
            print("Kill PBM and execute DEM")
            status = {'status':str(flag)}
            out_data = {'last timestep': str(new_timestep)}
            with open('PBM_status.json' , 'w') as pbmsf:
                json.dump(status, pbmsf)
            with open('PBM_output.json', 'w') as pbmsf:
                json.dump(out_data, pbmsf)
                #pbmsf.write(str(new_timestep))
        elif (flag == 2):
            print("Kill both DEM and PBM")
            status = {'status':str(flag)}
            out_data = {'last timestep': str(new_timestep)}
            with open('PBM_status.json' , 'w') as pbmsf:
                json.dump(status, pbmsf)
            with open('PBM_output.json', 'w') as pbmsf:
                json.dump(out_data, pbmsf)
                # pbmsf.write(str(new_timestep))
            

abcd = controllerPBMresourceMain(7,4,16,16,'/home/chai/Documents/git/CyberManufacturing/src/dummy_DEM_PBM/sample_copy',5)
abcd.main()






