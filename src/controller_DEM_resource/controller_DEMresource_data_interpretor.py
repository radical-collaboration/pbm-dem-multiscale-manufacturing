#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:32:29 2017

@author: chai
"""

import numpy as np
import controller_DEMresource_data_reader as DEMreader


class controller_DEM_resource_interpretor(object):

    def __init__(self, timestep, types, init_ts):
        self.timestep = timestep
        self.type = types
        self.init_timestep = init_ts
        self.num_of_particles = 0
        self.avg_impacts = 0
        self.tot_part_each_type = np.zeros(self.type)
        self.avg_vel_array = np.zeros(self.type)
        self.collision_matrix = np.zeros((self.type, self.type))
        self.init_impacts = 0
        self.init_collision_matrix = np.zeros((self.type, self.type))
        self.init_avg_vel = np.zeros(self.type)
        self.init_num_particles = 0

# method to average the velocity of the collision data
    def avg_all_data(self, ts): # here ts is the time step
        obj_data_reader = DEMreader.Controller_DEM_resource_reader(ts, self.type)
        temp_coll_data = np.zeros_like(self.collision_matrix)
        self.num_of_particles = obj_data_reader.number_of_particles
        for x in xrange(0, self.type):
            self.tot_part_each_type[x] = len(obj_data_reader.collision_data_acc_types[str(x+1)])
            temp1 = self.tot_part_each_type[x]
            temp2 = obj_data_reader.collision_data_acc_types[str(x+1)]
            vx = 0.0
            vy = 0.0
            vz = 0.0
            for i in xrange(0, int(temp1)):
                #print(obj_data_reader.collision_data_acc_types[str(x+1)][int(temp1)-1])
                #print(temp2[int(temp1 - 1)]['velocity'][0])
                vz = vz + float(temp2[int(i)]['velocity'][0][2])
                vx = vx + float(temp2[int(i)]['velocity'][0][0])
                vy = vy + float(temp2[int(i)]['velocity'][0][1])
                for j in xrange(0,self.type):
                    temp_coll_data[x][j] += int(temp2[int(i)]['collision'][0][j])
            v_avg = (vx ** 2 + vy ** 2 + vz ** 2) ** (0.5)
            self.avg_vel_array[x] = v_avg / self.tot_part_each_type[x]
        self.collision_matrix = temp_coll_data
        self.avg_impacts = obj_data_reader.number_of_impacts
        print("Recalculated for timestep %d"%ts)

# method to data of the initial time step
    def init_data_calculations(self):
        init_data_obj = DEMreader.Controller_DEM_resource_reader(self.init_timestep, self.type)
        self.init_num_particles = init_data_obj.number_of_particles
        i1 = self.avg_all_data(self.init_timestep)
        self.init_impacts = init_data_obj.number_of_impacts
        self.init_avg_vel = self.avg_vel_array
        self.init_collision_matrix = self.collision_matrix
        self.init_num_particles = self.num_of_particles
        print("Intial data at timestep %d calculated"%self.init_timestep)

# method to compare results
    def liggghts_data_comparison(self, ts):
        dem_timestep = 5e-7
        min_time_diff = 0.2
        flag = 0  # 0 keeps it running, 1 to change from DEM to PBM and 2 to reaches steady state and quit
        dump_difference = 50000
        min_timestep_diff = min_time_diff/dem_timestep
        if(ts - self.init_timestep > min_timestep_diff):
            #temp variables defined for comparison
            vai = sum(self.init_avg_vel) / self.type
            cai = sum(sum(self.init_collision_matrix))
            iai = float(self.init_impacts)
            a1 = self.avg_all_data(ts)
            va1 = sum(self.avg_vel_array) / self.type
            ca1 = sum(sum(self.collision_matrix))
            ia1 = float(self.avg_impacts)
            a2 = self.avg_all_data(ts + dump_difference)
            va2 = sum(self.avg_vel_array) / self.type
            ca2 = sum(sum(self.collision_matrix))
            ia2 = float(self.avg_impacts)
            a3 = self.avg_all_data(ts - dump_difference)
            va3 = sum(self.avg_vel_array) / self.type
            ca3 = sum(sum(self.collision_matrix))
            ia3 = float(self.avg_impacts)
            avg_vel = (va1 + va2 + va3) / 3
            avg_coll = (ca1 + ca2 + ca3) / 3
            avg_imp = (ia1 + ia2 + ia3) / 3
            vel_comp = (avg_vel - vai) / vai
            collision_comp = (avg_coll - cai) / cai
            impact_comp = (avg_imp - iai) / iai
            if(vel_comp > 0.1 or collision_comp > 0.1 or impact_comp > 0.1):
                flag = 1
            elif((ts - self.init_timestep) > (5 / dem_timestep)):
                flag = 2
            else:
                flag = 0
        return flag





# ------------------------------------------------------------------------------------------------

abcd = controller_DEM_resource_interpretor(500000, 16, 2000000)
# a1 = abcd.avg_all_data(5000000)
# print(abcd.num_of_particles)
# abcd.init_data_calculations()
# print(abcd.num_of_particles)
# a1 = abcd.avg_all_data(6000000)
# print(abcd.num_of_particles)
print(abcd.liggghts_data_comparison(5000000))