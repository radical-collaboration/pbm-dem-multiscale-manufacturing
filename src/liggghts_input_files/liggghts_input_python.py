#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 09:45:37 2017

@author: Chaitanya Sampat

############ LIGGGHTS input file creator in Python ###################
This file creates the liggghts input file using the inital timestep and
"""
import numpy as np
import numpy.matlib
import random
from scipy.stats import norm

class liggghts_input_creator(object):

    def __init__(self, init_timestep, final_timestep, min_dia, max_dia, types_of_particles, total_flow_rate, solid_density):
        self.init_timestep = init_timestep
        self.final_timestep = final_timestep
        self.min_dia = min_dia
        self.max_dia = max_dia
        self.types_of_particles = types_of_particles
        self.flow_rate = total_flow_rate
        self.density = solid_density
        self.main_writer()

    def isPrime(self, n):
        i = 2
        while (i * i <= n):
            if (n % i == 0):
                return False
            i += 1
        return True


    def main_writer(self):
        min_d = self.min_dia
        max_d = self.max_dia
        top = self.types_of_particles
        in_ts = self.init_timestep
        fin_ts = self.final_timestep
        if (top > 0 and min_d > 0 and max_d >= min_d and fin_ts > in_ts):
            radS1 = np.linspace(min_d, max_d, num=top)
            radS2 = np.linspace(min_d, max_d, num=top)
            radL1 = np.linspace(min_d, max_d, num=top)
            radG1 = np.linspace(min_d, max_d, num=top)
            totflowrate = self.flow_rate
            radS11, radS22 = np.meshgrid(radS1, radS2);
            mode = len(radS1)

            #if in_ts == 0:
            in_file = open("in.2_sim_new", "w")
            in_file.write("# This is the initial input file for the DEM LIGGGHTS simulation\n\n\n\n\n")
            #else:
            #in_file = open("in.restart_from_%d"%in_ts, "w")
            #in_file.write(" This the restart input fiel for the DEM LIGGGHTS simulation\n\n\n\n\n")
            in_file.write("modify_timing      on\n")
            in_file.write("atom_style         granular\n")
            in_file.write("atom_modify        map array\n")
            in_file.write("boundary           f f f\n")
            in_file.write("newton             off\n")
            in_file.write("communicate        single vel yes\n")
            in_file.write("units              si\n")
            in_file.write("region             reg block -0.01 0.450 -0.150 0.150 -0.080 0.100\n")
            #if in_ts == 0:
            in_file.write("processors         * * *\n")
            in_file.write("create_box         %d reg\n"%(mode+1))
            #else:
            #    in_file.write("read_restart   restart/granulator.%d.restart\n"%in_ts)
            in_file.write("neighbor           0.0005 bin\n")
            in_file.write("neigh_modify       delay 0 check yes page 10000000 one 30000\n")
            in_file.write("timestep           5e-7\n")
            in_file.write("fix                gravi all gravity 9.81 vector 0.0 -1.0 0.0\n\n\n")
            in_file.write("pair_style         gran model hertz tangential history\n")
            in_file.write("pair_coeff         * *\n\n\n")

            in_file.write("# Material properties required for new pair styles\n\n\n")
            ym = np.zeros(mode+1)
            ym[0:mode] = numpy.matlib.repmat(8e6, 1, mode)
            ym[16] = 9.9e8
            in_file.write("fix                m1 all property/global youngsModulus peratomtype ")
            for item in ym:
                in_file.write("%f "%item)
            in_file.write("\n")

            poissonsRatio = np.linspace(0.2, 0.2, num=mode+1)
            in_file.write("fix                m2 all property/global poissonsRatio peratomtype ")
            for item in poissonsRatio:
                in_file.write("%f "%item)
            in_file.write("\n")

            coef_restitution = np.linspace(0.4, 0.4, num=(mode+1)**2)
            in_file.write("fix                m3 all property/global coefficientRestitution peratomtypepair %d "%(mode+1))
            for item in coef_restitution:
                in_file.write("%f "%item)
            in_file.write("\n")

            coef_friction = np.linspace(0.5, 0.5, num=(mode+1)**2)
            in_file.write("fix                m4 all property/global coefficientFriction peratomtypepair %d "%(mode+1))
            for item in coef_friction:
                in_file.write("%f "%item)
            in_file.write("\n")

            roll_friction = np.linspace(0.02, 0.02, num=(mode+1)**2)
            in_file.write("fix                m5 all property/global coefficientRollingFriction peratomtypepair %d "%(mode+1))
            for item in roll_friction:
                in_file.write("%f "%item)
            in_file.write("\n\n\n\n\n")
            in_file.write("# Cad import and meshing\n")
            in_file.write("fix                cad1 all mesh/surface file impeller_coarse.stl element_exclusion_list read impeller heal auto_remove_duplicates precision 1e-9 type %d scale 1e-3 curvature 1.0e-10\n"%(mode+1))
            in_file.write("fix                cad2 all mesh/surface file shell_closed.stl element_exclusion_list read shell heal auto_remove_duplicates precision 1e-9 type %d scale 1e-3  curvature 1.0e-10\n\n\n"%(mode+1))
            in_file.write("fix                granulator all wall/gran model hertz tangential history mesh n_meshes 2 meshes cad2 cad1\n\n\n")
            in_file.write("#Rotating the impeller\n")
            in_file.write("variable           w equal -2000\n")
            in_file.write("fix                movecad1 all move/mesh mesh cad1 rotate/variable origin 0.193 0.00 0.00 axis 1. 0. 0. omega v_w\n\n\n")
            in_file.write("#Region and insertion\n")
            in_file.write("region             factory1 block 0.007 0.047 0.04 0.06 -0.02 0.02 units box\n\n\n")
            in_file.write("# Grouping\n")
            for x in xrange(0,mode):
                in_file.write("group              type%d_ type %d\n"%(x+1,x+1))
            in_file.write("\n\n\n\n")
            in_file.write("# Particle Specification\n")            
            minPrime = 10000
            maxPrime = 1000000
            cached_primes = [i for i in range(minPrime,maxPrime) if self.isPrime(i)]
            for x in xrange(0,mode):
                p = 10000
                q = 100000
                n = random.choice([i for i in cached_primes if p<i<q])
                in_file.write("fix                pts%d type%d_ particletemplate/sphere %d atom_type %d density constant %d radius constant %f\n"%(x+1,x+1,n,x+1,self.density,radS1[x]))
            in_file.write("\n\n\n\n")
            x = np.linspace(0, 1, mode)
            ploon = 0
            for i in x:
                ploon += (i - np.mean(x)) ** 2 / (mode - 1)
            ploon = np.sqrt(ploon)
            y = norm.pdf(x, (radS1[mode - 1] - radS1[1]), ploon)

            wt = [(j / sum(y)) for j in y]
            in_file.write("\n#Partcle Distribution for insertion\n")
            massflowrateinseconds = totflowrate / 3600.0
            st = ''
            for x in xrange(0,mode):
                st = st + 'pts%d'%(x+1) + ' ' + str(wt[x]) + ' '
            n = random.choice([i for i in cached_primes if minPrime<i<maxPrime])
            in_file.write("fix                pdd1 all particledistribution/discrete %d %d %s\n"%(n, mode, st))
            in_file.write("\n\n\n\n")
            n = random.choice([i for i in cached_primes if minPrime<i<maxPrime])
            in_file.write("# Particle insertion\n")
            in_file.write("fix                ins1 all insert/rate/region seed %d distributiontemplate pdd1 nparticles INF massrate %f insert_every 10000 overlapcheck yes vel constant 0.0 -1.0 -0.0001 region factory1 ntry_mc 1000\n\n\n\n"%(n,massflowrateinseconds))

            in_file.write("# Calculating particle wall collisions\n")
            in_file.write("compute            pwc all wall/gran/local id vel force contactArea delta # Calculate particle wall collision\n")
            in_file.write("compute            ppc all pair/gran/local id vel force contactArea delta # Calculate particle-particle collision\n")

            in_file.write("# Apply nve integration to all particles that are inserted as single particles\n")
            in_file.write("fix                integ all nve/sphere\n\n\n\n")

            in_file.write("# Collecting particle collision data \n")
            in_file.write("fix                fppacc all property/atom fppacc scalar yes yes yes 0\n")
            in_file.write("\n\n# Collecting particle-particle collision data\n")
            for x in xrange(0,mode):
                in_file.write("compute            cc_%d type%d_ contact/atom\n"%(x+1,x+1))
            in_file.write("\n\n\n\n# Dump files configurationn\n")
            in_file.write("run                1\n")
            in_file.write("dump               myDump  all custom 50000 post/dump*.atom id type x y z ix iy iz vx vy vz fx fy fz c_cc_1 c_cc_2 c_cc_3 c_cc_4 c_cc_5 c_cc_6 c_cc_7 c_cc_8 c_cc_9 c_cc_10 c_cc_11 c_cc_12 c_cc_13 c_cc_14 c_cc_15 c_cc_16 f_fppacc radius\n")
            in_file.write("dump               myDump2 all local 50000 post/coll*.atom c_pwc[1] c_pwc[2] c_pwc[3]\n")
            in_file.write("dump               myDump3 all local 50000 post/pcoll*.atom c_ppc[1] c_ppc[2] c_ppc[3] c_ppc[4]\n")
            in_file.write("dump               dumpstl1 all mesh/stl 50000 post/dump*.stl\n")
            in_file.write("restart            50000 restart/granulator.*.restart\n\n\n\n\n")
            in_file.write("\n# End run\n")
            in_file.write("run                %d upto\n"%(fin_ts + 1))
            in_file.write("\n# Unfixing data points to facilitate restart\n")
            in_file.write("unfix              ins1\n")
            for x in xrange(1,mode):
                in_file.write("unfix              pts%d\n"%(x+1))
            in_file.close()
# -------------------------------------------------------------------------------------------------

a = liggghts_input_creator(0,2000000,0.002,0.002,16,15,500)
#in_file.write("\n")
