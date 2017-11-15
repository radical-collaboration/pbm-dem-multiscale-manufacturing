#!/usr/bin/env python
from __future__ import division
__copyright__ = 'Copyright 2017-2018, http://radical.rutgers.edu'
__license__   = 'MIT'
__author__    = 'Ioannis Paraskevakos'


import os
import sys
import warnings
import radical.utils as ru
import radical.pilot as rp


class Executor(object):
    """
    The Executor class is responsible to implement the two way coupling
    of the DEM/PBM execution. A class object can either be initialized
    at creation by passing a config json file or via the config method.

    For more details on the contents of the json file, please see help
    of config method
    """

    def __init__ (self,config=None):
        """
        The constructor
        """

        # If there is no configuration file passed initialized everything
        # to None, otherwise to the values they have from the json file.
        # Attribute _session_name        : The name of the RADICAL-Pilot Session
        # Attribute _resource            : The resource that will be used for the execution
        # Attribute _cores               : The number of cores that will be used
        # Attribute _queue               : The queue that will be used 
        # Attribute _project             : The value of the project that will be charged
        # Attribute _pathtoLIGGHTSinputs : Absolute path to the input file for LIGGGHTS
        # Attribute _pathtoLIGGHTS       : Absolute path to LIGGGHTS executable to resource
        # Attribute _pathtoPBMexecutable : Absolute path to PBM executable to resource
        # Attribute _runtime             : Time resources are requested for. The value should be in minutes
        # Attribute _DEMcores            : The number of cores requested by DEM
        # Attribute _PBMcores            : The number of PBM cores
        # attribute _timesteps           : Number of  maximum timestamps that DEM run
        # attribute _diameter            : Diameter of the particles
        # Attribute _PBMs                : The number of concurrent PBM executions
        # Attribute _DEMs                : The number of concurrent DEM executions
        # Attribute _types               : Number of types of solids
        
        if not config:
            # Initialize everything to None when not a config file present
            self._session_name        = None
            self._resource            = None
            self._cores               = None
            self._queue               = None
            self._project             = None
            self._pathtoLIGGHTSinputs = None
            self._pathtoLIGGHTS       = None
            self._pathtoPBMexecutable = None
            self._runtime             = None
            self._DEMcores            = None
            self._PBMcores            = None
            self._timesteps           = None
            self._diameter            = None
            self._PBMs                = 1
            self._DEMs                = 1
        else:
            # Get the values of the attributes based on the config file
            try:
                conf = ru.read_json(config)
                self._session_name        = config['session']
                self._resource            = config['resource']
                self._cores               = config['cores']
                self._queue               = config['queue']
                self._project             = config['project']
                self._pathtoLIGGHTSinputs = config['pathtoLIGGHTSinputs']
                self._pathtoLIGGHTS       = config['pathtoLIGGHTS']
                self._pathtoPBMexecutable = config['pathtoPBMexecutable']
                self._DEMcores            = config['DEMcores']
                self._PBMcores            = config['PBMcores']
                self._timesteps           = config['timesteps']
                self._diameter            = config['diameter']
            except IOError as e:
                raise RuntimeError('%s is not a valid configuration file'%config)
            except KeyError as e:
                raise RuntimeError('Key %s does not exist in config file'%e)

            # If there is a key named PBMnumber use the assigned value, else spawn
            # as many as possible in the cores that DEM used.
            self._PBMs                = config['PBMnumber'] if conf.has_key('PBMnumber') else self._DEMcores/self._PBMcores
            self._DEMs                = 1
            self._types               = config['types'] if conf.has_key('types') else 16
            self._bins1               = config['bins1'] if conf.has_key('bins1') else 16
            self._bins2               = config['bins2'] if conf.has_key('bins2') else 16
            self._diff_DEM            = config['diff_DEM'] if conf.has_key('diff_DEM') else 50000

        # Initialize the Pilot Manager, Unit Manager and Session to None, until those objects are 
        # being created.
        self._pmgr = None
        self._umgr = None
        self._session = None

    def configure(self,config):
        """
        This method configures the executor. It needs a valid configuration file, otherwise
        it raises a runtime error. The configuration file should contain the following 
        dictionary:
        {
            'session'             : The name of the RADICAL-Pilot Session
            'resource'            : The resource that will be used for the execution
            'cores'               : The number of cores that will be used
            'queue'               : The queue that will be used 
            'project'             : The value of the project that will be charged
            'pathtoLIGGHTSinputs' : Absolute path to the input file for LIGGGHTS
            'pathtoLIGGHTS'       : Absolute path to LIGGGHTS executable to resource
            'pathtoPBMexecutable' : Absolute path to PBM executable to resource
            'runtime'             : Time resources are requested for. The value should be in minutes
            'DEMcores'            : The number of cores requested by DEM
            'PBMcores'            : The number of PBM cores
            'PBMs'                : The number of concurrent PBM executions
            'DEMs'                : The number of concurrent DEM executions
            'timesteps'           : Number of  maximum timestamps that DEM run
            'diameter'            : Diameter of the particles
        }

        Additional values are:
                * 'PBMnumber' : set's up the number of PBM instances
                * 'DEMnumber' : Set's up the number of DEM instances
        """

        try:
            conf = ru.read_json(config)
            self._session_name        = config['session']
            self._resource            = config['resource']
            self._cores               = config['cores']
            self._queue               = config['queue']
            self._project             = config['project']
            self._pathtoLIGGHTSinputs = config['pathtoLIGGHTSinputs']
            self._pathtoLIGGHTS       = config['pathtoLIGGHTS']
            self._pathtoPBMexecutable = config['pathtoPBMexecutable']
            self._DEMcores            = config['DEMcores']
            self._PBMcores            = config['PBMcores']
            self._timesteps           = config['timesteps']
            self._diameter            = config['diameter']
        except IOError as e:
            raise RuntimeError('%s is not a valid configuration file'%config)
        except KeyError as e:
            raise RuntimeError('Key %s does not exist in config file'%e)


        # If there is a key named PBMnumber use the assigned value, else spawn
        # as many as possible in the cores that DEM used.
        self._PBMs                = config['PBMnumber'] if conf.has_key('PBMnumber') else self._DEMcores/self._PBMcores
        self._DEMs                = 1
        self._types               = config['types'] if conf.has_key('types') else 16
        self._bins1               = config['bins1'] if conf.has_key('bins1') else 16
        self._bins2               = config['bins2'] if conf.has_key('bins2') else 16
        self._diff_DEM            = config['diff_DEM'] if conf.has_key('diff_DEM') else 50000

    def _start(self):
        """
        This method start the resource allocation and creates the necessary objects for
        executing on some resource.

        It does not submits any tasks to the resource for execution.
        """
        try:
            # Initially create a session with the session name the user provided.
            self._session = rp.Session(uid=self._session_name)

            # Create a pilot manager object for this session
            self._pmgr = rp.PilotManager(session=self._session)

            # Create a UnitManager object for this session
            self._umgr = rp.UnitManager(session=self._session)

            # Create the Pilot description based on the configuration that the user passed through
            # the config file.

            pd_desc = {
                        'resource'      : self._resource,
                        'runtime'       : self._runtime,  # pilot runtime (min)
                        'exit_on_error' : True            # In case of an error, the pilot exits
                        'project'       : self._project,
                        'queue'         : self._queue,
                        'cores'         : self._cores
                      }

            pdesc = rp.ComputePilotDescription(pd_init)

            # Launch the pilot.
            self._pmgr.submit_pilots(pdesc)
            # Register the ComputePilot in a UnitManager object.
            self._umgr.add_pilots(pilot)

            return 0
        except Exception as e:
            # Something unexpected happened in the pilot code above
            raise RuntimeError('caught Exception %s\n'%e)

        except (KeyboardInterrupt, SystemExit) as e:
            # the callback called sys.exit(), and we can here catch the
            # corresponding KeyboardInterrupt exception for shutdown.  We also catch
            # SystemExit (which gets raised if the main threads exits for some other
            # reason).
            warnings.warn('exit requested\n',UserWarning)

    def _start_dem_units(self,timestep=0,types=0,restart=False):
        """
        This method creates DEM units, their monitors and submits them for execution. 
        It also returns the unit objects so that they can be monitored. This method takes
        two arguments, timestep and type
        timestep : the timestep DEM monitoring should start
        type     : 
        restrt   : defines if the unit that will be created is restarting the DEM simulation
                   or not.
        """

        for i in range(self._DEMs):
            cud = rp.ComputeUnitDescription()
            cud.environment    = ['PATH='+self._pathtoLIGGHTS+':$PATH']
            cud.pre_exec       = ['mkdir CSVs','mkdir post','mkdir restart']
            cud.executable     = 'lmp_micstam'
            cud.arguments      = ['-in','in.*' ]
            cud.input_staging  = [{'source': self._pathtoLIGGHTSinputs+'/in.final_%fmm'%self._diameter,
                                   'target':'unit:///in.final_%fmm'%self._diameter,
                                   'action'  :rp.LINK},
                                  {'source': self._pathtoLIGGHTSinputs+'/shell_closed.stl',
                                   'target':'unit:///shell_closed.stl',
                                   'action'  :rp.LINK},
                                  {'source': self._pathtoLIGGHTSinputs+'/shell',
                                   'target':'unit:///shell',
                                   'action'  :rp.LINK},
                                  {'source': self._pathtoLIGGHTSinputs+'/impeller',
                                   'target':'unit:///impeller',
                                   'action'  :rp.LINK},
                                  {'source': self._pathtoLIGGHTSinputs+'/impeller_coarse.stl',
                                   'target':'unit:///impeller_coarse.stl',
                                   'action'  :rp.LINK},
                                ]
            cud.output_staging = [{'source': 'unit:///post/collision%d.atom'%self._timesteps,
                                   'target': 'pilot:///collision%d.atom'%self._timesteps,
                                   'action'  : rp.LINK},
                                  {'source': 'unit:///post/impact%d.atom'%self._timesteps,
                                   'target': 'pilot:///impact%d.atom'%self._timesteps,
                                   'action'  : rp.LINK},
                                   {'source': 'unit:///post/collision%d.atom'%self._timesteps-self._diff_DEM,
                                   'target': 'pilot:///collision%d.atom'%self._timesteps-self._diff_DEM,
                                   'action'  : rp.LINK},
                                  {'source': 'unit:///post/impact%d.atom'%self._timesteps-self._diff_DEM,
                                   'target': 'pilot:///impact%d.atom'%self._timesteps-self._diff_DEM,
                                   'action'  : rp.LINK},
                                   {'source': 'unit:///post/collision%d.atom'%self._timesteps-2*self._diff_DEM,
                                   'target': 'pilot:///collision%d.atom'%self._timesteps-2*self._diff_DEM,
                                   'action'  : rp.LINK},
                                  {'source': 'unit:///post/impact%d.atom'%self._timesteps-2*self._diff_DEM,
                                   'target': 'pilot:///impact%d.atom'%self._timesteps-2*self._diff_DEM,
                                   'action'  : rp.LINK},]
            cud.cores          = self._DEMcores
            cud.mpi            = True

            # Submit the first unit and wait until it staged its input files. Waiting is
            # needed so that we can get the path of the unit and pass it to the DEM monitor.
            dem_unit = self._umgr.submits_units(cud)
            # This line blocks the execution until the DEM unit has a path.
            self._umgr.wait_units(uids=dem_unit.uid,state=rp.states.AGENT_SCHEDULING_PENDING)

            #Now that we have a path we can continue
            cud2 = rp.ComputeUnitDescription()
            cud2.executable = 'python'
            cud2.arguments = ['controller_DEMresource_main.py',timestep,self._types,ru.Url(dem_unit.sandbox).path]
            cud2.input_staging = [{'source':'file:///controller_DEMresource_main.py'
                                  'target':'unit:///controller_DEMresource_main.py'
                                  'action': rp.TRANSFER},
                                  {'source':'file:///controller_DEMresource_data_reader.py'
                                  'target':'unit:///controller_DEMresource_data_reader.py'
                                  'action': rp.TRANSFER},
                                  {'source':'file:///controller_DEMresource_data_interpretor.py'
                                  'target':'unit:///controller_DEMresource_data_interpretor.py'
                                  'action': rp.TRANSFER}]
            cud2.output_staging = [{'source': 'unit:///DEM_status.json',
                                    'target': 'file:///DEM_status.json',
                                    'action'  : rp.TRANSFER}]

            # Submit the monitor unit and return
            dem_monitor_unit = self._umgr.submits_units(cud2)


            self._dem_unit = dem_unit
            self._dem_monitor_unit = dem_monitor_unit

    def _check_DEM_status(self,status_file):
        """
        This method checks the state of the DEM execution and based on the status code
        reports whether or not the simulation should continue with  starting PBM or 
        stop everything.

        Inputs:
        status_file: The name of a json file that contains the status and other entries
                     necessary for the PBM execution

        Returns:
        cont: A boolean value that shows if the simulation should continue or not
        arguments: A list that contains the necessary arguments for the PBM execution.
        """

        status_fid = open(status_file)
        status_dict = json.load(status_fid)
        status_fid.close()

        if status_dict['status'] == 1:
            cont = True
        else:
            cont = False
        
        pbm_timestep      = status_dict['DEM_time_step']
        pbm_init_timestep = status_dict['PBM_init_time_step']
        pbm_mixing_time   = status_dict['mixing_times']

        return cont,pbm_timestep,pbm_init_timestep,pbm_mixing_time

    def _check_PBM_status(self,status_file):
        """
        This method checks the state of the DEM execution and based on the status code
        reports whether or not the simulation should continue with  starting PBM or 
        stop everything.

        Inputs:
        status_file: The name of a json file that contains the status and other entries
                     necessary for the PBM execution

        Returns:
        cont: A boolean value that shows if the simulation should continue or not
        arguments: A list that contains the necessary arguments for the PBM execution.
        """

        status_fid = open(status_file)
        status_dict = json.load(status_fid)
        status_fid.close()
        if status_dict['status'] == 1:
            cont = True
        else:
            cont = False

        return cont,status_dict['last_time_step']


    def _start_pbm_units(self,timestep=0, init_timestep=0,\
                         mixing_time=0):

        try:
            # Get the path to the DEM folder. This will allow correct staging
            # input for the PBMs
            dem_monitor_path = ru.Url(self._dem_monitor_unit.sandbox).path

            # Create a list of PBM units and insert every description in that list
            pbm_cud_list = list()
            for i in range(self._PBMs):
                stage_in = [{'source': 'pilot:///collision%d.atom'%timesteps,
                             'target': 'unit:///sampledumpfiles/collision%d.%d_%d'%(timesteps,self._DEMcores,self._diameter),
                             'action'  : rp.LINK},
                            {'source': 'pilot:///impact%d.atom'%timesteps,
                             'target': 'unit:///sampledumpfiles/impact%d.%d_%f'%(timesteps,self._DEMcores,self._diameter),
                             'action'  : rp.LINK},
                            {'source': 'pilot:///collision%d.atom'%timesteps-self._diff_DEM,
                             'target': 'unit:///sampledumpfiles/collision%d.%d_%d'%(timesteps-self._diff_DEM,self._DEMcores,self._diameter),
                             'action'  : rp.LINK},
                            {'source': 'pilot:///impact%d.atom'%timesteps-self._diff_DEM,
                             'target': 'unit:///sampledumpfiles/impact%d.%d_%f'%(timesteps-self._diff_DEM,self._DEMcores,self._diameter),
                             'action'  : rp.LINK},
                            {'source': 'pilot:///collision%d.atom'%timesteps-2*self._diff_DEM,
                             'target': 'unit:///sampledumpfiles/collision%d.%d_%d'%(timesteps-2*self._diff_DEM,self._DEMcores,self._diameter),
                             'action'  : rp.LINK},
                            {'source': 'pilot:///impact%d.atom'%timesteps-2*self._diff_DEM,
                             'target': 'unit:///sampledumpfiles/impact%d.%d_%f'%(timesteps-2*self._diff_DEM,self._DEMcores,self._diameter),
                              'action'  : rp.LINK}]
        
                cud = rp.ComputeUnitDescription()
                cud.environment    = ['PATH='+self._pathtoLIGGHTS+':$PATH']
                cud.pre_exec       = ['mkdir csvDump','mkdir txtDumps']
                cud.executable     = 'model.out'
                cud.arguments      = [dem_monitor_path+'PBM_input.in',self._DEMcores,self._diameter]
                cud.input_staging  = stage_in
                cud.output_staging = [{'source': 'unit:///csvDump.tar.gz',
                                       'target': 'client:///csvDump.tar.gz',
                                       'action'  : rp.TRANSFER}]

                cud.cores          = self._PBMcores
                cud.mpi            = True
                pbm_cud_list.append(cud)

            # Submit them to the agent and wait until all have a path
            pbm_uids = self._umgr.submits_units(pbm_cud_list)
            self._umgr.wait_units(uids=pbm_uids.uid,state=rp.states.AGENT_SCHEDULING_PENDING)

            pbm_monitor_cud_list = list()
            for i in range(self._PBMs):
                cud = rp.ComputeUnitDescription()
                cud.input_staging = [{'source':'file:///controllerPBMresourceMain.py'
                                       'target':'unit:///controllerPBMresourceMain.py'
                                       'action': rp.TRANSFER},
                                      {'source':'file:///controllerPBMresourceDataReader.py'
                                       'target':'unit:///controllerPBMresourceDataReader.py'
                                       'action': rp.TRANSFER},
                                      {'source':'file:///controllerPBMresourceDataInterpretor.py'
                                       'target':'unit:///controllerPBMresourceDataInterpretor.py'
                                       'action': rp.TRANSFER}]
                cud.output_staging = [{'source': 'unit:///PBM_status.json',
                                       'target': 'file:///PBM_status.json',
                                       'action'  : rp.TRANSFER}]
                cud.executable     = 'python'
                cud.arguments      = ['controllerPBMresourceMain.py',\
                                      init_timestep,self._bins1,self._bins2,\
                                      ru.Url(pbm_uids[i].sandbox).path+'csvDump',mixing_time]
                cud.cores          = 1
                cud.mpi            = False
                pbm_cud_list.append(cud)

            # Submit the monitor unit and return
            pbm_monitor_unit = self._umgr.submits_units(cud2)

            self._pbm_unit = pbm_uids
            self._pbm_monitor_unit = pbm_monitor_unit

    def _shutdown(self):
        self._session.close()

    def run(self):


        # Start the RP session and setup the selected resource.
        self._start()


        cont = True
        #Does DEM restart or not
        restart = False
        #Timestep input for DEM.
        ts=0

        while cont:

            self._start_dem_units(timestep=ts,restart=restart)

            self._umgr.wait_units(uid=self._dem_monitor_unit.uid)

            # Check DEM status returns whether the execution should continue
            # or not. It also returns the timestep PBM should start.
            cont, pbm_timestep, pbm_init_timestep, pbm_mixing_time = self._check_DEM_status('DEM_status.json')

            if self._dem_unit.state() != rp.FINAL:
                self._umgr.cancel_units(uid=self._dem_unit.uid)

            if cont == True:
                self._start_pbm_units(timestep=pbm_timestep, init_timestep=pbm_init_timestep,\
                                      mixing_time=pbm_mixing_time)

                # Waits for all the PBM units to finish. Should wait only for the
                # first
                self._umgr.wait_units(uid=[cu.uid for cu in self._pbm_monitor_units])

                self._umgr.cancel_units(uid=[cu.uid for cu in self._pbm_units])

                cont, ts = self._check_PBM_status('PBM_status.json')

            if cont == True:
                restart = True

        self._shutdown()
    