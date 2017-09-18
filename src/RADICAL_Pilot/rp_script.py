#!/usr/bin/env python
from __future__ import division
__copyright__ = 'Copyright 2013-2014, http://radical.rutgers.edu'
__license__   = 'MIT'
__author__    = 'Ioannis Paraskevakos'


import os
import sys
import argparse
import numpy as np
import textwrap
# ------------------------------------------------------------------------------
#
# READ the RADICAL-Pilot documentation: http://radicalpilot.readthedocs.org/
#
# ------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config',help=textwrap.dedent('JSON file that contains the following information:\n\
        resource, cores,runtime,queue, path to LIGGHTS,path to LIGGHTS inputs,path to PBM executable,session,\n\
        project,timesteps, diameter, and the number of cores for DEM and PBM. Check config.json for more details'))
    parser.add_argument('--verbose', help="Verbosity level",default='REPORT')
    parser.add_argument('--profile', help="Profile value", default='FALSE')
    args = parser.parse_args()


    os.environ['RADICAL_PILOT_VERBOSE'] = args.verbose
    os.environ['RADICAL_PILOT_PROFILE'] = args.profile

    import radical.pilot as rp
    import radical.utils as ru

    # we use a reporter class for nicer output
    report = ru.LogReporter(name='radical.pilot', level=args.verbose)
    report.title('Getting Started (RP version %s)' % rp.version)

    # read the config used for resource details
    report.info('read config')
    config = ru.read_json(args.config)
    report.ok('>>ok\n')

    # Create a new session. No need to try/except this: if session creation
    # fails, there is not much we can do anyways...
    session = rp.Session(uid=config['session'])
    cfg = session.get_resource_config(config['resource'])
    new_cfg=rp.ResourceConfig(config['resource'],cfg)
    new_cfg.mpi_launch_method='IBRUN'
    session.add_resource_config(new_cfg)

    # all other pilot code is now tried/excepted.  If an exception is caught, we
    # can rely on the session object to exist and be valid, and we can thus tear
    # the whole RP stack down via a 'session.close()' call in the 'finally'
    # clause...
    try:

        report.header('submit pilots')

        # Add a Pilot Manager. Pilot managers manage one or more ComputePilots.
        pmgr = rp.PilotManager(session=session)

        # Define an [n]-core local pilot that runs for [x] minutes
        # Here we use a dict to initialize the description object
        pd_init = {
                'resource'      : config['resource'],
                'runtime'       : config['runtime'],  # pilot runtime (min)
                'exit_on_error' : True,
                'project'       : config['project'],
                'queue'         : config['queue'],
                'cores'         : config['cores']
                }
        pdesc = rp.ComputePilotDescription(pd_init)

        # Launch the pilot.
        pilot = pmgr.submit_pilots(pdesc)


        report.header('submit units')

        # Register the ComputePilot in a UnitManager object.
        umgr = rp.UnitManager(session=session)
        umgr.add_pilots(pilot)

        # Create a workload of ComputeUnits. 
        # Each compute unit runs a MPI test application.
        report.info('create DEM unit description\n\t')


        # create a new CU description, and fill it.
        # Here we don't use dict initialization.
        cud = rp.ComputeUnitDescription()
        cud.environment    = ['PATH='+config['pathtoLIGGHTS']+':$PATH']
        cud.pre_exec       = ['mkdir CSVs','mkdir post','mkdir restart']
        cud.executable     = 'lmp_micstam'
        cud.arguments      = ['-in','in.*' ]
        cud.input_staging  = [{'source': config['pathtoLIGGHTSinputs']+'in.final_%fmm'%config['diameter'],
                               'target':'unit:///in.final_%fmm'%config['diameter'],
                               'action'  :rp.LINK},
                              {'source': config['pathtoLIGGHTSinputs']+'shell_closed.stl',
                               'target':'unit:///shell_closed.stl',
                               'action'  :rp.LINK},
                              {'source': config['pathtoLIGGHTSinputs']+'shell',
                               'target':'unit:///shell',
                               'action'  :rp.LINK},
                              {'source': config['pathtoLIGGHTSinputs']+'impeller',
                               'target':'unit:///impeller',
                               'action'  :rp.LINK},
                              {'source': config['pathtoLIGGHTSinputs']+'impeller_coarse.stl',
                               'target':'unit:///impeller_coarse.stl',
                               'action'  :rp.LINK},
                            ]
        cud.output_staging = [{'source': 'unit:///post/collision%d.atom'%config['timesteps'],
                               'target': 'pilot:///collision%d.atom'%config['timesteps'],
                               'action'  : rp.LINK},
                              {'source': 'unit:///post/impact%d.atom'%config['timesteps'],
                               'target': 'pilot:///impact%d.atom'%config['timesteps'],
                               'action'  : rp.LINK}]
        cud.cores          = config['DEMcores']
        cud.mpi            = True
        report.progress()
        report.ok('>>ok\n')

        # Submit the previously created ComputeUnit descriptions to the
        # PilotManager. This will trigger the selected scheduler to start
        # assigning ComputeUnits to the ComputePilots.
        units = umgr.submit_units(cud)

        # Wait for all compute units to reach a final state (DONE, CANCELED or FAILED).
        report.header('gather results')
        umgr.wait_units()
    
        report.info('\n')
        collision = {'source': 'pilot:///collision%d.atom'%config['timesteps'],
                     'target': 'unit:///sampledumpfiles/collision%d.%d_%f'%(config['timesteps'],config['cores'],config['diameter'])
                     'action'  : rp.LINK}

        impact = {'source': 'pilot:///impact200000.atom',
                  'target': 'unit:///sampledumpfiles/impact%d.%d_%f'%(config['timesteps'],config['cores'],config['diameter'])
                  'action'  : rp.LINK}
        # create a new CU description, and fill it.
        # Here we don't use dict initialization.
        cud = rp.ComputeUnitDescription()
        cud.environment    = ['PATH='+config['pathtoPBMexecutable']+':$PATH']
        cud.executable     = 'model.out'
        cud.arguments      = [config['DEMcores'],config['diameter']]
        cud.input_staging  = [collision,impact]
        cud.post_exec      = ['tar cfz csvDump.tar.gz csvDump']
        cud.output_staging = [{'source': 'unit:///csvDump.tar.gz',
                               'target': 'client:///csvDump.tar.gz',
                               'action'  : rp.TRANSFER}]

        cud.cores          = config['PBMcores']
        cud.mpi            = True
        report.progress()
        report.ok('>>ok\n')
        
        units = umgr.submit_units(cud)
        report.header('gather results')
        umgr.wait_units()


    except Exception as e:
        # Something unexpected happened in the pilot code above
        report.error('caught Exception: %s\n' % e)
        raise

    except (KeyboardInterrupt, SystemExit) as e:
        # the callback called sys.exit(), and we can here catch the
        # corresponding KeyboardInterrupt exception for shutdown.  We also catch
        # SystemExit (which gets raised if the main threads exits for some other
        # reason).
        report.warn('exit requested\n')

    finally:
        # always clean up the session, no matter if we caught an exception or
        # not.  This will kill all remaining pilots.
        report.header('finalize')
        session.close()

    report.header()


#-------------------------------------------------------------------------------

