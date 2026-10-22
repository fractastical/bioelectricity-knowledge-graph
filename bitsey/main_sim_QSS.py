#!/usr/bin/env python3

# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

#    import pdb; pdb.set_trace()

import numpy as np
import sim as sim
import eplot as eplt
import edebug as edb
import re

#####################################################
# For an overview of the network data structures, see sim.py

# This is the function we're using for lab #1
def setup_QSS1 (p):
    p.adaptive_timestep = False		# Slow but accurate mode
    p.sim_dump_interval = 10		# Because the lab asks for long-dump
    p.sim_long_dump_interval = 10	# info at t=.05 sec
    num_cells=4
    n_GJs = 0
    sim.init_big_arrays (num_cells, n_GJs, p)

    # By default, we have D_Na = D_K = 1e-18
    Na = sim.ion_i['Na']; K = sim.ion_i['K']
    sim.Dm_array[Na,0]= 20.0e-18
    sim.Dm_array[Na,1]=  5.0e-18
    sim.Dm_array[K,2] = 10.0e-18
    sim.Dm_array[K,3] = 160.0e-18

def post_QSS1 (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots)

# This function shows that the QSS Vmem is insensitive to initial Vmem.
# 4 cells. Each one has Dm_Na set so that QSS Vmem will be -57mV.
# Tweak initial [Na] in cells 1-3 so that each starts out 5mV higher than
# the previous one. Nonetheless, they all converge to the same 30mV very
# quickly.
def setup_QSS2 (p):
    p.adaptive_timestep = False		# Slow but accurate mode
    p.sim_dump_interval = 10		# Because the lab asks for long-dump
    p.sim_long_dump_interval = 10	# info at t=.05 sec
    num_cells=4
    n_GJs = 0
    sim.init_big_arrays (num_cells, n_GJs, p)

    # By default, we have D_Na = D_K = 1e-18
    Na = sim.ion_i['Na']; K = sim.ion_i['K']
    sim.Dm_array[Na,:] = 1e-18	# Cell #2 Na isn't set, so use the default
    sim.Dm_array[K,:] = 0	# Put in the actual number here!

    # Offset the initial concentrations.
    sim.cc_cells[Na,0] -= .025
    sim.cc_cells[Na,1] -= 0	# Put in the actual number here!
    sim.cc_cells[Na,2] -= 0	# Put in the actual number here!
    sim.cc_cells[Na,3] -= 0	# Put in the actual number here!

def post_QSS2 (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots)


######################################################################
######################################################################

def command_line ():
    # If no command-line arguments are given, prompt for them.
    import sys
    if (len (sys.argv) <= 1):
        args = 'main.py ' + input('Arguments: ')
        sys.argv = re.split (' ', args)

    if (len(sys.argv) != 3):
       raise SyntaxError('Usage: python3 main.py test-name-to-run sim-end-time')

    end_time = float (sys.argv[2])
    name = sys.argv[1]

    # Run whatever test got chosen on the command line.
    GP = sim.Params()
    eval ('setup_'+name+'(GP)')

    # Initialize Vmem -- typically 0, or close to that.
    Vm = sim.compute_Vm (sim.cc_cells)
    assert (np.abs(Vm)<.5).all(), \
            "Your initial voltages are too far from charge neutral"
    np.set_printoptions (formatter={'float': '{: 6.3g}'.format}, linewidth=90)
    print ('Initial Vm   ={}mV\n'.format(1000*Vm))

    print ('Starting main simulation loop')
    t_shots, cc_shots = sim.sim (end_time)
    print ('Simulation is finished.')

    # We often want a printed dump of the final simulation results.
    np.set_printoptions (formatter={'float': '{:.6g}'.format}, linewidth=90)
    edb.dump (end_time, sim.cc_cells, edb.Units.mol_per_m3s, True)
    #edb.dump (end_time, sim.cc_cells, edb.Units.mV_per_s, True)

    # Run a post-simulating hook for plotting, if the user has made one.
    import sys
    if ('post_'+name in dir(sys.modules[__name__])):
        eval ('post_'+name+'(GP, t_shots, cc_shots)')

command_line()
