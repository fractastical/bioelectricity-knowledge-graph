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

# This is the function we're using for simulation-lab #1 (SS bioelectricity)
# for the first simulation.
def setup_SS1 (p):
    p.adaptive_timestep = True		# Fast (but sometimes unstable) mode
    p.sim_dump_interval = 1000000	# Every 1M*.005=5K seconds

    num_cells=4		# 4 cells

    # Gap junctions connect cells (we'll learn more soon).
    # No gap junctions means that the 4 cells are independent.
    n_GJs = 0

    # We've declared how many cells and GJs the simulation has. Create them
    # all, with default values.
    sim.init_big_arrays (num_cells, n_GJs, p)

    # To make it easier to set per-ion values below.
    Na = sim.ion_i['Na']; K = sim.ion_i['K']; Cl = sim.ion_i['Cl']

    # We want this lab to use a little different default values than the usual
    # ones. Specifically, a different diffusion constant for the K ion channels.
    # Note that this is the baseline value; simulation #3 will double it.
    # The notation [K,:] means to set the K ion channel in all four cells.
    sim.Dm_array[K,:] = 10.0e-18

def setup_SS2 (p):
    p.adaptive_timestep = True		# Fast (but sometimes unstable) mode
    p.sim_dump_interval = 1000000	# Every 1M*.005=5K seconds
    num_cells=4				# 4 cells
    n_GJs = 0				# No gap junctions.
    sim.init_big_arrays (num_cells, n_GJs, p)	# Create the network.

    # To make it easier to set per-ion values below.
    Na = sim.ion_i['Na']; K = sim.ion_i['K']; Cl = sim.ion_i['Cl']

    # We want this lab to use a little different default values than the usual
    # ones. Specifically, a different diffusion constant for the K ion channels.
    # Note that this is the baseline value; simulation #3 will double it.
    sim.Dm_array[K,:] = 10.0e-18

    # Simulation #2 asks you to change the initial cell-internal concentrations.
    # Here's where you can do that.
    # The defaults (set in sim.init_big_arrays) are [Na+]=10, [K+]=125,
    # [Cl-]=55, [other-]=80 (units are all moles/m^3). Note that they add up to
    # charge neutral. Make sure to keep the cell interior pretty close to charge
    # neutral, or else the simulation will explode!
    # For example, if we want to double [Na+], it would go from 10 to 20. We
    # will also increase [Cl] by 10 moles/m^3 to preserve charge neutrality.
    # You should do similar things for the other cells, as per the instructions.
    sim.cc_cells[Na,1] += 10	# Cell #1: [Na] is 20 rather than 10
    sim.cc_cells[Cl,1] += 10	# [Cl] changes to keep charge neutrality.
    #... do something similar for cells #2 and #3.

def setup_SS3 (p):
    p.adaptive_timestep = True		# Fast (but sometimes unstable) mode
    p.sim_dump_interval = 1000000	# Every 1M*.005=5K seconds
    num_cells=4				# 4 cells
    n_GJs = 0				# No gap junctions.
    sim.init_big_arrays (num_cells, n_GJs, p)	# Create the network.

    # To make it easier to set per-ion values below.
    Na = sim.ion_i['Na']; K = sim.ion_i['K']; Cl = sim.ion_i['Cl']

    # We want this lab to use a little different default values than the usual
    # ones. Specifically, a different diffusion constant for the K ion channels.
    # Note that this is the baseline value; simulation #3 will double it.
    sim.Dm_array[K,:] = 10.0e-18

    # Simulation #3 asks you to alter the ion-channel conductances D_Na, D_K
    # and D_Cl. By default, sim.init_big_arrays() sets all three to 1e-18, and
    # then we overrode D_K above to 10e-18. Here's where we override them for
    # simulation #3.
    sim.Dm_array[Na,1] *= 2	 # Cell #1 doubles D_Na
    #... do something similar for cells #2 and #3.

# "Post_" functions are the standard way in BITSEY to produce plots or other
# information from post-simulation data.
# GP:	an object with lots of basic parameters such as cell size and
#	physical constants (you can find the details in sim.py).
# t_shots,cc_shots: BITSEY doesn't call post_lab1() until *after* the simulation
#	is over. So the only to plot, say, a cell voltage *over time* is if the
#	simulation has saved lots of data at lots of timepoints during the sim.
#	That's what t_shots and cc_shots are. T_shots is a list of which
#	timepoints we have saved data for. Cc_shots is also a list, with one
#	item for each timepoint in t_shots. Each of these items is the 2D array
#	of the concentration of every ion in every cell, which is enough data to
#	compute all Vmem values (and indeed most anything else we would want
#	to plot).
def post_SS1 (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots)
    #eplt.plot_ion (t_shots, cc_shots, 'Na', title="Na in all cells")

def post_SS2 (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots)
    eplt.plot_ion (t_shots, cc_shots, 'Na', title="Na in all cells")
    eplt.plot_ion (t_shots, cc_shots, 'K',  title="K in all cells")
    eplt.plot_ion (t_shots, cc_shots, 'Cl', title="Cl in all cells")

def post_SS3 (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots)
    #eplt.plot_ion (t_shots, cc_shots, 'Na', title="Na in all cells")


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
