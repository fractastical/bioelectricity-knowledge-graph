#!/usr/bin/env python3

# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

#    import pdb; pdb.set_trace()

import numpy as np
import sim as sim
import eplot as eplt
import edebug as edb
import re

# A nice central place to set the three tau values we will use.
M_TAU=1
H_TAU=10
N_TAU=10

######################################################
# For an overview of the network data structures, see sim.py

######################################################
# HH gating
######################################################

# The big picture:
# We gate the Na and K ICs with the usual HH model; Na via m^3*h and K via n^4.
# This is the functions HH_Na_gating() and HH_K_gating().
# But where do m, h and n come from?

# We'll have three ligands m, h and n. They are not "real" ligands; they exist
# only to calculate how to gate the Na and K ion channels. Thus, they all have
# Z=0 and they don't diffuse through ICs or GJs (so they stay in whatever cell
# they started in).
# Each cell creates its m, h and n via gen/decay. Like any gen/decay gate, we
# thus have gating functions -- HH_M_gen_decay(), HH_H_gen_decay() and
# HH_N_gen_decay() -- that return the net generation of their species at any
# point in time. They do this calculation by the classic method of first
# calculating (e.g., for M) M_tau and M_inf and then returning (M_inf-M)/M_tau.
# However, instead of taking the published equations for M_inf and M_tau, I
# tweaked my own (I couldn't get the standard ones to work).

# Finally, we implement a little Na-generation gate that kicks off our one and
# only AP at about t=200.

# The HH gating functions that gate the Na and K ion channels.
# The Na IC is gated by m^3 * h; the K IC by n^4. Since we have m, h and n
# simply being ligands, then Na and K are gated with simple ligand-gated ion
# channels. Here are the gating functions (the ones that are actually called
# at runtime).
# We return an array[num_cells] of numbers in [0,1] -- the gating numbers.
class HH_Na_gating (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_IC, out_ion=sim.ion_i['Na'])

  def func (self, cc, Vm_ignore, deltaV_GJ_ignore, t_ignore):
    m=sim.ion_i['m']; h=sim.ion_i['h']
    gate = np.power(cc[m,:], 3.0) * cc[h,:]	# m^3 * h
    # In principle, m and h are always in [0,1]. In practice, if m_tau is very
    # small then the numerical integrator can set [m] slightly <0. But since
    # m^3*h scales the GHK fluxes, then a negative [m] can *reverse* the flux
    # direction of Na, and suddenly our nice Na flux that GHK uses as a
    # stabilizing force instead makes the system run away :-(.
    return (np.clip (gate, 0, 1))

class HH_K_gating (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_IC, out_ion=sim.ion_i['K'])

  def func (self, cc, Vm_ignore, deltaV_GJ_ignore, t_ignore):
    n=sim.ion_i['n']
    return (np.power (cc[n,:],4.0))		# n^4

# Next, the HH gen/decay gating functions that control the concentrations of
# m, h and n. There are two steps for each:
# 1. Create (e.g., for m) m_inf and m_tau as functions of Vmem. They are
#    typically given by mathematical functions. Instead, we build them here with
#    simple linear interpolation.
# 2. dm/dt = (m_inf-m)/tau_m. And dm/dt is the net-generation rate that our
#    gen/decay gate must return.
class HH_M_gen_decay (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_GD, out_ion=sim.ion_i['m'])
  def func (self, cc, Vm, deltaV_GJ_ignore, t_ignore):
    M = cc[sim.ion_i['m'],:]
    M_inf = np.interp (Vm, [-1,-.050, .010, 1], [.1,.1,1,1])
    M_tau = np.ones_like (Vm)*M_TAU
    #print ("M_inf={}, tau={:.6g}".format(M_inf, M_tau[0]))
    return ((M_inf - M) / M_tau)

class HH_H_gen_decay (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_GD, out_ion=sim.ion_i['h'])
  def func (self, cc, Vm, deltaV_GJ_ignore, t_ignore):
    H = cc[sim.ion_i['h'],:]
    H_inf = np.interp (Vm, [-1,-.050, -.045, 1], [1,1,.001,.001])
    H_tau = np.ones_like (Vm)*H_TAU
    return ((H_inf - H) / H_tau)

class HH_N_gen_decay (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_GD, out_ion=sim.ion_i['n'])
  def func (self, cc, Vm, deltaV_GJ_ignore, t_ignore):
    N = cc[sim.ion_i['n'],:]
    N_inf = np.interp (Vm, [-1,-.050, .010, 1], [.1,.1,1,1])
    N_tau = np.ones_like (Vm)*N_TAU
    return ((N_inf - N) / N_tau)

# The gating function for Na gen/decay. We use it to slowly ramp up Vmem so as
# to trigger an AP at about t=200.
class HH_delayed_gen (sim.Gate):
  def __init__ (self):
    sim.Gate.__init__ (self, dest=self.GATE_GD, out_ion=sim.ion_i['Na'])

  def func (self, cc, Vm, deltaV_GJ_ignore, t):
    arr = np.zeros_like(Vm)
    if ((t > 200) and (t<205)):
        arr[0] = .03
    return (arr)

# The HH entry point.
def setup_HH (p):
    num_cells=1
    n_GJs = 0
    p.use_implicit = True	# Use the implicit integrator
    # Max_step ensures that the implicit integrator doesn't completely miss the
    # skinny input pulse of [Na] that kicks off the AP. Rel_tol and abs_tol get
    # tightened because otherwise we get nonphysical results such as Vmem
    # dropping down to -250mV even when the lowest Vnernst is -70mV.
    # Also note that the AP is quite skinny -- so if you run the sim for too
    # long, the plotting granularity becomes too wide to catch the AP!
    p.max_step=2; p.rel_tol=1e-5; p.abs_tol=1e-8

    # Each cell gets 3 new 'species' m, h and n. They are neutral and stay in
    # their original cells.
    sim.init_big_arrays (num_cells, n_GJs, p, ['m', 'h', 'n'])
    m=sim.ion_i['m']; h=sim.ion_i['h']; n=sim.ion_i['n']
    Na=sim.ion_i['Na']; K=sim.ion_i['K']; Cl=sim.ion_i['Cl']; P=sim.ion_i['P']
    sim.z_array[[m,h,n]] = 0
    sim.GJ_diffusion[[m,h,n]] = 0	# No diffusion through GJs or ICs, so
    sim.Dm_array[[m,h,n],:]= 0		# m, n, & h just stay in their cells.

    # This test case was built around D_Na = 1e-18 and D_K = 20e-18 (which
    # makes our initial cell concentrations pretty stable). But that was before
    # we added our HH gating into the Na and K ICs.
    # At steady-state baseline (i.e., waiting for an AP), we have m^3h=1e-3
    # and n^4=1e-4. So change D_Na and D_K to make them right at baseline.
    sim.Dm_array[Na,:] = 1.0e-15
    sim.Dm_array[K,:] = 20.0e-14

    # Initial concentrations. Our sim will settle to whatever is stable. We
    # speed that up by seeding the final settled values from a previous sim.
    # This seeding means that we can kick off the AP nice and early, without
    # needing to wait a long time for things to stabilize first. Perhaps more
    # importantly, it means that we don't start off the sim at a Vmem that
    # would kick off an AP.
    sim.cc_cells[K, :] = 87.572
    sim.cc_cells[Cl,:] = 15.669
    sim.cc_cells[P,:] = 79.652

    # Our desired Vmem is -58.5mV. Pick Na by first setting it for initial
    # charge neutrality...
    sim.cc_cells[Na,:]  = sim.cc_cells[Cl,:]+sim.cc_cells[P,:]-sim.cc_cells[K,:]
    assert ((sim.cc_cells[Na,:]>0).all())

    # ... then offset [Na] a bit to get to -58.5mV initially.
    scale = (p.cell_sa*p.cm) / (p.F * p.cell_vol)
    sim.cc_cells[Na,:] = sim.cc_cells[Na,:] -.0585*scale

    # Initial baseline values of m,h,n. Again, these make sure that we start
    # the sim at baseline rather than in the middle of an AP.
    sim.cc_cells[m,:] = .1
    sim.cc_cells[h,:] = 1
    sim.cc_cells[n,:] = .1

    # Instantiate the gen/decay gates that generate [m,h,n] with a rate of
    # (e.g., for m) dm/dt=(m_inf-m)/m_tau.
    HH_M_gen_decay()
    HH_H_gen_decay()
    HH_N_gen_decay()

    # The IC gates: gate D_Na by m^3*h and D_K by n^4.
    HH_Na_gating()
    HH_K_gating()

    # Create a trickle of Na+ into the cell via generation to kick off the AP
    HH_delayed_gen()

def post_HH (GP, t_shots, cc_shots):
    eplt.plot_Vmem (t_shots, cc_shots, title=f"{M_TAU=}, {H_TAU=}, {N_TAU=}")

    # Look at m, h and n for debugging.
    eplt.plot_onecell (t_shots, cc_shots, 0, ['m','h','n'],
                       title=f"m,h,n for {M_TAU=}, {H_TAU=}, {N_TAU=}")

    # Look at m^3*h and n^4, the final gates on the Na & K ICs.
    Na=sim.ion_i['Na']; K=sim.ion_i['K']
    m=sim.ion_i['m']; n=sim.ion_i['n']; h=sim.ion_i['h']
    Gna=[cc[m]**3 * cc[h]* sim.Dm_array[Na] for cc in cc_shots]
    Gk =[cc[n]**4 *sim.Dm_array[K, :] for cc in cc_shots]
    eplt.plot_onecell (t_shots, cc_shots, 0, [(Gna,'NaGate'), (Gk,'KGate')],
                       title=f"m^3h and n^4 for {M_TAU=}, {H_TAU=}, {N_TAU=}")

    # Get end-of-sim Vnernst values as a sanity check.
    # edb.analyze_equiv_network (cc_shots[-1],GP)

    # In principle, the APs should be so fast the [Na], [K] and [Cl] don't
    # change much. This plot shows that actually, they do change a bit.
    # eplt.plot_onecell (t_shots, cc_shots, 0, 'Na','K', 'Cl')

######################################################
# End of Hodgkin Huxley.
######################################################

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
