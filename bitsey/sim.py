# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

'''
Big picture: what does the network look like?
    The network is just a bunch of cells and GJs. However, there is no "class
    cell" or "class GJ".
    Instead, our data structure is mostly just a bunch of arrays.
    Some (like the membrane voltage, Vm) are 1D arrays[N_CELLS], with 1 element
    per cell. Others (like the ion concentration, cc_cells) are 2D
    arrays[N_IONS,N_CELLS], and can thus track (e.g.,) the concentration of
    every ion in every cell independently. By keeping the information in arrays,
    we get to use canned matrix-math routines in numpy for efficient simulation.

Which arrays track the cells?
    The first set of arrays are all [N_IONS,N_CELLS]; they have one row per ion,
    and one column per cell.
	cc_cells[i,c] = current concentration of ion #i in cell #c. It is given
	  in moles/m3, which is equivalent to mmoles/liter. The elements of
	  cc_cells are the actual state variables of the simulation.
	Dm_array[i,c] = membrane ion-channel diffusion constant of cell #c for
	  ion #i (m2/s). One might think that Dm would just be the diffusion
	  constant of an ion itself, and would remain constant and be
	  mostly identical across all cells. Instead, Dm merges this with the
	  fraction of the cell's surface area that is covered by channels. So
	  if 1/3 of a cell membrane were covered by ion channels, then we would
	  cut that cell's Dm down by 3x. Also note that most ion channels use
          facilitated diffusion, and so these diffusion constants may be
          unrelated to an ion's "normal" diffusion constant in water.
    The next set are all arrays[N_CELLS]:
	z_array[i] = valence of ion #i. This array remains constant.
	cc_env [i] = concentration of ion #i in the ECF (moles/m3). In fact, it
	  remains constant also, since we assume the ECF is too big to be
	  affected by anything.
	Vm[i] = current Vmem of cell #i in volts. It is not, strictly speaking,
	  a state variable since we recompute it frequently from cc_cells and
	  z_array.

    ion_i is a dict that maps an ion name into its row index in the 2D
	arrays above; e.g., ion_i['K'] => 1.

What about gap junctions (GJs)?
    To get a network more complex than just isolated individual cells, we need
    GJs. These, too, are just stored in a big array:
	GJ_connects is a structured array [N_GJs]. Each element is one GJ.
          ['from']  and ['to'] specify the indices of its input and output cells
          (i.e., the column in the big arrays above).
    In addition to the network connectivity of the GJs, we must also set how
    quickly the various ions travel through them. GJ_diffusion[n_ions][n_GJs]
    gives the basic diffusion constant for each ion type through each GJ in
    m^2/s. The 'n_ions' index lets a GJ conduct different ions differently;
    and the n_GJs index gives us the flexibility of having different GJs each
    having a different diffusion constant for the same ion (perhaps because
    there are multiple GJ types in the sim, or perhaps because different cells
    have a different fraction of their cell-membrane area occupied by GJs).
    Init_big_arrays() creates each ion with (among its other parameters) its
    D_GJ. By default, D_GJ=(same as regular diffusion constant) for Na, K and
    Cl; and 0 for protein; and 1e-18 for any additional ions. Whichever the
    case, GJ_diffusion[ion,:] = that ion's D_GJ. However, you are welcome to
    change this value (e.g., to scale it up or down) anytime after
    init_big_arrays() is done.

The gating system
    The gating system is how we implement gated ion channels, gated GJs,
    generation and decay. We can gate based on either ligand concentration or
    Vmem. In fact, for instructional purposes, we allow gatings that are
    completely "un-physical." E.g., the concentration of a ligand in one cell
    can affect an ion-channel diffusion constant in a completely different cell!

    The gating system is easily extensible. There are a few built-in gates
    (e.g., a ligand-gated Hill function and a proportional metabolite decay),
    but you can add as many new ones as you like (e.g., Hodgkin-Huxley gating
    for ion channels in neurons).

    The heart of the system is Gate objects (i.e., objects of the Gate class).
    Each Gate object gates an entire vector of ICs, GJs or gen/decay ligands;
    we thus need only instantiate a small number of Gate objects in any sim. We
    may have one to do Hill-model gating of all Na channels, another to do the
    same for all K channels, and a third to implement M decay across all cells.

    This small number of Gates are kept in three lists: the ion-channel gates
    (reducing flux through ICs), the GJ gates (ditto, for GJs), and the
    gen/decay gates (which create an ion's generation or decay, potentially
    depending on other inputs to, e.g., build a GRN).

    Each Gate object contains:
    - the gating function. As noted, you can easily add new gating functions.
      Each gating function must be able to compute an entire vector in parallel.
    - a parameter set. E.g., a Hill function might have Km, N and KVmax. The
      format of this is up to the particular Gate subclass; usually it is a 2D
      array with one row per parameter and one column per cell or GJ (whichever
      is appropriate).

    As noted, the single gating function that a Gate object contains will
    compute the flux attenuation for an entire vector in parallel.
    Specifically, an IC or gen/decay function must deal with every cell at once;
    a GJ gating function must deal with every GJ at once. This vectoring is what
    makes the system run fast -- numpy vector operations are much faster than
    scalars. It is also more biologically realistic than having a Gate for each
    individual cell -- e.g., a given type of ion channel is typically present to
    some degree in each cell. However, it does mean that if you want, e.g., to
    gate Na channels in only a few cells, then you must be clever; typically by
    picking parameter values for other cells that have the effect of creating a
    gate of unity value that does nothing.

    The parameter set in a Gate object is thus actually a 2D array with one row
    per parameter, and one column per cell (or GJ). Thus, while all the cells in
    a single Gate object share the same gating function, they may have different
    parameters per cell. (And again, note that if you build your own Gate
    subclass, you can do whatever you like with the parameters -- but every one
    so far does follow this format).

    An end user typically builds a Gate object that gates an entire vector of
    cells or GJs) in one of two ways. Consider, e.g., a Hill gating. You first
    call make_Hill_gates() to make an entire gate all at once. It takes a
    parameter for the gate destination (IC, GJ or gen/decay), and numerous
    (in this case) Hill parameters. In most cases, this is the one-stop shopping
    function for making a gate. However, it makes an entire vector (i.e., the
    gating for every cell or GJ) all at once. You can give make_Hill_gates() a
    scalar for any given parameter (which assigns that value to every cell) or a
    vector (which can assign each cell its own value). However, the latter may
    not be convenient for the caller to set up, since you must assemble the
    entire 2D parameter vector in advance.

    The other option is to (as always) first call make_Hill_gates(), typically
    giving it parameters that create unity gating (i.e., a do-nothing gate).
    Then call make_Hill_gate() as many times as needed, with each call
    overriding the existing do-nothing parameters for one cell.

    What about the output of a Gate.func()? At each sim timepoint, sim_slopes()
    calls Gate.func() for each gate object. Ion-channel Gates are used in
        f_GHK[g.out_ion] *= g.func (cc_cells, Vm, deltaV_GJ, t)
    Typically, g.func() returns a vector[num_cells], where each element
    of the vector scales the IC flux for its respective cell for the single
    ion whose IC it gates.

    Gen/decay gates are used in
        slew_cc[g.out_ion] += g.func (cc_cells, Vm, deltaV_GJ, t)
    Typically, g.func() returns a vector[num_cells], where each element
    of the vector adds a flux in moles/(m3*s) to its respective cell.

    Finally, GJ gates are used in
        f_GJ *= g.func (cc_cells, Vm, deltaV_GJ, t)
    Typically, g.func() returns a vector[num_GJs], where each element of the
    vector scales the GJ flux for its respective GJ -- for every single ion.

    However, any given Gate subclass can return whatever shape it wants to,
    which can make the above code work with the magic of broadcasting.
    For example, an IC gate could actually return a vector[NUM_IONS,NUM_CELLS]
    that would scale *all* ion-channel fluxes in one fell swoop (regardless of
    how biologically realistic that may be).
    A GJ Gate could do the same, so as to scale different ion fluxes differently
    (though again, that seems unlikely to be biologically realistic).
    Or an IC Gate could return a single scalar, which would scale that Gate's
    IC in every cell equally.

    We can also play tricks with the inputs to g.func(), especially using
    advanced indexing. Sim_slopes() passes cc_cells[] and Vm[] to an IC
    g.func(). Then a line in g.func() might say
        conc = cc[self.in_ion],
    which would grab a vector[N_CELLS] of the concentrations for one ion in all
    cells, as a prelude to each cell computing a function of its own ion
    concentration to gate its own ion channel. But we can also do non-physical
    operations for instructional purposes.

    The trick is to use advanced indexing. If self.in_ion in this code is,not
    just a scalar, but the tuple (2,[0,3,2,1]), then the same code will pull
    ion #2 and then pick/rearrange the values to give each cell the conc for a
    *different* cell, so that (e.g.,) an ion channel in one cell gates based on
    a ligand in a *different* cell.

    A GJ Gate.func() not only gets cc_cells[] and Vm[], but also deltaV_GJ[]
    (which is an array[N_GJ] of the current deltaV across each GJ that
    sim_slopes() computes).
'''
import numpy as np
import scipy
import math	# To get pi.
import sim_toolbox as stb
import operator
import edebug as edb

#    import pdb; pdb.set_trace()

# Stuff all of the fundamental constants and commonly used parameters into a
# class. Any instance of this class will thus have all the constants.
class Params(object):
    def __init__(self):
        self.F = 96485  # Faraday constant [C/mol]
        self.R = 8.314  # Gas constant [J/K*mol]
        self.eo = 8.854e-12  # permittivity of free space [F/m]
        self.kb = 1.3806e-23  # Boltzmann constant [m2 kg/ s2 K1]
        self.q = 1.602e-19  # electron charge [C]
        self.tm = 7.5e-9  # thickness of cell membrane [nm]
        self.cm = 0.05  # patch capacitance of membrane [F/m2]
        self.T = 310  # temperature in Kelvin
        self.deltaGATP = -37000 # free energy released in ATP hydrolysis [J/mol]
        self.cATP = 1.5  # ATP concentration (mol/m3)
        self.cADP = 0.15  # ADP concentration (mol/m3)
        self.cPi = 0.15  # Pi concentration (mol/m3)
        self.alpha_NaK =1.0e-7 # max rate constant Na-K ATPase/unit surface area
        self.KmNK_Na = 12.0  # NaKATPase enzyme ext Na half-max sat value
        self.KmNK_K = 0.2  # NaKATPase enzyme ext K half-max sat value
        self.KmNK_ATP = 0.5  # NaKATPase enzyme ATP half-max sat value
        self.cell_r = 5.0e-6  # radius of single cell (m)
        self.GJ_len = 100e-9  # distance between two GJ connected cells (m)
        self.cell_sa = (4 * math.pi * self.cell_r ** 2)  # cell surface area
        self.cell_vol = ((4 / 3) * math.pi * self.cell_r ** 3)  # cell volume
        self.k26mV_inv = self.F / (self.R*self.T)

        # Simulation control. Sim_*_dump_interval is multiples of 
        # timestep, and is only relevant for explicit integration.
        self.time_step = .005		# min time step for explicit integration
        self.sim_dump_interval=10	# e.g., every 10 time steps
        self.sim_long_dump_interval=100
        self.no_dumps = False		# turn off all dumps if you want to

        # Parameters for adaptive explicit integration. These place a limit on
        # how much any cell's Vmem, or any ion concentration, can change in one
        # timestep.
        self.adaptive_timestep = True	# So that these params get used.
        self.sim_integ_max_delt_Vm = .0001	# Volts/step
        # .001 means that no [ion] can change by more than .1% in one timestep.
        self.sim_integ_max_delt_cc = .001
        self.use_implicit = False	# Use the implicit integrator
        # Parameters for implicit integration
        self.max_step = 1e6		# big number
        self.rel_tol = 1e-3		# these are the default rel_tol and
        self.abs_tol = 1e-6		# abs_tol for solve_ivp()

# Create the big simulation arrays: cc_cells, cc_env, Dm_array, z_array,
# GJ_diffusion, GJ_connects. And a few minor ones too.
def init_big_arrays (n_cells, n_GJs, p, extra_ions=[]):
    global cc_cells, cc_env, Dm_array, z_array, ion_i, \
           GJ_diffusion, GJ_connects, GP, pump_scaling
    GP = p

    ### TEMPORARY CODE FOR PUMPS?
    pump_scaling = np.ones ((n_cells))

    # ion properties (Name, base membrane diffusion [m2/s], valence
    #	initial concentration inside cell [mol/m3],
    #	fixed concentration outside cell [mol/m3], 
    # These are temporary structures. We use them to provide initial values for
    # the big arrays we are about to build, and to specify the order of which
    # row represents which ion in those arrays.
    Na={'Name':'Na', 'D_mem':1e-18, 'D_GJ':1e-18, 'z':1, 'c_in':10, 'c_out':145}
    K ={'Name':'K',  'D_mem':1e-18, 'D_GJ':1e-18, 'z':1, 'c_in':125,'c_out':5}
    Cl={'Name':'Cl', 'D_mem':1e-18, 'D_GJ':1e-18, 'z':-1,'c_in':55, 'c_out':140}
    P= {'Name':'P',  'D_mem':0,     'D_GJ':0,     'z':-1,'c_in':80, 'c_out':10}

    # stack the above individual dictionaries into a list to make it easier to
    # process them in the loop below.
    ions_vect = [Na, K, Cl, P]

    # Any particular sim may want to declare extra ions.
    for ion in extra_ions:
        ions_vect.append ({'Name':ion, 'D_mem':0.0, 'D_GJ':1e-18,
                           'z':0, 'c_in':0,  'c_out':0})
    n_ions = len(ions_vect)

    cc_cells = np.empty ((n_ions, n_cells))
    Dm_array = np.empty ((n_ions, n_cells))
    z_array  = np.empty ((n_ions))
    cc_env   = np.empty ((n_ions))
    GJ_diffusion = np.empty ((n_ions,n_GJs))

    ion_i    = {}

    # Push the parameters of the above ions into the various arrays.
    for row, ion_obj in enumerate(ions_vect):
        cc_cells[row,:] = ion_obj['c_in']	# initial cell conc
        cc_env  [row]   = ion_obj['c_out']	# fixed environmental conc
        Dm_array [row] = ion_obj['D_mem']	# initial membrane diff coeff
        z_array[row] = ion_obj['z']		# fixed ion valence
        GJ_diffusion[row] = ion_obj['D_GJ']	# diffusion rate through GJs
        ion_i[ion_obj['Name']] = row		# map ion name -> its row

    # Create default arrays for GJs.
    GJ_connects=np.zeros((n_GJs), dtype=[('from','i4'),('to','i4')])

    # The global lists of gating objects.
    global IC_gates, GJ_gates, GD_gates
    IC_gates=[]; GJ_gates=[]; GD_gates=[]

# The main "do-it" simulation function for both explicit and implicit 
# integration.
# It takes the current cc_cells[n_ions,n_cells], does all of the physics work,
# and returns an array of concentration slew rates [n_ions,n_cells]; i.e.,
# moles/m3 per second.
# Explicit simulation only calls sim_slopes() from sim.sim(); the implicit
# simulator calls if from sim.sim_implicit().
def sim_slopes (t):
    global cc_env, Dm_array, z_array, ion_i, Vm, GJ_connects, GP, \
           cc_cells, pump_scaling
    Vm = compute_Vm (cc_cells)

    # 1. We will return slew_cc as (moles/m3) per sec. However, we first
    #    accumulate most of the fluxes as (moles/m2) per sec, and then multiply
    #    by the cell surface area later on.
    # 2. For (moles/m2)/s: m2 of what area? You might think that for, e.g., ion
    #    channels, it should be per m2 of ion-channel area -- but it's not.
    #    *All* fluxes, including those through GJs, are per m2 of cell-membrane
    #    area; that's what makes the computation consistent for not only ICs,
    #    but also GJs and pumps.
    #    Thus, the (e.g.,) diffusion rate through ion channels must be scaled
    #    down by the fraction of membrane area occupied by channels (and similar
    #    for ion pumps and GJs).
    slew_cc = np.zeros (cc_cells.shape) # Temporarily in (moles/m2) per sec.

    # Run the Na/K-ATPase ion pump in each cell.
    # Returns two 1D arrays[N_CELLS] of fluxes; units are (moles/m2)/s
    pump_Na,pump_K,_ = stb.pumpNaKATP(cc_cells[ion_i['Na']],cc_env[ion_i['Na']],
                                      cc_cells[ion_i['K']], cc_env[ion_i['K']],
                                      Vm, GP.T, GP, pump_scaling)

    # Kill the Na-K pumps on some worm-interior cells. The hidden rule is that
    # if a cell has no ion channels (i.e., Dm=0 for all ions) then kill the
    # pumps for that cell.
    keep_pumps = np.any (Dm_array>0, 0)	# array[n_cells]
    pump_Na *= keep_pumps
    pump_K  *= keep_pumps

    # Accumulate the pump results into our running-total array slew-cc.
    slew_cc[ion_i['Na']] = pump_Na
    slew_cc[ion_i['K']]  = pump_K

    # deltaV across each GJ.
    deltaV_GJ = (Vm[GJ_connects['from']] - Vm[GJ_connects['to']]) # [n_GJs]

    # Ion-channel GHK fluxes. It returns an array of the same shape as cc_cells,
    # but typically only the Na, K and Cl rows are nonzero.
    f_GHK = stb.GHK (cc_env, cc_cells, Dm_array, z_array, Vm, GP)

    # Apply the IC gates; each Gate attenuates all cells for one ion.
    for g in IC_gates:
        f_GHK[g.out_ion] *= g.func (cc_cells, Vm, deltaV_GJ, t)
    slew_cc += f_GHK	# again, still in (moles/m2) per sec

    # Gap-junction computations.

    # Get the gap-junction Norton-equivalent circuits for all ions at once.
    # Units of Ith are mol/(m2*s); units of Gth are mol/(m2*s) per Volt.
    (GJ_Ith, GJ_Gth) = GJ_norton(deltaV_GJ)	# Both are [n_ions,n_GJs].
    f_GJ = GJ_Ith + deltaV_GJ*GJ_Gth	# [n_ions,n_GJs]

    # g.func returns scale[n_GJs]; f_GJ is [n_ions,n_GJs]
    for g in GJ_gates:
        f_GJ *= g.func (cc_cells, Vm, deltaV_GJ, t)

    # Update cells with GJ flux:
    # Note that the simple slew_cc[ion_index, GJ_connects['to']] += f_GJ
    # doesn't actually work in the case of two GJs driving the same 'to'
    # cell. Instead, we use np.add.at().
    for ion_index in range (cc_cells.shape[0]):	# for each ion
        np.add.at (slew_cc[ion_index,:], GJ_connects['from'], -f_GJ[ion_index])
        np.add.at (slew_cc[ion_index,:], GJ_connects['to'],    f_GJ[ion_index])

    # The current slew_cc units are moles/(m2*s), where the m2 is m2 of
    # cell-membrane area. To convert to moles/s entering the cell, we multiply
    # by the cell's surface area. Then, to convert to moles/m3 per s entering
    # the cell, we divide by the cell volume.
    slew_cc *= (GP.cell_sa / GP.cell_vol)

    # Next, do generation and decay. These are already natively in moles/(m3*s).
    for g in GD_gates:
        slew_cc[g.out_ion] += g.func (cc_cells, Vm, deltaV_GJ, t)

    return (slew_cc)	# Moles/m3 per second.

# Given: per-cell, per-ion charges in moles/m3.
# First: sum them per-cell, scaled by valence to get "signed-moles/m3"
# Next: multiply by F to convert moles->coulombs. Multiply by cell volume/
# surface area to get coulombs/m2, and finally divide by Farads/m2.
# The final scaling factor is F * p.cell_vol / (p.cell_sa*p.cm),
# or about 3200 mV per (mol/m3)
def compute_Vm (Cc_cells):
    global GP

    # Calculate Vmem from scratch via the charge in the cells.
    # rho_cells is [n_cells] and is coulombs/m3.
    rho_cells = (Cc_cells * z_array[:,np.newaxis]).sum(axis=0) * GP.F

    # Scale by vol/sa to get coulombs/m2-cell-surface-area, and then use V=Q/C.
    return (rho_cells * GP.cell_vol / (GP.cell_sa*GP.cm))

# The main entry point into this file. For implicit integration, it just
# farms the work to sim_implicit().
# Printout during simulation:
#	For explicit simulation, we print at intervals based on
#	sim_dump_interval and sim_long_dump_interval. For implicit simulation,
#	we don't print at all.
# Return values: (t_shots, cc_shots).
#	For both explicit and implicit simulation, we always take 100 snapshots
# during time [0,50] and then 200 shots over the entire rest of the simulation.
def sim (end_time):
    global cc_cells, Vm, GP, GJ_diffusion

    n_ions = cc_cells.shape[0]

    if (GP.use_implicit):
        return (sim_implicit (end_time))

    # Save snapshots of core variables for plotting.
    t_shots=[]; cc_shots=[]; last_shot=-100;

    # Sanity check the GJ sizes. With a concentration dC across a GJ, then the
    # resulting flux is Ithev=(D_GJ*/GJ_len)*dC moles/m2s, which gives a change
    # in concentration of (Ithev*cell_SA/cell_vol)*time_step. If this is so
    # negative that if flips dC in the opposite direction and increases its
    # magnitude, we're in trouble.
    for ion_index in range(n_ions):
        if (GP.time_step * GJ_diffusion[ion_index] * GP.cell_sa > .25 * GP.cell_vol*GP.GJ_len).any():
            print ("WARNING: your GJ diffusion rate may be too big for",
		"explicit integration at the current timestep for ion",ion_index)
            import pdb; pdb.set_trace()

    # Run the simulation loop:
    i=0; t=0
    last_print=-1000000			# Last time we printed a dump.
    last_long_print=-1000000		# Last time we long-printed a dump.
    time_step = GP.time_step		# seconds
    while (t < end_time):
        slew_cc = sim_slopes(t)

        if (GP.adaptive_timestep):
            # Compute Vmem slew (in Volts/s). Essentially, it's just slew_Q/C.
            # Slew_cc is slew-flux in moles/m3 per second.
            # First, in each cell, sum all ions weighted by valence.
            # Then mpy by...
            #   ...(cell_vol/cell_sa) ->moles/(m2 of cell-memb cross-sec area)/s
            #   ...F -> Coul /(m2 of cell-membrane cross-sec area)/s
            #   ...1/C_mem -> Volts/s.
            mult = (GP.cell_vol / GP.cell_sa) * (GP.F/ GP.cm)
            slew_Vm = (slew_cc * z_array[:,np.newaxis]).sum(axis=0) * mult

            # Timestep control.
            # max_volts / (volts/sec) => max_time
            max_t_Vm = GP.sim_integ_max_delt_Vm / (np.absolute (slew_Vm).max())
            # (moles/m3*sec) / (moles/m3) => fractional_change / sec
            frac_cc = slew_cc/(cc_cells+.00001)
            max_t_cc = GP.sim_integ_max_delt_cc / np.absolute(frac_cc).max()
            n_steps = max (1, int (min (max_t_Vm, max_t_cc) / time_step))
            #print ("S0L[1]=",cc_cells[6,1])
            #print ('At t={}: max_t_Vm={}s, max_t_cc={}s => {} steps'.format(t, max_t_Vm, max_t_cc, n_steps))
            #print ('steps_Vm=', (.001/(time_step*np.absolute (slew_Vm))).astype(int))
        else:
            n_steps = 1

        # Dump out status occasionally during the simulation.
        # Note that this may be irregular; numerical integration could, e.g.,
        # repeatedly do i += 7; so if sim_dump_interval=10 we would rarely dump!
        if ((i >= last_print + GP.sim_dump_interval) and not GP.no_dumps):
            long = (i >= last_long_print + GP.sim_long_dump_interval)
            if (long):
                last_long_print=i
            last_print = i
            edb.dump (t, cc_cells, edb.Units.mol_per_m3s, long)
            #edb.analyze_equiv_network (cc_cells, GP)
            #edb.dump_gating ()

        if (t>9999999):		# A hook to stop & debug during a sim.
            edb.debug_print_GJ (GP, cc_cells, 1)
            import pdb; pdb.set_trace()
            print (sim_slopes (t))

        cc_cells += slew_cc * n_steps * time_step
        i += n_steps
        t = i*time_step

        # Save information for plotting at sample points. Early on (when things
        # are changing quickly) save lots of info. Afterwards, save seldom so
        # as to save memory (say 100 points before & 200 after)
        boundary=min (50,end_time);
        before=boundary/100; after=(end_time-boundary)/200
        interval = (before if t<boundary else after)
        if (t > last_shot+interval):
            t_shots.append(t)
            cc_shots.append(cc_cells.copy())
            last_shot = t

    return (t_shots, cc_shots)

# Replacement for sim(); it uses scipy.integrate.solve_ivp()
# Like sim(), it returns (t_shots, cc_shots).
def sim_implicit (end_time):
    import scipy.integrate
    global cc_cells, Vm, GP
    num_ions, num_cells = cc_cells.shape

    def wrap (t, y):
        global cc_cells
        #print ('----------------\nt={:.9g}'.format(t))
        #np.set_printoptions(formatter={'float':'{:6.2f}'.format},linewidth=120)
        #print ('y={}'.format(y))
        #print ("Vm = ", compute_Vm (y.reshape((num_ions,num_cells))))
        cc_cells = y.reshape(num_ions,num_cells)       # So sim_slopes() sees it
        slew_cc = sim_slopes (t) 	# moles/(m3*s)
        slew_cc = slew_cc.reshape (num_ions*num_cells)
        #np.set_printoptions(formatter={'float':'{:7.2g}'.format},linewidth=120)
        #print ('slews={}'.format(slew_cc))
        return (slew_cc)

    # Save information for plotting at sample points. Early on (when things
    # are changing quickly) save lots of info. Afterwards, save seldom so
    # as to save memory. So, 100 points in t=[0,50], then 200 in [50, end_time].
    boundary=min (50,end_time)
    t_eval = np.linspace (0,boundary,50,endpoint=False)
    if (end_time>50):
        t_eval = np.append (t_eval, np.linspace (boundary, end_time, 2000))

    # run the simulation loop
    y0 = cc_cells.reshape (num_ions*num_cells)
    bunch = scipy.integrate.solve_ivp (wrap, (0,end_time), y0, method='BDF',
                                       t_eval=t_eval, max_step=GP.max_step,
				       rtol=GP.rel_tol, atol=GP.abs_tol)
    print ('{} func evals, status={} ({}), success={}'.format \
             (bunch.nfev, bunch.status, bunch.message, bunch.success))
    if (not bunch.success):
        raise ValueError
    t_shots = t_eval.tolist()
    # bunch.y is [n_ions*n_cells, n_timepoints]
    cc_shots = [y.reshape((num_ions,num_cells)) for y in bunch.y.T]
    cc_cells = cc_shots[-1]
    return (t_shots, cc_shots)

# Builds and returns a Norton equivalent model for all GJs.
# Specifically, two arrays GJ_Ith and GJ_Gth of [n_ions,n_GJ].
# - GJ_Ith[i,g] is the diffusive flux of ion #i in the direction of
#   GJ[g].from->to, and has units (mol/m2*s). It ignores ion valence completely,
#   since it's just diffusion.
# - GJ_Gth*(Vfrom-Vto) is the drift flux of particles (in mol/m2*s) in the
#   from->to direction. It takes into account the ion being positive or negative
#   (so as to predict the direction of the flux, which tells us which of the
#   two cell's [ion] to use), and also uses the ion's valence because drift flux
#   is proportional to valence.
#   GJ_Gth has units (mol/m2*s) per Volt. It is posi(nega)tive for
#   posi(nega)tive ions, so that the direction of the flux works out correctly.
# The caller uses them to compute total flux (in moles/(m2*s), in the direction
# from GJ.from to GJ.to, as flux = GJ_Ith + (Vfrom-Vto)*GJ_Gth.
def GJ_norton (deltaV_GJ):
    global cc_cells, GP, GJ_diffusion
    n_GJ = GJ_connects.size
    n_ions = cc_env.size

    # Compute ion drift and diffusion through GJs. Assume fixed GJ spacing
    # of GJ_len between connected cells.
    # First, compute d_conc/dx (assume constant conc in cells, and constant
    # gradients in the GJs).
    GJ_from = GJ_connects['from']	# Arrays of [n_GJ]
    GJ_to   = GJ_connects['to']
    deltaC_GJ=(cc_cells[:,GJ_to]-cc_cells[:,GJ_from])/GP.GJ_len #[n_ions,n_GJ]
    D = np.zeros ((n_ions, n_GJ))	# diffusion constants.
    conc = np.zeros ((n_ions, n_GJ))	# the concentration we use for drift

    for ion_index in range(n_ions):
        # Drift flux = velocity * conc. But a GJ connects two cells -- which
        # cell's concentration do we use? We originally used the average, but
        # that resulted in occasional negative concentrations. Now we use the
        # value of the cell that "sources" the ions (e.g., for a positive ion,
        # the cell with the more positive Vmem).
        # forw[n_GJ]: element #i is True if GJ[i] drifts forwards.
        forw = ((deltaV_GJ < 0) == (z_array[ion_index] < 0))
        # cell_idx[i] is the x such that GJ[i] uses cc[ion_index,x] for its conc
        cell_idx = np.where (forw, GJ_from, GJ_to)	# [n_GJ]
        # Now use advanced indexing to make conc[] for this ion, all GJs
        conc[ion_index] = cc_cells[ion_index,cell_idx]		# [n_GJ]

    # GJ drift flux.
    u = GJ_diffusion/(GP.kb*GP.T)	# drift mobility, from Einstein relation
    # [n_ions,n_GJ] * scalar * [n_ions] * [n_ions,n_GJ] = [n_ions,n_GJ]
    GJ_Gth = conc * GP.q * z_array.reshape((n_ions,1)) * u / GP.GJ_len
    GJ_Ith = -GJ_diffusion * deltaC_GJ	# diffusion in mol/m2s
    return (GJ_Ith, GJ_Gth)

# Global arrays of Gate objects
IC_gates=[]
GJ_gates=[]
GD_gates=[]

# The runtime data structure for a single Gate. Initialization parameters:
# - dest: says whether this gate controls an ion channel, GJ or gen/decay. It's
#   one of the three constants just above, and controls whether the gate gets
#   appended to IC_gates,GJ_gates or GD_gates.
# - out_ion: if the gate affects an IC or gen/decay, then which ion it affects.
#   However, a Gate that affects GJs will attenuate all ions equally, and hence
#   doesn't use this field (in which case we don't set it, and nobody should
#   look at it).
class Gate:
  # Constants saying which of the 3 global arrays above to put a new gate into.
  GATE_IC=0; GATE_GJ=1; GATE_GD=2

  def __init__(self, dest, out_ion):
    if (dest != self.GATE_GJ):	# Which ion does this Gate affect? (And see
        self.out_ion = out_ion	# the comment just above).
    self.width=(GJ_connects.size if (dest==self.GATE_GJ) else cc_cells.shape[1])
    gate_arrays = (IC_gates, GJ_gates, GD_gates)	# Stick this Gate into
    gate_arrays[dest].append (self)			# the appropriate list.

######################################################
# Hill gating; i.e., using Hill-function(ligand concentration) as a gate.
######################################################
# A Hill_gate object is a Gate that implements Hill-function gating; it 
# creates a new Hill gating for one output ion across all cells. It would most
# commonly be used for ion-channel gating or for generation/decay.
# You can supply kM, N and inv in many ways:
# - not at all, in which they default to unity gating; then fill the cells in
#   one by one later with change_Hill_gate(). Specifically, if you let 'inv'
#   default, then we get unity gates everywhere.
# - single numbers, which apply to all cells
# - a cell-by-cell array to gate each cell differently.
class Hill_gate (Gate):
  def __init__(self, gtype,out_ion,in_ion,inv=None,kM=1,N=1,kVMax=1,offset=0.0):
    Gate.__init__ (self, gtype, out_ion)
    self.in_ion=in_ion	# Tell Hill_gating() which ligand controls the gate.
    if (inv==None):	# default everything to unity gating
        inv=True; kM=1e30	# i.e., kVMax * (1 / (1 - 0**N)) with kVMax=1.
        inv=True; kM=1; N=0; kVMax=2 # alternate version, a bit more elegant.
    self.params = np.zeros ((5,self.width))	# (5 params) x (# cells or GJs)
    self.params[0],self.params[1],self.params[2],self.params[3],\
	self.params[4]=inv,kM,N,kVMax,offset

  # Changes the parameter of an existing gate 'g' in just one 'cell'.
  def change_Hill_gate (self, cell, inv, kM=1, kVMax=1, N=1, offset=0.0):
    self.params[0:5,cell] = [inv, kM, N, kVMax, offset]

  # If a ligand in one cell is to affect a gate in a *different* cell, then
  # Hill_gate_unphysical_inputs() uses advanced indexing to let us do so.
  # The trick: normally in_ion is simply the number of which ion gates us, and
  # Hill_gating() gets its vector of inputs with cc_cells[in_ion]. So if
  # in_ion=1, then cc_cells[1] is [K] for all cells.
  # However, Hill_gate_unphysical_inputs() makes in_ion into a tuple; e.g.,
  # (1, [2,3,3,0,0]) for a 5-cell system. Then cc_cells[(1, [2,3,3,0,0])] gives
  # us cc_cells[1][2,3,3,0,0], which first gets cc_cells[1] and then uses the
  # trick of indexing a vector with another vector to give us a *reordering* of
  # cc_cells[1]: in this case, it is [1,2], [1,3], [1,3], [1,0], [1,0].
  # Unphysical gating is *needed* anytime we want to gate GJs rather than ICs
  # (see the comment at the top of the file for more details).
  def Hill_gate_unphysical_inputs (self, inputs):
    self.in_ion=tuple([self.in_ion,inputs])

  # The function that actually computes the gating vector at runtime.
  def func (self, cc, Vm_ignore, deltaV_GJ_ignore, t_ignore):
    # self.params[] is [5, # cells or GJs]. Assign each of its rows (i.e., a
    # vector of one parameter across all cells) to one variable.
    inv,kM,N,kVMax,offset = self.params #Each var is a vector (one val per cell)
    conc = cc[self.in_ion]	# vector [N_CELLS]
    out = 1 / (1 + ((conc/kM)**N))
    out = kVMax*((inv*out + (1-inv)*(1-out)) + offset)	# any better way to do self?
    return (out)


######################################################
# Vmem gating; i.e., use Vmem to gate things via the Hill-like function
#	kVMax / (1 + np.exp (N * (V-kM)))
######################################################
# The general-usage comments above for Hill-function gating still apply here.

# Force the attenuation to 1 with kM=big and kVMax=1
class Vmem_gate (Gate):
  def __init__ (self, gtype, out_ion, kM=100, N=1, kVMax=1):
    Gate.__init__ (self, gtype, out_ion)
    self.in_cells = []
    self.params = np.zeros ((3,self.width))
    self.params[0], self.params[1], self.params[2] = kM, N, kVMax

  def change_Vmem_gate (self, cell, kM, N=1, kVMax=1):
    self.params[0:3,cell] = kM, N, kVMax

  def Vmem_gate_unphysical_inputs (self, inputs):
    self.in_cells = inputs

  def func (self, cc_ignore, Vm, deltaV_GJ_ignore, t_ignore):
    kM,N,kVMax = self.params	# assign each row to one variable
    V = Vm[self.in_cells]
    out = kVMax / (1 + np.exp (N * (V-kM)))
    return (out)

######################################################
# Const-gen/decay gating. Implement generation or decay of one ion.
# Generation is a simple constant generation rate in moles/(m^3 * s);
# Decay is simple linear decay (constant decay rate in 1/s).
# As usual, the constants can be per-Gate or per-cell. But if you want
# Hill-function generation, then this Gate isn't right for you. You could use
# a Hill-function gate with gtype=GATE_GD; then instantiate a GD_const_gating
# with gen=0 to get the decay.
######################################################
# The general-usage comments above for Hill-function gating still apply here.

class GD_const_gate (Gate):
  def __init__ (self, gtype, out_ion, gen=0, decay=0):
    assert (gtype==self.GATE_GD)	# There's no other sensible choice.
    Gate.__init__ (self, dest=gtype, out_ion=out_ion)
    self.in_ion=out_ion
    self.params = np.zeros ((2,self.width))
    self.params[0], self.params[1] = gen, decay

    # Sanity check: if k_d is too big, then explicit integration can oscillate
    # wildly. For an ion concentration I in moles/m3, I_next = I - I*k_d*dt
    # = I*(1-k_d*dt). This can explode if |1-k_d*dt| > 1. Since k_d is always>0,
    # it explodes if k_d*dt > 2.
    global GP
    if (decay * GP.time_step > 2) or (decay<0):
        print ("WARNING: your degradation k_d may be too big for explicit",
		"integration at the current timestep")

  def change_GD_const_gate (self, cell, gen=0, decay=0):
    self.params[0:2, cell] = [gen, decay]

  def func (self, cc, Vm_ignore, deltaV_GJ_ignore, t_ignore):
    # Units are moles/m3.
    gen_rate,decay_rate = self.params	# Each is vector[N_CELLS]
    out = gen_rate - cc[self.in_ion]*decay_rate
    return (out)
