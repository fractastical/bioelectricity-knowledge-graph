# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

import numpy as np
import matplotlib	# These two lines are sometimes needed to plot
#matplotlib.use('Agg')	# to a file from a process with no X display.
import matplotlib.pyplot as plt
import re
from matplotlib.collections import LineCollection
import sim
#import pdb; pdb.set_trace()

# Various simple functions to plot voltage, ion concentration, etc.
#   plot_Vmem (t_shots, cc_shots, cells=None, title=None, filename=None)
#   plot_ion  (t_shots, cc_shots, ions, cells=None, filename=None)
#	Plot voltage or ion concentration vs. time in one or more cells (one
#	cell per plot-line).
#	'cells' is a list of which cells to plot (if not given, plot all cells).
#	'ions' is a single string (e.g., 'Na') or a list of strings. It can also
#		be a mixed list such as ['Na[10]','Na[13]','K[5]'] for the
#		ultimate in flexibility; i.e., any ion that is given with a
#		subscript will use that rather than 'cells'.
#   plot_onecell (t_shots, cc_shots, cell, what_to_plot, filename=None)
#	Plot Vmem *and* ion concentrations and other optional quantities in
#	just one cell, vs. time.
#	Vmem always gets plotted. 'What_to_plot' says what else to plot; it is
#	an array of specifiers. Most simply, each is just an ion name (e.g.,
#	'Na'). So to plot Vmem, [Na] and [Cl] in cell #4 you would call
#		plot_onecell (t_shots, cc_shots, 4, ['Na','Cl'])
#	You can also get fancy -- any specifier in what_to_plot can also be a
#	tuple of ([data_vals, label]). So you would prepare, on your own, an
#	array data_vals (which is a list of numbers to plot, and for which each
#	number in the list corresponds to one of the times in t_shots), and
#	'label' is what to call the data in this tuple.
#	In any case, the horizontal axis is time. One vertical axis for Vmem,
#	and another for everything else.
#   plot_worm (cc_cells, what_to_plot, title=None, filename=None):
#	Plot one quantity (Vmem or a single [ion]) across every cell in a worm.
#	'What_to_plot' can be 'V' or, e.g., 'Na', or even ['Na','K'].
#	We only plot one single timepoint -- the x axis is cell number, not
#	time.
#   plot_worm_times (ion, t_shots, cc_shots, n_views,
#			t_start=0, t_end=1, title=None, filename=None):
#	Very similar to plot_worm(), but plots multiple snapshots in time on
#	the same graph. N_views, t_start and t_end work together to pick the
#	timepoints; if the sim went to t=100, then (4,0,1) would plot at
#	t=0,25,75,100; (3,.3,.5) would plot at t=30,40,50.
#	This function can only plot a single ion concentration, not Vmem.

def plot_Vmem (t_shots, cc_shots, cell_idxs=None, title=None, filename=None):
    # If cell_idxs isn't given, make it an array of range(n_cells).
    if (cell_idxs==None):
        cell_idxs = range(cc_shots[0].shape[1])
    if (title==None):
        title = f'Vmem in cell{list(cell_idxs)}'

    Vm_shots = [sim.compute_Vm (c)*1000 for c in cc_shots]
    for cell in cell_idxs:
        data = [shot[cell] for shot in Vm_shots]
        plt.plot (t_shots, data, label=f'Vmem[{cell}]')
    plot_common('Vmem(mV)', title, filename)

def plot_ion (t_shots,cc_shots, ions,cells=None, title=None, filename=None):
    if (not type(ions) == list):	# If 'ions' is a scalar, turn it into a
        ions = [ions]			# single-element list for consistency.
    # If cells isn't given, make it an array of range(n_cells).
    if (cells==None):
        cells = range(cc_shots[0].shape[1])

    for ion in ions:
        match = re.match (r'([^[]+)\[([0-9]+)\]$', ion)
        if (match):
            ion=match.group(1)
            idxs = [int(match.group(2))]
        else:
            idxs = cells
        ion_i = ion_idx (ion)
        for cell_idx in idxs:
            data = [cc[ion_i,cell_idx] for cc in cc_shots]
            plt.plot(t_shots, data, linewidth=2.0, label=f'{ion}[{cell_idx}]')

    plot_common ("mol/m3", title, filename)

# Plot the x and y labels and the legend (i.e., the list of what things we're
# graphing and which color line each uses).
# Then call plot_it() to output the plot.
def plot_common (ylabel, title=None, filename=None):
    plt.xlabel('Time (s)', fontsize = 20)
    plt.ylabel(ylabel)

    # Now the legend.
    leg = plt.legend(loc='best',ncol=2,mode="expand",shadow=True,fancybox=True)
    leg.get_frame().set_alpha(0.5)

    plot_it(title, filename)

# Plot lots of things in just one cell. 
def plot_onecell (t_shots,cc_shots, cell,what_to_plot,title=None,filename=None):
    fig,ax1 = plt.subplots()	# First make an axis *just for Vmem*
    ax1.set_xlabel('Time (s)', fontsize = 20)
    ax1.set_ylabel('Vmem (mV)')

    # Plot Vmem in the axis we just made. It's always a black dotted line.
    Vm_shots=np.asarray([sim.compute_Vm (c)[cell]*1000 for c in cc_shots])
    ax1.plot(t_shots, Vm_shots, ':k', linewidth=2.0, label="Vmem")
    print ('Vmem is in [',Vm_shots.min(),':',Vm_shots.max(),']')

    # Now a second axis for everything else. It's mostly for ion concentrations;
    # so label it as mM/m^3, even if we plot other stuff on it too.
    ax2 = ax1.twinx()
    ax2.set_ylabel('mM/m^3')

    for what in what_to_plot:
        if (isinstance(what,str)):	# A string is an ion name
            ion = ion_idx (what)
            data = np.asarray([shot[ion,cell] for shot in cc_shots])
            print (what,'for cell #',cell,'is in [',data.min(),
							':',data.max(),']')
            ax2.plot(t_shots, data, linewidth=2.0, label='['+what+']')
        elif (type(what) is tuple):	# (data array to plot, label)
            plt.plot (t_shots, what[0], label=what[1])
        else:
            raise ValueError("Unknown argument",what,"to eplot.plot_onecell()")

    # Now for a nuisance. We made two axes -- and by default each has its own
    # legend (the little diagram listing what we graphed). So merge the two
    # into one single legend.
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    lines = lines_1+lines_2; labels = labels_1+labels_2
    legend = ax1.legend(lines, labels, loc='best',ncol=2,shadow=True)
    legend.get_frame().set_alpha(0.5)

    plot_it(title, filename)

######################################################
# Plot Vmem or a set of [ion] across every cell in a worm. I.e., what_to_plot
# can be 'V' or something like ['M', 'P']
# We only plot one single timepoint -- the x axis is cell number, not time.
######################################################
def plot_worm (cc_cells, what_to_plot, title=None, filename=None):
    plt.figure()

    if (what_to_plot == 'V'):
        data = sim.compute_Vm (cc_cells)
        plt.plot(data*1000, linewidth=2.0, marker='.', label='Vmem (mV)')
        plt.ylabel('Worm voltage in mv')
    else:	# plot an ion
        for ion in what_to_plot:
            data = cc_cells[ion_idx(ion)]
            plt.plot(data, linewidth=2.0, marker='.', label=ion)
        plt.ylabel (f"[{','.join(what_to_plot)}] in moles/m^3")

    plt.legend (loc="upper right")
    plt.xlabel('Cell # (tail is cell 0)', fontsize = 20)

    plot_it(title, filename)

# Graph one [ion] across the worm. But do it at multiple timepoints to see the
# curve shape evolve.
# 'Ion' is which ion to plot (text); t_shots and cc_shots are the usual summary
# of the simulation. N_views, t_start and t_end work together to pick the
# timepoints; if the sim went to t=100, then
#	(4,0,1) would plot at t=0,25,75,100
#	(3,.3,.5) would plot at t=30,40,50
def plot_worm_times (ion, t_shots, cc_shots, n_views,
			t_start=0, t_end=1, title=None, filename=None):
    plt.figure()

    # A bit of tricky coding here. First, 'times' will be the *requested*
    # timepoints (e.g., [30,40,50]). Then 't_indices' will be the indices (into
    # t_shots and cc_shots) of the snapshots closest to those times.
    # We use interp(), with the PWL function given by x values in t_shots and
    # y values just counting 0,1,..., to find the index closest to desired.
    times = np.linspace (t_start, t_end, n_views) * t_shots[-1]
    t_indices = np.interp (times, t_shots, range(len(t_shots)))
    t_indices = np.rint(t_indices).astype(int)	# round to integers

    ion_index = ion_idx(ion)
    for i in t_indices:
        data = cc_shots[i][ion_index]
        plt.plot(data, linewidth=1.5, marker='.', label=f't={t_shots[i]}')

    plt.xlabel('Cell # along the worm')
    plt.ylabel('['+ion+'] (mol/m3)')
    plt.legend (loc="upper right")

    plot_it(title, filename)

######################################################
# Pretty plotting.
# A pretty plot is a pretty representation of the network of cells at one point
# in time.
# - Each cell is drawn as a circle on the screen, with the cell's Vmem
#   written inside. Furthermore, the cells are colored according to their Vmem
#   (like a heat map).
# - GJs are drawn as lines between the appropriate cells.
######################################################

# Set the shape of a future output plot of the cell network.
# Inputs: shape
#	'shape' is a two-element array that sets the shape of an eventual
#	pretty plot. Each cell will be drawn as a circle; the circles will be
#	in a rectangular grid, with the GJs draw as lines connecting the cells.
#	Shape[0] tells how many rows in the plot; [1] tells how many columns.
#	With two rows of three columns, the cells will be drawn as
#	3 4 5
#	0 1 2
# Outputs: none, but the side effect of saving shape for a later plot call.
g_shape = None
def set_network_shape (shape):
    global g_shape
    g_shape = shape

# Draw a network with each cell as a circle, and each GJ as a line connecting
# its two cells.
# Label each cell with its index
# Color each cell based on whatever data we're given.
# Inputs: data[num_cells]
#	A 1D array of numbers, one per cell, saying what to plot for each cell.
def pretty_plot (data):
    global g_shape
    num_layers = g_shape[0]
    cells_per_layer = g_shape[1]

    # Assign a plot-worthy x,y coordinate pair to each cell.
    # Specifically, build xypts[N_CELLS,2]: each row #r is the (x,y) coordinates
    # of where to plot cell #r. Each layer of cells is a row in the plot, and
    # row #0 is at the bottom (with cell #0 at the left).
    # So if there are 2 cells per layer and 3 layers, then xypts is
    #	[[0. 0.] [1. 0.] [0. 1.] [1. 1.] [0. 2.] [1. 2.] ]
    # I.e., cell #0 is drawn at the lower left; then go left to right across
    # the bottom row, then left to right one row up, etc.
    # And num_layers(cells_per_layer) is the number of rows(column)
    y = np.linspace(0, num_layers-1, num_layers)
    x = np.linspace(0, cells_per_layer-1, cells_per_layer)
    X, Y = np.meshgrid(x, y)
    xypts = np.asarray([X.ravel(), Y.ravel()]).T

    # Line segments defining gap junctions (for plotting)
    # After this fancy indexing, GJ_segs is a 3D array [N_GJ,2,2]
    # Any piece GJ_segs[g] is a 2x2 array that describes how to draw GJ #g
    # as a line segment. I.e., the 1st row of the 2x2 array is the (x,y)
    # location of the cell for the GJ's input, and the 2nd is the (x,y) location
    # of the cell for the GJ's output.
    GJ_from_to = np.stack((sim.GJ_connects['from'],sim.GJ_connects['to']),1)
    GJ_segs = xypts[GJ_from_to]

    plt.figure()
    ax = plt.subplot(111)

    # Each cell is a circle:
    # - 'c' gives the colors; plt.colorbar() picks a color for each element of
    #   Vm (i.e., maps from our Vmem range to nice colors)
    # - 's' is the area of the circles.
    plt.scatter(xypts[:,0], xypts[:,1], c=data, s=500)
    # Draw a bar showing the mapping from Vmem to color
    plt.colorbar(orientation='horizontal')

    # Label each cell with its index number.
    for i, (xi, yi) in enumerate(xypts):
        label = "Cell #{:d}\n({:.2f})".format (i, (data[i]))
        print (label)
        plt.text(xi, yi, label, fontsize=14, fontweight='bold')

    # Draw the gap junctions
    GJs = LineCollection(GJ_segs)
    ax.add_collection(GJs)

    plt.axis('equal')	# ensure that circles are circular and not ellipses.
    plt.axis('off')
    plt.show()

def ion_idx (name):
    ion = sim.ion_i.get (name, -1)
    if (ion == -1):
        raise ValueError('Unknown ion "'+name+'" given to plot routines')
    return (ion)

def plot_it(title, filename):
    if (title != None):
        plt.title (title)
    if (filename is None):
        plt.show()
    else:
        plt.savefig(filename)
