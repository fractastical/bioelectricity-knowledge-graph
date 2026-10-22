ecg_filename = "scope_out.csv"
PQRST = (1,1.23,1.45,1.67,2.11)	# Put your own 5 numbers here

import csv
import numpy as np, math
import matplotlib.pyplot as plt

# Read an ASCII file where each line is lead1, lead2 (and remember that lead1
# is LA-RA and lead2 is LL-RA).
# Returns (lead1, lead2); if the file has 1000 lines, then lead1 and lead2
# are each a Numpy array of 1000 integers.
def read_file (ecg_filename):
    lines = []
    with open(ecg_filename, mode='r') as file:
        csvFile = csv.reader(file)
        lines = [l for l in csvFile]
    # Check for formatting
    for line in lines:
        assert (len(line)==2) or (len(line)==3 and line[2]==""), \
		f"Wrong number of fields on line *{lines[0]}*"

    lead1 = [int(line[0]) for line in lines]
    lead2 = [int(line[1]) for line in lines]

    lead1 = np.array (lead1)		# LA - RA
    lead2 = np.array (lead2)		# L - RA
    return (lead1, lead2)

# Given 'leads': a list of arrays, where each array is one lead.
# Find the lead with the biggest deltaY (i.e., the largest sample minus the
# smallest sample for that array) and return that deltaY. This will let us plot
# all the leads with the same Y-axis bounds.
def deltaY (leads):
    max_range = 0
    for lead in leads:
        range = lead.max() - lead.min()
        max_range = max (max_range, range)

    # Round up to the nearest 100 to keep the Y axis having nice round numbers.
    max_range = (math.ceil(max_range/100)*100)
    return (max_range)

# Plot one lead as its own subplot.
#  - Ax is the subplot to draw on.
#  - Lead_data is the data array
#  - Title is the plot title for this subplot.
#  - PQRST is a list of 5 "x" locations for P, Q, R, S and T
#  - dY is the deltaY to use for all six plots; this plot's Y range will be
#    from (its min Y) to (its min Y + dY).
def draw_lead (ax, lead_data, title, PQRST, dY):
    # The time axis is 2ms per sample point
    t = .002 * np.arange (lead_data.shape[0])

    # Set the Y min/max limits so that each of the six plots has the same scale.
    Ymin = lead_data.min()
    ax.set_ylim (Ymin, Ymin+dY)
    print (f"{title} range=[{Ymin}-{Ymin+dY}]")

    # Draw the main plot & set title
    ax.plot (t, lead_data)
    ax.set_title (title)

    # Draw dots at the P,Q,R,S and T locations.
    for t in PQRST:
        idx = (int) (t / .002)		# Again, .002 sec/sample
        y = lead_data[idx]		# The plot value at this T
        ax.plot ([t], [y], marker='o', markerfacecolor='red', markersize=5)

def main():
    (lead1, lead2) = read_file (ecg_filename)

    # Compute lead3. That's easy.
    lead3 = lead2 - lead1		# LL - LA

    # Define RA=0 by convention. Then lead1=LA and lead2=LL.
    # Wilson's central terminal = W = (LA +RA + LL)/3 = (lead1+0+lead2)/3
    LA=lead1; RA=np.zeros_like(lead1); LL=lead2
    W = (LA+RA+LL)/3

    aVR = 1.5 * (RA-W)
    aVL = 1.5 * (LA-W)
    aVF = 1.5 * (LL-W)

    # Put all six plots on the same Y scale. dY is the largest maxY-minY over
    # all six leads.
    dY = deltaY ((lead1,lead2,lead3,aVR,aVL,aVF))
    print (f"{dY=}")

    # sharex means that all subplots always share the same x axis; so when you
    # interactively zoom into one then the others all follow.
    fig,axes = plt.subplots (2,3, sharex=True)
    fig.suptitle ("All six leads")

    draw_lead (axes[0,0], lead1, "lead 1, LA-RA (0 CW)",    PQRST, dY)
    draw_lead (axes[0,1], lead2, "lead 2, LL-RA (60 deg)",  PQRST, dY)
    draw_lead (axes[0,2], aVF,   "aVF, LL-WCT (90 deg)",    PQRST, dY)
    draw_lead (axes[1,0], lead3, "lead 3, LL-LA (120 deg)", PQRST, dY)
    draw_lead (axes[1,1], aVR,   "aVR, RA-WCT (-150 deg)",  PQRST, dY)
    draw_lead (axes[1,2], aVL,   "aVL, LA-WCT (-30 deg)",   PQRST, dY)

    axes[0,0].tick_params (labelbottom=True)
    axes[0,1].tick_params (labelbottom=True)
    axes[0,2].tick_params (labelbottom=True)
    plt.show()

main()
