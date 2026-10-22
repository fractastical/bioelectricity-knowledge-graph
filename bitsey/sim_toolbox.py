#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter


# Toolbox of functions used in the Simulator class to calculate key bioelectric properties.

# This is the new/improved version of GHK(). The only difference is that it
# works on entire 2D arrays at once, rather than one ion at a time.
def GHK (cExt,cIn,D,ion_z,Vmem,p):
    """
    Goldman-Hodges-Katz between two connected volumes.
    Return the flux moving *into* the cell

    This function simplifies to regular diffusion if Vba == 0.0.

    This function takes numpy matrix values as input.

    This is the Goldman Flux/Current Equation (not to be confused with the
    Goldman Equation). Note: the Nernst-Planck equation has been trialed in
    place of this, and it does not reproduce proper reversal potentials.

    Parameters
    ----------
    cExt[n_ions]	concentration in the ECF [moles/m3]
    cIn[n_ions,n_cells]	concentration in the ICF [moles/m3]
    D[n_ions,n_cells]	Diffusion constant of the ion  [m2/s]
    ion_z[n_ions]	valence of the ion
    Vmem[n_cells]	voltage difference vIn - vOut
    p			an instance of the Parameters class

    Returns
    --------
    flux[n_ions,n_cells] Chemical flux magnitude flowing *into* the cell [mol/s]

    """
    # Multiply z*Vmem. But Vmem changes cell by cell (and is thus []), and z is
    # [i]. So the result is [c,i] and has potentially a different number in each
    # entry. Now, multiplying [i]*[c] is illegal. But [i,1]*[c,1], by the rules
    # of broadcasting, correctly becomes [i,c] as desired.
    # Second issue: GHK becomes 0/0 when z*Vmem=0. We avoid this by adding a
    # small FLOAT_DELTA to every element of ZVrel.
    FLOAT_DELTA = 1.0e-25
    i = ion_z.size; c=Vmem.size
    ZVrel = (np.reshape(ion_z,(i,1))*np.reshape(Vmem*p.k26mV_inv,(1,c))) \
          + FLOAT_DELTA

    exp_ZVrel = np.exp(-ZVrel)		# exp (-Z*Vrel)
    deno = -np.expm1(-ZVrel)		# 1 - exp(-Z*Vrel), the GHK denominator
    # tm = thickness of the cell membrane [meters]
    P = D / p.tm 			# permeability of a channel = D/L; [i,c]

    # -P*Z*Vrel * (cIn - cExt*exp(-Z*Vrel)) / (1 - exp(-Z*Vrel))
    cExt_exp = (cExt.T * exp_ZVrel.T).T
    flux = -P*ZVrel*((cIn - cExt_exp)/deno)

    return flux

# This is the old (and unused) version of GHK. It takes 1D inputs.
def GHK_old(cExt,cIn,D,zc,Vmem,p):
    """
    Goldman-Hodges-Katz between two connected volumes.
    Return the flux moving *into* the cell

    This function defaults to regular diffusion if Vba == 0.0.

    This function takes numpy matrix values as input. All inputs must be
    matrices of the same shape.

    This is the Goldman Flux/Current Equation (not to be confused with the
    Goldman Equation). Note: the Nernst-Planck equation has been trialed in
    place of this, and it does not reproduce proper reversal potentials.

    Parameters
    ----------
    cExt	concentration in the ECF [moles/m3]
    cIn		concentration in the ICF [moles/m3]
    D		Diffusion constant of the ion  [m2/s]
    zc		valence of the ion
    Vmem	voltage difference vIn - vOut
    p           an instance of the Parameters class

    Returns
    --------
    flux        Chemical flux magnitude flowing *into* the cell [mol/s]

    """

    # GHK becomes 0/0 when Vmem=0. Avoid this by trying to never encounter
    # Vmem=0. I'm not sure why Alexis added FLOAT_NONCE to zc also.
    FLOAT_NONCE = 1.0e-25
    Vmem += FLOAT_NONCE
    zc += FLOAT_NONCE

    ZVrel = zc*Vmem* p.k26mV_inv	# Z * (Vmem/26mV), or Z*Vrel
    exp_ZVrel = np.exp(-ZVrel)		# exp (-Z*Vrel)
    deno = -np.expm1(-ZVrel)		# 1 - exp(-Z*Vrel), the GHK denominator
    # tm= thickness of the cell membrane [m]
    P = D / p.tm			# permeability of a channel = D/L

    # -P*Z*Vrel * (cIn - cExt*exp(-Z*Vrel)) / (1 - exp(-Z*Vrel))
    flux = -P*ZVrel*((cIn -cExt*exp_ZVrel)/deno)

    return flux

def pumpNaKATP(cNai,cNao,cKi,cKo,Vm,T,p,block, met = None):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object

    met             A "metabolism" vector containing concentrations of ATP, ADP and Pi


    Returns
    -------
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    deltaGATP_o = p.deltaGATP  # standard free energy of ATP hydrolysis reaction in J/(mol K)

    cATP = p.cATP
    cADP = p.cADP
    cPi  = p.cPi

    # calculate the reaction coefficient Q:
    Qnumo = (cADP*1e-3)*(cPi*1e-3)*((cNao*1e-3)**3)*((cKi*1e-3)**2)
    Qdenomo = (cATP*1e-3)*((cNai*1e-3)**3)*((cKo*1e-3)** 2)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-15

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-(deltaGATP_o / (p.R * T) - ((p.F * Vm) / (p.R * T))))

    # calculate the enzyme coefficient:
    numo_E = ((cNai/p.KmNK_Na)**3) * ((cKo/p.KmNK_K)**2) * (cATP/p.KmNK_ATP)
    denomo_E = (1 + (cNai/p.KmNK_Na)**3)*(1+(cKo/p.KmNK_K)**2)*(1+(cATP/p.KmNK_ATP))

    fwd_co = numo_E/denomo_E

    f_Na = -3*block*p.alpha_NaK*fwd_co*(1 - (Q/Keq))  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    f_K = -(2/3)*f_Na          # flux as [mol/m2s]

    return f_Na, f_K, -f_Na  # FIXME get rid of this return of extra -f_Na!!
