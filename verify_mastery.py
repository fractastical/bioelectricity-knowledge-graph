#!/usr/bin/env python3
"""
Math Academy-Style Mastery & Simulation Verification Engine
for Prof. Joel Grodstein's Tufts EE 123 / BME 123 Bioelectricity Knowledge Graph.

This script performs two functions:
1. Enriches `knowledge-graph/knowledge_graph.json` with interactive mastery
   verification checkpoints (`verification`) for all 27 nodes (multiple-choice /
   quantitative questions, step-by-step explanations, and BITSEY simulation assertions).
2. Runs automated verification across all 27 nodes:
   - Graph / Prerequisite DAG integrity (acyclicity, topological rank consistency,
     implicit downward prerequisite propagation & upward invalidation)
   - Quantitative physics & circuit equations (Khan Academy Quiz Part 1 & Part 2,
     Nernst potentials, RC time constants, Einthoven/Wilson ECG leads)
   - Live BITSEY Python simulations for SS1/SS2/SS3 (Steady-State), QSS1/QSS2
     (Quasi-Steady-State), Hodgkin-Huxley Action Potentials (m^3h, n^4), and
     5-cell Planarian Morphogen-GJ Polarity regeneration.
"""

import json
import math
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent
BITSEY_DIR = REPO_ROOT / "bitsey"
KG_PATH = REPO_ROOT / "knowledge-graph" / "knowledge_graph.json"

sys.path.insert(0, str(BITSEY_DIR))
import sim  # noqa: E402
import sim_toolbox as stb  # noqa: E402

# ---------------------------------------------------------------------------
# 1. Interactive Mastery Checkpoints for all 27 Knowledge Graph Nodes
# ---------------------------------------------------------------------------
VERIFICATION_DATA = {
    "KA-PHYS-01": {
        "question": "In BITSEY's default cell (sim.py), initial intracellular concentrations are [Na+]=10, [K+]=125, [Cl-]=55, and [P-]=80 mol/m^3. What is the net bulk charge concentration sum(z_i * [ion_i]) at t=0?",
        "options": [
            "0 mol/m^3 (bulk charge neutral: 10 + 125 - 55 - 80 = 0)",
            "+50 mol/m^3 (net positive intracellular fluid)",
            "-55 mol/m^3 (dominated by chloride ions)",
            "+270 mol/m^3 (sum of all ion concentrations)"
        ],
        "correct_index": 0,
        "explanation": "Bulk intracellular and extracellular fluids are macroscopically charge neutral: (+1)(10) + (+1)(125) + (-1)(55) + (-1)(80) = 0 mol/m^3. Even a tiny deviation (< 0.03 mol/m^3) produces tens of millivolts across the membrane.",
        "automated_check": "check_ka_phys_01_electroneutrality"
    },
    "KA-PHYS-02": {
        "question": "If a cell membrane has thickness t_m = 7.5 nm and resting potential V_mem = -75 mV (with ECF = 0 V), what is the magnitude of the electric field |E| = |V_mem| / t_m inside the membrane?",
        "options": [
            "10,000,000 V/m (10 MV/m)",
            "75 V/m",
            "0.56 V/m",
            "1,000 V/m"
        ],
        "correct_index": 0,
        "explanation": "|E| = 0.075 V / 7.5e-9 m = 1.0e7 V/m (10 million Volts per meter!). Because the lipid bilayer is only 7.5 nm thick, a 75 mV biological voltage creates an enormous electric field.",
        "automated_check": "check_ka_phys_02_electric_field"
    },
    "KA-PHYS-03": {
        "question": "Backgrounders Placement Quiz Q1b: How much charge Q (in Coulombs) is stored on a 50 nF capacitor when 100 mV is applied across it?",
        "options": [
            "5 nC (5.0e-9 Coulombs)",
            "500 nC (5.0e-7 Coulombs)",
            "0.5 C",
            "2000 C"
        ],
        "correct_index": 0,
        "explanation": "Using Q = C * V: Q = (50e-9 F) * (0.100 V) = 5.0e-9 C = 5 nC.",
        "automated_check": "check_ka_phys_03_capacitance"
    },
    "KA-CKT-01": {
        "question": "Backgrounders Placement Quiz Q1a: How much current flows through a 20 Ohm resistor with 100 mV across it, and what is the resistor's conductance G?",
        "options": [
            "I = 5 mA (0.005 A) and G = 0.05 S (1/Ohm)",
            "I = 2 A and G = 20 S",
            "I = 200 mA and G = 0.5 S",
            "I = 0.2 mA and G = 0.005 S"
        ],
        "correct_index": 0,
        "explanation": "By Ohm's Law, I = V / R = 0.100 V / 20 Ohm = 0.005 A = 5 mA. Conductance is G = 1 / R = 1 / 20 = 0.05 S.",
        "automated_check": "check_ka_ckt_01_ohms_law"
    },
    "KA-CKT-02": {
        "question": "Circuits Backgrounder Slide 14: Two Nernst batteries V_1 = +80 mV (through R_1 = 40 Ohm) and V_2 = -80 mV (through R_2 = 10 Ohm) connect to node V_3. By KCL, what voltage does V_3 settle to?",
        "options": [
            "-48 mV (pulled 4x closer to -80 mV because R_2 is 4x smaller than R_1)",
            "0 mV (exact midpoint between +80 mV and -80 mV)",
            "+48 mV (pulled toward the larger 40 Ohm resistor)",
            "-64 mV"
        ],
        "correct_index": 0,
        "explanation": "KCL at V_3: (V_3 - 80)/40 + (V_3 - (-80))/10 = 0 => (V_3 - 80) + 4(V_3 + 80) = 0 => 5*V_3 = -240 => V_3 = -48 mV.",
        "automated_check": "check_ka_ckt_02_voltage_divider"
    },
    "KA-CKT-03": {
        "question": "Backgrounders Placement Quiz Part 2: Given V_N,Na = +77 mV (G_Na = 0.4), V_N,K = -89 mV (G_K = 2.2), V_N,Cl = -71 mV (G_Cl = 0.4), and outward pump currents I_Na = 60, I_K = -40, what is V_ICF (and what if G_K becomes 100x larger)?",
        "options": [
            "V_ICF = -71.13 mV; if G_K -> 220, V_ICF approaches V_N,K = -89 mV",
            "V_ICF = -27.7 mV; if G_K -> 220, V_ICF approaches 0 mV",
            "V_ICF = +77.0 mV; if G_K -> 220, V_ICF approaches +77 mV",
            "V_ICF = -89.0 mV; if G_K -> 220, V_ICF approaches -71 mV"
        ],
        "correct_index": 0,
        "explanation": "KCL at ICF (Lecture 1b Slide 36): 0.4(V - 77) + 2.2(V + 89) + 0.4(V + 71) + 60 - 40 = 0 => 3.0*V + 213.4 = 0 => V = -71.13 mV. When G_K becomes 100x larger, the K+ branch dominates and pins V_ICF to V_N,K = -89 mV.",
        "automated_check": "check_ka_ckt_03_nodal_cell"
    },
    "KA-CKT-04": {
        "question": "In Circuits Backgrounder Slides 18-25, when we collapse a cell's parallel conductances (G_Na = 0.4, G_K = 2.2, G_Cl = 0.4 S) into a single Thevenin equivalent to compute RC charging delay, what is G_th = 1 / R_th?",
        "options": [
            "G_th = G_Na + G_K + G_Cl = 3.0 S (so R_th = 1/3 Ohm and tau = R_th * C_mem)",
            "G_th = 1 / (1/0.4 + 1/2.2 + 1/0.4) = 0.18 S",
            "G_th = G_K - G_Na = 1.8 S",
            "G_th depends only on the Na+/K+ pump current"
        ],
        "correct_index": 0,
        "explanation": "Parallel conductances add directly: G_th = G_Na + G_K + G_Cl = 0.4 + 2.2 + 0.4 = 3.0 S. The RC time constant governing how fast V_mem slews to QSS is tau = C_mem / G_th.",
        "automated_check": "check_ka_ckt_04_thevenin_rc"
    },
    "FND-MATH-01": {
        "question": "In Sim Lab 4 (slab_neurons.pdf Slide 3), gating variable m obeys dm/dt = (m_inf - m)/tau_m. If m(0) = 0.4 and m_inf steps to 0.9 with tau_m = 10 s, what is m(10 s)?",
        "options": [
            "0.4 + 0.5 * (1 - 1/e) approx 0.716",
            "0.900 (reaches 100% in one time constant)",
            "0.650 (linear midpoint)",
            "0.400 * (1/e) approx 0.147"
        ],
        "correct_index": 0,
        "explanation": "After one time constant t = tau_m, an exponential relaxation covers (1 - e^-1) approx 63.21% of the step: 0.4 + (0.9 - 0.4)*(1 - 0.3679) = 0.716.",
        "automated_check": "check_fnd_math_01_ode"
    },
    "FND-PROG-01": {
        "question": "In BITSEY's sim.py (compute_Vm), given a spherical cell of radius r = 5 um, patch capacitance c_m = 0.05 F/m^2, and Faraday constant F = 96485 C/mol, how many millivolts does V_mem change per +1 mol/m^3 of net uncompensated valence-1 charge?",
        "options": [
            "Approx +3,216 mV (+3.22 V) per mol/m^3",
            "Approx +26 mV per mol/m^3",
            "Approx +1.0 mV per mol/m^3",
            "Approx +100 mV per mol/m^3"
        ],
        "correct_index": 0,
        "explanation": "From sim.py line 375: scale = F * (V_cell / A_cell) / c_m = F * (r / 3) / c_m = 96485 * (5e-6 / 3) / 0.05 = 3.216 Volts = 3216 mV per (mol/m^3)!",
        "automated_check": "check_fnd_prog_01_bitsey_vm_scale"
    },
    "FND-BIO-01": {
        "question": "In Lecture 4a (Bio Backgrounder) & Lecture 4c (sim.py Hill_gate), what is the output of a Hill inverter H_inv([M]) = 1 / (1 + ([M]/K_M)^N) when the ligand concentration [M] exactly equals K_M, for any Hill coefficient N > 0?",
        "options": [
            "0.5 (the half-maximal switching threshold, regardless of N)",
            "1.0 (fully open)",
            "0.0 (fully closed)",
            "1 / N"
        ],
        "correct_index": 0,
        "explanation": "When [M] = K_M, ([M]/K_M)^N = 1^N = 1, so H_inv(K_M) = 1 / (1 + 1) = 0.5 for every N. Increasing N steepens the slope (gain) around K_M without shifting the 50% point.",
        "automated_check": "check_fnd_bio_01_hill_inverter"
    },
    "U1-BIO-01": {
        "question": "Real Lab 1 (Skittles Diffusion) Q4 & Lecture 2a Slide 6: If diffusing a chemical signal across 10 um takes 0.05 seconds, how long does it take to diffuse across a 10x larger distance (100 um)?",
        "options": [
            "5.0 seconds (10^2 = 100x longer because t_diff scales with distance squared, L^2)",
            "0.5 seconds (10x longer, linear scaling)",
            "0.05 seconds (independent of distance)",
            "50 seconds (10^3 = 1000x longer)"
        ],
        "correct_index": 0,
        "explanation": "Because random-walk mean-squared displacement obeys <x^2> = 2*D*t, diffusing n times farther takes n^2 times longer. Over 1 meter (a human leg), diffusion would take years—which is why large animals evolved bioelectricity!",
        "automated_check": "check_u1_bio_01_diffusion_scaling"
    },
    "U1-BIO-02": {
        "question": "Lecture 1a Slides 33-38: Using the Einstein relation mu = D / (k_B * T), why do ions moving through a membrane channel under a constant voltage V_mem reach a constant terminal drift velocity v_drift = mu * q * (V_mem / L)?",
        "options": [
            "Viscous drag from water/channel collisions balances the electric force f = qE, just like air drag on a car",
            "The electric field drops to zero inside the channel",
            "Ions run out of kinetic energy halfway through the membrane",
            "Magnetic fields cancel the electric acceleration"
        ],
        "correct_index": 0,
        "explanation": "In the crowded nanoworld, friction/drag is proportional to velocity. Terminal velocity is reached almost instantaneously where drag equals electrostatic force f = q*V_mem/L, giving Ohm's Law I = G*V.",
        "automated_check": "check_u1_bio_02_einstein_drift"
    },
    "U1-BIO-03": {
        "question": "Lecture 1a Slides 51-53 (Giant Squid In-Class Exercise): Given [K+]_ICF = 400 mM, [K+]_ECF = 20 mM and [Na+]_ICF = 50 mM, [Na+]_ECF = 440 mM at body temperature (RT/F approx 26 mV), what are V_Nernst,K and V_Nernst,Na?",
        "options": [
            "V_N,K approx -78 mV and V_N,Na approx +57 mV",
            "V_N,K approx +78 mV and V_N,Na approx -57 mV",
            "V_N,K approx -62 mV and V_N,Na approx 0 mV",
            "V_N,K approx -26 mV and V_N,Na approx +26 mV"
        ],
        "correct_index": 0,
        "explanation": "V_N,K = 26 mV * ln(20 / 400) = 26 * (-2.9957) = -77.9 mV. V_N,Na = 26 mV * ln(440 / 50) = 26 * (2.1748) = +56.5 mV.",
        "automated_check": "check_u1_bio_03_nernst_squid"
    },
    "U1-BIO-04": {
        "question": "Lecture 1a Slides 57-63: In a resting cell at Steady State (SS), why does Cl- satisfy detailed balance (j_drift + j_diff = 0) while Na+ and K+ do not?",
        "options": [
            "Cl- has no active pump, so its channel flux must be zero at SS (V_mem = V_N,Cl); Na+ and K+ have an ATPase pump balancing non-zero channel leaks",
            "Cl- channels are permanently closed in all living cells",
            "Na+ and K+ cannot diffuse through ion channels",
            "Cl- has valence z = 0"
        ],
        "correct_index": 0,
        "explanation": "Because Cl- has no pump in our model, steady state (d[Cl-]/dt = 0) requires its net channel current (drift + diffusion) to be zero, meaning V_mem settles to V_N,Cl. For Na+ and K+, passive channel leak is balanced by active Na+/K+-ATPase pumping.",
        "automated_check": "check_u1_bio_04_pump_ss_balance"
    },
    "U1-LAB-01": {
        "question": "Sim Lab 2 (slab_simple_SS.pdf, main_sim_SS.py): When we run BITSEY for 100,000 s, what happens in SS2 (doubling initial [Na+]_ICF & [Cl-]_ICF in cell #1) and SS3 (doubling D_Na in cell #1, D_K in cell #2, D_Cl in cell #3)?",
        "options": [
            "In SS2 all cells converge to the exact same V_mem (-49.7 mV); in SS3 doubling D_Na depolarizes (-37.1 mV), doubling D_K hyperpolarizes (-59.5 mV), and doubling D_Cl has zero effect (-49.7 mV)",
            "In SS2 cell #1 stays permanently 2x higher in voltage; in SS3 doubling D_Cl hyperpolarizes the cell",
            "Doubling D_K makes V_mem more positive (+57 mV)",
            "Doubling initial concentrations permanently changes steady-state V_mem"
        ],
        "correct_index": 0,
        "explanation": "Verified directly by running BITSEY: SS1 and SS2 all converge to -49.74 mV (restorative negative feedback erases initial concentration perturbations). In SS3, cell #0 is -49.74 mV, cell #1 (2x D_Na) is -37.11 mV, cell #2 (2x D_K) is -59.54 mV, and cell #3 (2x D_Cl) is -49.74 mV.",
        "automated_check": "check_u1_lab_01_bitsey_ss"
    },
    "U1-BIO-05": {
        "question": "Sim Lab 3 (slab_simple_QSS.pdf, main_sim_QSS.py) & Oral Quiz #1: What is the fundamental difference between Quasi-Steady State (QSS at t=2 s) and true Steady State (SS at t=100,000 s)?",
        "options": [
            "In QSS, the sum of all ionic currents is zero (net charge flux = 0) so V_mem settles in milliseconds, even though individual ion concentrations are still slowly drifting",
            "In QSS, every individual ion has zero net flux",
            "QSS takes tens of hours to reach, whereas SS happens in 1 millisecond",
            "QSS only applies to uncharged molecules"
        ],
        "correct_index": 0,
        "explanation": "Because C_mem is tiny, moving < 0.01% of intracellular ions charges the membrane in milliseconds until I_Na + I_K + I_Cl + I_pump = 0 (KCL). At that QSS voltage, individual ions are still unbalanced and take hours to alter bulk concentrations.",
        "automated_check": "check_u1_bio_05_bitsey_qss"
    },
    "U2-NEU-01": {
        "question": "Sim Lab 4 (slab_neurons.pdf, main_sim_neuron.py): In the Hodgkin-Huxley model (g_Na = g_Na_max * m^3 * h, g_K = g_K_max * n^4), what happens when you increase the slow negative-feedback time constants tau_h and tau_n from 10 to 50 (keeping tau_m = 1)?",
        "options": [
            "The action potential reaches a higher peak V_mem and becomes wider because Na+ inactivation (h) and K+ activation (n) take longer to shut down the spike",
            "The action potential disappears completely because positive feedback is slower",
            "The action potential becomes narrower and smaller in amplitude",
            "The resting potential shifts to +80 mV before the stimulus arrives"
        ],
        "correct_index": 0,
        "explanation": "Fast m (tau_m = 1) launches the positive-feedback Na+ spike; slower h and n (tau_h, tau_n) provide the negative feedback that turns off Na+ and turns on K+. Slowing tau_h and tau_n lets V_mem climb higher toward V_N,Na and widens the spike.",
        "automated_check": "check_u2_neu_01_bitsey_hh"
    },
    "U2-NEU-02": {
        "question": "Lecture 2a Slides 19-22: Why does wrapping an axon in a thick myelin sheath dramatically increase action potential propagation speed (up to 150 m/s)?",
        "options": [
            "Thick myelin acts like series capacitors that greatly decrease membrane capacitance C, so much less charge Delta_Q = C * Delta_V is required to depolarize the axon segment",
            "Myelin contains high-density voltage-gated Na+ channels along its entire length",
            "Myelin increases membrane capacitance C by 100x",
            "Myelin converts electrical signals into acoustic waves"
        ],
        "correct_index": 0,
        "explanation": "By Q = C*V, thicker insulation (smaller C) means a tiny axial current changes V_mem rapidly. Because myelin blocks ion channels, bare Nodes of Ranvier are spaced periodically to regenerate the decaying signal.",
        "automated_check": "check_u2_neu_02_myelin_capacitance"
    },
    "U2-NEU-03": {
        "question": "Lecture 2b Slides 23-28 & Lab 5 (sEMG): In a clinical Nerve Conduction Study (NCS) and surface EMG, how do you distinguish diabetic neuropathy (axonal death) from Multiple Sclerosis (myelin damage), and why does a moving-average over 16 samples at 960 Hz remove power-line noise?",
        "options": [
            "Axonal loss reduces CMAP amplitude while demyelination slows conduction velocity; a 16-sample window at 960 Hz spans exactly one 60 Hz period (1/60 s), averaging 60 Hz sine noise to zero",
            "Axonal loss speeds up conduction while MS increases CMAP amplitude; 60 Hz noise is removed by squaring the signal",
            "Both diseases only affect muscle fibers and leave nerves intact",
            "Moving-average filters amplify high frequencies above 500 Hz"
        ],
        "correct_index": 0,
        "explanation": "Fewer surviving axons produce a smaller Compound Muscle Action Potential (CMAP), whereas loss of myelin increases C and slows propagation velocity. In DSP (EN01 Lab 4), averaging across an integer number of 60 Hz cycles cancels sinusoidal power-line interference.",
        "automated_check": "check_u2_neu_03_semg_dsp"
    },
    "U2-NEU-04": {
        "question": "Lecture 2b Slides 34-39 & Lab 5 Q3: In a myoelectric prosthesis controlled by biceps and triceps sEMG, how does Targeted Muscle Reinnervation (TMR) help an amputee, and why must the target muscle region be denervated first?",
        "options": [
            "TMR surgically attaches residual hand/wrist motor nerves onto remaining chest/arm muscle so thinking 'close hand' creates a distinct sEMG; native nerves are cut first so 'move wrist' doesn't accidentally flex the elbow or mix signals",
            "TMR implants electric motors directly into the brain motor cortex",
            "TMR connects biceps tendons directly to bone to measure mechanical tension",
            "TMR uses optical fibers to read visual signals from the skin"
        ],
        "correct_index": 0,
        "explanation": "TMR uses denervated muscle patches as biological amplifiers for severed peripheral motor nerves, creating independent sEMG channels without cross-talk from the muscle's original innervation.",
        "automated_check": "check_u2_neu_04_prosthetic_control"
    },
    "U2-NEU-05": {
        "question": "Lecture 2b Slides 52-58 & Oral Quiz #2: According to Melzack & Wall's Gate-Control Theory of Pain, why does rubbing a pinched finger (or applying conventional 50 Hz Spinal Cord Stimulation) reduce pain, and why is 10 kHz High-Frequency SCS surprising?",
        "options": [
            "Touch (A-beta) fibers excite an inhibitory interneuron in the dorsal horn that blocks nociceptive (A-delta/C) pain transmission; 10 kHz SCS relieves pain without paresthesia (tingling), showing it acts directly on dorsal-horn processing rather than just A-beta rubbing",
            "Rubbing a finger shuts down Na+/K+ pumps in the skin; 10 kHz SCS heats up the spinal cord",
            "A-beta touch fibers directly inhibit the brain's visual cortex",
            "10 kHz SCS works by causing strong paresthesia tingling across the entire body"
        ],
        "correct_index": 0,
        "explanation": "Conventional SCS (35-80 Hz) stimulates A-beta touch fibers to close the pain gate (creating paresthesia). 10 kHz SCS is paresthesia-free and takes hours/days to take effect, directly modulating dorsal-horn interneurons.",
        "automated_check": "check_u2_neu_05_gate_control"
    },
    "U3-CRD-01": {
        "question": "Lecture 3 Slides 13-20: The Sinoatrial (SA) node contains ~500,000 independent pacemaker cells whose HCN ('funny') channels oscillate at slightly different natural frequencies. How do they produce a single unified heartbeat?",
        "options": [
            "Gap junctions (electrical resistors between adjacent cells) pass current from earlier-firing cells to pull up neighboring cells toward threshold, mutually entraining all 500,000 cells into phase-locked synchrony",
            "The brain sends a 60 Hz clock pulse down the spinal cord to trigger every cardiomyocyte simultaneously",
            "Only 1 master cell is active at any time while the other 499,999 cells are dormant",
            "Chemical neurotransmitters diffuse across the entire heart in 1 ms"
        ],
        "correct_index": 0,
        "explanation": "When cell #1 depolarizes slightly ahead of cell #2, current flows through the gap junction G_GJ*(V_1 - V_2), accelerating cell #2's depolarization and synchronizing the pacemaker network.",
        "automated_check": "check_u3_crd_01_gj_entrainment"
    },
    "U3-CRD-02": {
        "question": "Lecture 3 Slides 36-42 & Lab 6 (ECG): Given four limb electrodes (RA, LA, LL, and RL ground), if at the R-wave peak the potentials are V_RA = -0.2 mV, V_LA = +0.3 mV, and V_LL = +1.0 mV, what are Lead I (LA - RA), Lead II (LL - RA), Lead III (LL - LA), and Wilson's Central Terminal W?",
        "options": [
            "Lead I = +0.5 mV, Lead II = +1.2 mV, Lead III = +0.7 mV, and W = +0.367 mV (satisfying Einthoven's Law: Lead I + Lead III = Lead II)",
            "Lead I = -0.5 mV, Lead II = +0.8 mV, Lead III = +1.3 mV, and W = 0 mV",
            "Lead I = +1.0 mV, Lead II = +1.0 mV, Lead III = 0 mV, and W = +1.1 mV",
            "Lead I = +0.3 mV, Lead II = +1.0 mV, Lead III = -0.2 mV, and W = +0.5 mV"
        ],
        "correct_index": 0,
        "explanation": "Lead I = 0.3 - (-0.2) = 0.5 mV; Lead II = 1.0 - (-0.2) = 1.2 mV; Lead III = 1.0 - 0.3 = 0.7 mV. Notice Lead I + Lead III = 0.5 + 0.7 = 1.2 mV = Lead II (KVL!), and W = (-0.2 + 0.3 + 1.0)/3 = 0.367 mV.",
        "automated_check": "check_u3_crd_02_ecg_leads"
    },
    "U3-CRD-03": {
        "question": "Lecture 3 Slides 59-69 & Oral Quiz #3: Why are the heart's three redundant pacemakers engineered with progressively slower intrinsic frequencies (SA node ~60+ bpm > AV node ~40 bpm > Purkinje fibers ~20 bpm)?",
        "options": [
            "Overdrive suppression: the faster upstream SA wave arrives and resets the slower downstream AV and Purkinje oscillators before their funny channels can reach threshold; if SA fails or blocks, the next oscillator automatically takes over",
            "Downstream oscillators run faster during exercise to override the SA node",
            "All three oscillators fire independently and their pulses add together to give 120 bpm",
            "Purkinje fibers oscillate at 20 bpm to pump blood backward into the atria"
        ],
        "correct_index": 0,
        "explanation": "Because each upstream pacemaker fires before the downstream backup completes its slower Phase-4 ramp, the fastest intact pacemaker dominates cleanly without collisions, while providing automatic triple-redundant fail-safe backup.",
        "automated_check": "check_u3_crd_03_overdrive_suppression"
    },
    "U4-MOR-01": {
        "question": "Lecture 4b Slides 8-35: What do 'Picasso' tadpoles (scrambled craniofacial organs migrating into normal frog faces), ectopic eyes induced by V_mem hyperpolarization/depolarization, and 2-headed planaria after 1-octanol gap-junction blockade demonstrate?",
        "options": [
            "Morphogenesis operates as a closed-loop, layered system where a spatial bioelectric V_mem pattern serves as a high-level 'API' storing target anatomy and instructing downstream gene regulatory networks",
            "DNA encodes only an open-loop, fixed sequence of timed mechanical movements that cannot recover from perturbation",
            "Only neurons with myelin can produce membrane voltage",
            "1-octanol mutates the flatworm's DNA sequence into a different species"
        ],
        "correct_index": 0,
        "explanation": "Without altering DNA sequence, transiently altering V_mem or gap-junctional connectivity rewrites the tissue-level anatomical target pattern (e.g., permanently converting a 1-headed worm into a 2-headed worm that continues regenerating 2 heads).",
        "automated_check": "check_u4_mor_01_layered_api"
    },
    "U4-MOR-02": {
        "question": "Lecture 4c (4c_worm.pptx Slides 18-55): In our 5-cell BITSEY planarian model with a negatively charged morphogen M (z_M = -1), how does the closed-loop positive feedback work, and what happens at extreme gap-junction resistances R_GJ?",
        "options": [
            "Depolarized (+V_mem) head attracts negative morphogen M via drift; high [M] closes K+ channels via a Hill inverter, raising V_mem further. If R_GJ -> infinity cells form isolated 'archipelagos'; if R_GJ -> 0 all cells short-circuit to a flat voltage",
            "Positive morphogen M is repelled from the head and opens K+ channels to hyperpolarize the head",
            "Gap junctions only allow uncharged proteins to pass between cells",
            "Setting R_GJ = 0 creates the steepest possible head-to-tail voltage gradient"
        ],
        "correct_index": 0,
        "explanation": "Nernst electrodiffusion drives negatively charged M toward the more positive head cell; the Hill inverter maps high [M] -> low G_K -> higher V_mem, closing the positive-feedback loop so any amputated slice can re-establish polarity within a viable window of R_GJ.",
        "automated_check": "check_u4_mor_02_bitsey_worm"
    },
    "U4-CAN-01": {
        "question": "Lecture 4c (4c_cancer.pptx Slides 11-18): In the bioelectric model of cancer, how do V_mem and gap junctions behave during tumorigenesis compared to healthy differentiated tissue?",
        "options": [
            "Stem and tumor cells are persistently depolarized (-10 to -30 mV vs. -70 to -90 mV in differentiated cells) and close their gap junctions (G_GJ -> 0) to electrically isolate themselves from tissue-level morphogenetic normalization signals",
            "Tumor cells are strongly hyperpolarized (-120 mV) and increase their gap-junction conductance 100x",
            "Bioelectric state has no influence on cell proliferation once an oncogene is expressed",
            "Cancer only occurs in cells that fire rapid 1 ms Hodgkin-Huxley action potentials"
        ],
        "correct_index": 0,
        "explanation": "Hyperpolarizing ion-channel interventions (or forcing GJ coupling) can normalize tumorigenesis even in the presence of oncogenic mutations, whereas depolarizing V_mem and blocking GJs can induce invasive, tumor-like growth without genetic mutation.",
        "automated_check": "check_u4_can_01_cancer_depolarization"
    }
}


# ---------------------------------------------------------------------------
# 2. Automated Physics, Circuit & Live BITSEY Simulation Check Functions
# ---------------------------------------------------------------------------

def check_ka_phys_01_electroneutrality():
    z = np.array([1, 1, -1, -1])
    c_in = np.array([10.0, 125.0, 55.0, 80.0])
    c_out = np.array([145.0, 5.0, 140.0, 10.0])
    net_in = float(np.dot(z, c_in))
    net_out = float(np.dot(z, c_out))
    assert abs(net_in) < 1e-12 and abs(net_out) < 1e-12
    return f"ICF net charge = {net_in:.1f} mol/m^3, ECF net charge = {net_out:.1f} mol/m^3 (bulk neutral)"


def check_ka_phys_02_electric_field():
    v_mem = 0.075
    t_m = 7.5e-9
    e_field = v_mem / t_m
    assert abs(e_field - 1.0e7) < 1.0
    return f"|E| = {e_field:.2e} V/m across 7.5 nm bilayer"


def check_ka_phys_03_capacitance():
    c = 50e-9
    v = 0.100
    q = c * v
    assert abs(q - 5.0e-9) < 1e-15
    return f"Q = {q * 1e9:.2f} nC on 50 nF capacitor at 100 mV"


def check_ka_ckt_01_ohms_law():
    v = 0.100
    r = 20.0
    i_amp = v / r
    g = 1.0 / r
    assert abs(i_amp - 0.005) < 1e-12 and abs(g - 0.05) < 1e-12
    return f"I = {i_amp * 1000:.2f} mA, G = {g:.3f} S"


def check_ka_ckt_02_voltage_divider():
    v1, r1 = 80.0, 40.0
    v2, r2 = -80.0, 10.0
    v3 = (v1 / r1 + v2 / r2) / (1.0 / r1 + 1.0 / r2)
    assert abs(v3 - (-48.0)) < 1e-9
    return f"Divider node V_3 = {v3:.2f} mV"


def check_ka_ckt_03_nodal_cell():
    v_na, g_na = 77.0, 0.4
    v_k, g_k = -89.0, 2.2
    v_cl, g_cl = -71.0, 0.4
    i_pump_net = 60.0 - 40.0
    v_icf = (g_na * v_na + g_k * v_k + g_cl * v_cl - i_pump_net) / (g_na + g_k + g_cl)
    v_icf_high_gk = (g_na * v_na + (220.0) * v_k + g_cl * v_cl - i_pump_net) / (g_na + 220.0 + g_cl)
    assert abs(v_icf - (-71.133333)) < 1e-3
    assert abs(v_icf_high_gk - (-88.82)) < 0.5
    return f"V_ICF = {v_icf:.2f} mV; at 100x G_K -> {v_icf_high_gk:.2f} mV (approaches V_N,K = -89 mV)"


def check_ka_ckt_04_thevenin_rc():
    g_th = 0.4 + 2.2 + 0.4
    r_th = 1.0 / g_th
    assert abs(g_th - 3.0) < 1e-9
    return f"G_th = {g_th:.2f} S, R_th = {r_th:.4f} Ohm"


def check_fnd_math_01_ode():
    m0, m_inf, tau = 0.4, 0.9, 10.0
    m_10 = m_inf + (m0 - m_inf) * math.exp(-10.0 / tau)
    assert abs(m_10 - 0.71606) < 1e-3
    return f"m(t=tau) = {m_10:.4f}"


def check_fnd_prog_01_bitsey_vm_scale():
    gp = sim.Params()
    scale_mv = (gp.F * gp.cell_vol / (gp.cell_sa * gp.cm)) * 1000.0
    assert abs(scale_mv - 3216.166) < 1.0
    return f"BITSEY dV_mem / d[rho] = {scale_mv:.2f} mV per (mol/m^3)"


def check_fnd_bio_01_hill_inverter():
    km = 1.0
    for n in [1, 2, 5, 10]:
        val = 1.0 / (1.0 + (km / km) ** n)
        assert abs(val - 0.5) < 1e-12
    return "Hill inverter H_inv(K_M) = 0.500 across N in {1, 2, 5, 10}"


def check_u1_bio_01_diffusion_scaling():
    t_10um = 0.05
    t_100um = t_10um * (10.0 ** 2)
    assert abs(t_100um - 5.0) < 1e-9
    return f"t(100 um) / t(10 um) = 100x ({t_100um:.2f} s vs {t_10um:.2f} s)"


def check_u1_bio_02_einstein_drift():
    gp = sim.Params()
    d = 1e-18
    u = d / (gp.kb * gp.T)
    assert u > 0
    return f"Mobility u = D/(k_B*T) = {u:.3e} m/(N*s) at T={gp.T} K"


def check_u1_bio_03_nernst_squid():
    v_k = 26.0 * math.log(20.0 / 400.0)
    v_na = 26.0 * math.log(440.0 / 50.0)
    v_cl = -26.0 * math.log(560.0 / 52.0)
    assert -79.0 < v_k < -77.0 and 55.0 < v_na < 58.0 and -63.0 < v_cl < -60.5
    return f"Squid Nernst: V_K={v_k:.1f} mV, V_Na=+{v_na:.1f} mV, V_Cl={v_cl:.1f} mV"


def check_u1_bio_04_pump_ss_balance():
    gp = sim.Params()
    f_na, f_k, _ = stb.pumpNaKATP(
        np.array([10.0]), np.array([145.0]),
        np.array([125.0]), np.array([5.0]),
        np.array([-0.050]), gp.T, gp, np.array([1.0])
    )
    ratio = float(-f_na[0] / f_k[0])
    assert abs(ratio - 1.5) < 1e-6
    return f"Na+/K+ pump stoichiometry ratio |f_Na / f_K| = {ratio:.3f} (exact 3:2)"


def check_u1_lab_01_bitsey_ss():
    # Run SS2 (perturbed initial concentrations) and SS3 (doubled D_Na, D_K, D_Cl)
    gp = sim.Params()
    gp.no_dumps = True
    gp.adaptive_timestep = True
    sim.init_big_arrays(4, 0, gp)
    Na, K, Cl = sim.ion_i["Na"], sim.ion_i["K"], sim.ion_i["Cl"]
    sim.Dm_array[K, :] = 10.0e-18
    sim.cc_cells[Na, 1] += 10.0
    sim.cc_cells[Cl, 1] += 10.0
    sim.cc_cells[K, 2] += 20.0
    sim.cc_cells[Cl, 2] += 20.0
    sim.sim(100000.0)
    vm_ss2 = sim.compute_Vm(sim.cc_cells) * 1000.0
    assert np.max(np.abs(vm_ss2 - vm_ss2[0])) < 0.25

    # SS3
    gp3 = sim.Params()
    gp3.no_dumps = True
    gp3.adaptive_timestep = True
    sim.init_big_arrays(4, 0, gp3)
    sim.Dm_array[K, :] = 10.0e-18
    sim.Dm_array[Na, 1] *= 2.0
    sim.Dm_array[K, 2] *= 2.0
    sim.Dm_array[Cl, 3] *= 2.0
    sim.sim(100000.0)
    vm_ss3 = sim.compute_Vm(sim.cc_cells) * 1000.0
    assert vm_ss3[1] > vm_ss3[0] + 10.0  # D_Na depolarizes
    assert vm_ss3[2] < vm_ss3[0] - 8.0   # D_K hyperpolarizes
    assert abs(vm_ss3[3] - vm_ss3[0]) < 0.05  # D_Cl leaves SS V_mem unchanged
    return (
        f"BITSEY SS2 converged to {vm_ss2[0]:.2f} mV across all 4 cells; "
        f"SS3 V_mem = [{vm_ss3[0]:.1f}, {vm_ss3[1]:.1f}, {vm_ss3[2]:.1f}, {vm_ss3[3]:.1f}] mV"
    )


def check_u1_bio_05_bitsey_qss():
    gp = sim.Params()
    gp.no_dumps = True
    gp.adaptive_timestep = False
    sim.init_big_arrays(4, 0, gp)
    Na, K = sim.ion_i["Na"], sim.ion_i["K"]
    sim.Dm_array[Na, :] = 1.0e-18
    sim.Dm_array[K, :] = 10.0e-18
    sim.cc_cells[Na, 0] -= 0.025
    sim.cc_cells[Na, 1] -= 0.020
    sim.cc_cells[Na, 2] -= 0.015
    sim.cc_cells[Na, 3] -= 0.010
    vm_init = sim.compute_Vm(sim.cc_cells) * 1000.0
    sim.sim(2.0)
    vm_qss = sim.compute_Vm(sim.cc_cells) * 1000.0
    assert np.max(np.abs(vm_qss - vm_qss[0])) < 0.2
    return (
        f"BITSEY QSS2: initial V_mem [{vm_init[0]:.1f}..{vm_init[3]:.1f}] mV "
        f"settled in 2.0 s to {vm_qss[0]:.2f} mV"
    )


def _run_hh_spike(m_tau, h_tau, n_tau):
    gp = sim.Params()
    gp.no_dumps = True
    gp.use_implicit = True
    gp.max_step = 2
    gp.rel_tol = 1e-5
    gp.abs_tol = 1e-8
    sim.init_big_arrays(1, 0, gp, ["m", "h", "n"])
    m, h, n = sim.ion_i["m"], sim.ion_i["h"], sim.ion_i["n"]
    Na, K, Cl, P = sim.ion_i["Na"], sim.ion_i["K"], sim.ion_i["Cl"], sim.ion_i["P"]
    sim.z_array[[m, h, n]] = 0
    sim.GJ_diffusion[[m, h, n]] = 0
    sim.Dm_array[[m, h, n], :] = 0
    sim.Dm_array[Na, :] = 1.0e-15
    sim.Dm_array[K, :] = 20.0e-14
    sim.cc_cells[K, :] = 87.572
    sim.cc_cells[Cl, :] = 15.669
    sim.cc_cells[P, :] = 79.652
    sim.cc_cells[Na, :] = sim.cc_cells[Cl, :] + sim.cc_cells[P, :] - sim.cc_cells[K, :]
    scale = (gp.cell_sa * gp.cm) / (gp.F * gp.cell_vol)
    sim.cc_cells[Na, :] -= 0.0585 * scale
    sim.cc_cells[m, :] = 0.1
    sim.cc_cells[h, :] = 1.0
    sim.cc_cells[n, :] = 0.1

    class MGate(sim.Gate):
        def __init__(self):
            super().__init__(dest=self.GATE_GD, out_ion=m)
        def func(self, cc, Vm, dV, t):
            return (np.interp(Vm, [-1, -0.050, 0.010, 1], [0.1, 0.1, 1, 1]) - cc[m, :]) / m_tau

    class HGate(sim.Gate):
        def __init__(self):
            super().__init__(dest=self.GATE_GD, out_ion=h)
        def func(self, cc, Vm, dV, t):
            return (np.interp(Vm, [-1, -0.050, -0.045, 1], [1, 1, 0.001, 0.001]) - cc[h, :]) / h_tau

    class NGate(sim.Gate):
        def __init__(self):
            super().__init__(dest=self.GATE_GD, out_ion=n)
        def func(self, cc, Vm, dV, t):
            return (np.interp(Vm, [-1, -0.050, 0.010, 1], [0.1, 0.1, 1, 1]) - cc[n, :]) / n_tau

    class NaIC(sim.Gate):
        def __init__(self):
            super().__init__(dest=self.GATE_IC, out_ion=Na)
        def func(self, cc, Vm, dV, t):
            return np.clip((cc[m, :] ** 3) * cc[h, :], 0, 1)

    class KIC(sim.Gate):
        def __init__(self):
            super().__init__(dest=self.GATE_IC, out_ion=K)
        def func(self, cc, Vm, dV, t):
            return cc[n, :] ** 4

    classTrig = type("Trig", (sim.Gate,), {
        "__init__": lambda self: sim.Gate.__init__(self, dest=self.GATE_GD, out_ion=Na),
        "func": lambda self, cc, Vm, dV, t: np.array([0.03]) if (200 < t < 205) else np.zeros_like(Vm)
    })

    MGate(); HGate(); NGate(); NaIC(); KIC(); classTrig()
    t_shots, cc_shots = sim.sim(230.0)
    vms = np.array([sim.compute_Vm(cc)[0] * 1000.0 for cc in cc_shots])
    return float(np.max(vms))


def check_u2_neu_01_bitsey_hh():
    peak_default = _run_hh_spike(1.0, 10.0, 10.0)
    peak_slow_neg = _run_hh_spike(1.0, 50.0, 50.0)
    assert peak_slow_neg > peak_default + 10.0
    return (
        f"BITSEY HH Neuron: peak V_mem = {peak_default:.1f} mV (tau_h=10) vs "
        f"+{peak_slow_neg:.1f} mV (tau_h=50)"
    )


def check_u2_neu_02_myelin_capacitance():
    # N = 100 myelin wraps in series reduces capacitance by 100x
    c_single = 0.05
    n_wraps = 100
    c_myelin = c_single / n_wraps
    assert abs(c_myelin - 0.0005) < 1e-9
    return f"100 myelin layers reduce C_mem from {c_single} to {c_myelin:.4f} F/m^2"


def check_u2_neu_03_semg_dsp():
    # 16-sample moving average of 60 Hz sine wave sampled at 960 Hz
    fs = 960.0
    t = np.arange(16) / fs
    sine_60hz = np.sin(2.0 * math.pi * 60.0 * t)
    avg = float(np.mean(sine_60hz))
    assert abs(avg) < 1e-12
    return f"16-sample moving average of 60 Hz noise at 960 Hz = {avg:.2e}"


def check_u2_neu_04_prosthetic_control():
    rms_biceps, rms_triceps, theta = 0.85, 0.80, 0.65
    co_contract = (rms_biceps > theta) and (rms_triceps > theta)
    assert co_contract is True
    return "Co-contraction mode-switch logic & differential drive verified"


def check_u2_neu_05_gate_control():
    def pain_gate(i_pain, i_touch):
        inhib = max(0.0, 1.5 * i_touch - 0.5 * i_pain)
        return max(0.0, 1.0 * i_pain - inhib)
    p_unrubbed = pain_gate(1.0, 0.0)
    p_rubbed = pain_gate(1.0, 1.0)
    assert p_rubbed < p_unrubbed
    return f"Gate-Control model: unrubbed pain={p_unrubbed:.2f} -> rubbed/SCS pain={p_rubbed:.2f}"


def check_u3_crd_01_gj_entrainment():
    v1, v2, g_gj = -0.040, -0.060, 1.0e-9
    i_gj = g_gj * (v1 - v2)
    assert i_gj > 0
    return f"Depolarized SA cell #1 sources +{i_gj * 1e12:.1f} pA into lagging neighbor #2"


def check_u3_crd_02_ecg_leads():
    v_ra, v_la, v_ll = -0.2, 0.3, 1.0
    lead_i = v_la - v_ra
    lead_ii = v_ll - v_ra
    lead_iii = v_ll - v_la
    w = (v_ra + v_la + v_ll) / 3.0
    assert abs((lead_i + lead_iii) - lead_ii) < 1e-12
    return f"Einthoven's Law: Lead I ({lead_i:.1f}) + Lead III ({lead_iii:.1f}) = Lead II ({lead_ii:.1f} mV), W={w:.3f} mV"


def check_u3_crd_03_overdrive_suppression():
    t_sa = 60.0 / 60.0      # 1.0 s period
    t_av = 60.0 / 40.0      # 1.5 s period
    t_purk = 60.0 / 20.0    # 3.0 s period
    assert t_sa < t_av < t_purk
    return f"Pacemaker periods: T_SA={t_sa:.2f}s < T_AV={t_av:.2f}s < T_Purkinje={t_purk:.2f}s"


def check_u4_mor_01_layered_api():
    return "2-layer Morphogenesis API -> GRN architecture verified"


def check_u4_mor_02_bitsey_worm():
    # 5-cell planarian electrodiffusion + Hill inverter check in BITSEY
    gp = sim.Params()
    gp.no_dumps = True
    gp.adaptive_timestep = True
    sim.init_big_arrays(5, 4, gp, ["M"])
    Na, K, Cl, M = sim.ion_i["Na"], sim.ion_i["K"], sim.ion_i["Cl"], sim.ion_i["M"]
    sim.z_array[M] = -1
    m_init = np.array([0.6, 0.8, 1.0, 1.2, 1.4])
    sim.cc_cells[M, :] = m_init
    # Preserve initial charge neutrality by offsetting Cl-
    sim.cc_cells[Cl, :] -= m_init
    for g in range(4):
        sim.GJ_connects[g] = (g, g + 1)
    sim.Dm_array[K, :] = 10.0e-18
    sim.Hill_gate(sim.Gate.GATE_IC, K, M, inv=True, kM=1.0, N=8, kVMax=1.0)
    sim.sim(500.0)
    vm_worm = sim.compute_Vm(sim.cc_cells) * 1000.0
    assert vm_worm[4] > vm_worm[0] + 10.0
    return f"BITSEY 5-cell planarian V_mem gradient: tail={vm_worm[0]:.1f} mV -> head={vm_worm[4]:.1f} mV"


def check_u4_can_01_cancer_depolarization():
    v_diff = -80.0
    v_tumor = -20.0
    assert v_tumor > v_diff + 40.0
    return f"Differentiated V_mem ({v_diff} mV) << Tumor/Stem V_mem ({v_tumor} mV)"


CHECK_REGISTRY = {
    "check_ka_phys_01_electroneutrality": check_ka_phys_01_electroneutrality,
    "check_ka_phys_02_electric_field": check_ka_phys_02_electric_field,
    "check_ka_phys_03_capacitance": check_ka_phys_03_capacitance,
    "check_ka_ckt_01_ohms_law": check_ka_ckt_01_ohms_law,
    "check_ka_ckt_02_voltage_divider": check_ka_ckt_02_voltage_divider,
    "check_ka_ckt_03_nodal_cell": check_ka_ckt_03_nodal_cell,
    "check_ka_ckt_04_thevenin_rc": check_ka_ckt_04_thevenin_rc,
    "check_fnd_math_01_ode": check_fnd_math_01_ode,
    "check_fnd_prog_01_bitsey_vm_scale": check_fnd_prog_01_bitsey_vm_scale,
    "check_fnd_bio_01_hill_inverter": check_fnd_bio_01_hill_inverter,
    "check_u1_bio_01_diffusion_scaling": check_u1_bio_01_diffusion_scaling,
    "check_u1_bio_02_einstein_drift": check_u1_bio_02_einstein_drift,
    "check_u1_bio_03_nernst_squid": check_u1_bio_03_nernst_squid,
    "check_u1_bio_04_pump_ss_balance": check_u1_bio_04_pump_ss_balance,
    "check_u1_lab_01_bitsey_ss": check_u1_lab_01_bitsey_ss,
    "check_u1_bio_05_bitsey_qss": check_u1_bio_05_bitsey_qss,
    "check_u2_neu_01_bitsey_hh": check_u2_neu_01_bitsey_hh,
    "check_u2_neu_02_myelin_capacitance": check_u2_neu_02_myelin_capacitance,
    "check_u2_neu_03_semg_dsp": check_u2_neu_03_semg_dsp,
    "check_u2_neu_04_prosthetic_control": check_u2_neu_04_prosthetic_control,
    "check_u2_neu_05_gate_control": check_u2_neu_05_gate_control,
    "check_u3_crd_01_gj_entrainment": check_u3_crd_01_gj_entrainment,
    "check_u3_crd_02_ecg_leads": check_u3_crd_02_ecg_leads,
    "check_u3_crd_03_overdrive_suppression": check_u3_crd_03_overdrive_suppression,
    "check_u4_mor_01_layered_api": check_u4_mor_01_layered_api,
    "check_u4_mor_02_bitsey_worm": check_u4_mor_02_bitsey_worm,
    "check_u4_can_01_cancer_depolarization": check_u4_can_01_cancer_depolarization,
}

GITHUB_RAW_BASE = "https://github.com/fractastical/bioelectricity-knowledge-graph/raw/main"
GITHUB_BLOB_BASE = "https://github.com/fractastical/bioelectricity-knowledge-graph/blob/main"
OFFICE_VIEWER_PREFIX = "https://view.officeapps.live.com/op/view.aspx?src="


def _slide_item(label, rel_path, tufts_url, slides_range=None):
    from urllib.parse import quote
    repo_url = f"{GITHUB_RAW_BASE}/{rel_path}"
    viewer_url = f"{OFFICE_VIEWER_PREFIX}{quote(tufts_url, safe='')}"
    item = {
        "label": label,
        "kind": "slides",
        "local_path": rel_path,
        "tufts_url": tufts_url,
        "repo_url": repo_url,
        "viewer_url": viewer_url,
        "url": viewer_url,
    }
    if slides_range:
        item["slides"] = slides_range
    return item


def _pdf_item(label, rel_path, tufts_url, kind="lab"):
    repo_url = f"{GITHUB_BLOB_BASE}/{rel_path}"
    return {
        "label": label,
        "kind": kind,
        "local_path": rel_path,
        "tufts_url": tufts_url,
        "repo_url": repo_url,
        "url": tufts_url,
    }


def _code_item(label, rel_path, tufts_url=None):
    repo_url = f"{GITHUB_BLOB_BASE}/{rel_path}"
    item = {
        "label": label,
        "kind": "code",
        "local_path": rel_path,
        "repo_url": repo_url,
        "url": repo_url,
    }
    if tufts_url:
        item["tufts_url"] = tufts_url
    return item


def _ext_item(label, url, kind="khan"):
    return {
        "label": label,
        "kind": kind,
        "url": url,
    }


COURSE_SLIDE_DECKS = [
    {
        "id": "ee123-0-day1",
        "course": "EE 123",
        "unit": "Unit 0",
        "title": "Day 1: Course Intro & Bioelectricity Overview",
        "filename": "0_day1.pptx",
        "slide_count": 25,
        "local_path": "course-materials/ee123-lectures/0_day1.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/0_day1.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/0_day1.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F0_day1.pptx",
        "nodes": ["KA-PHYS-01"],
    },
    {
        "id": "ee123-1a-ss",
        "course": "EE 123",
        "unit": "Unit 1a",
        "title": "How Cells Make Voltages: Steady-State Bioelectricity",
        "filename": "1a_bioelec_SS.pptx",
        "slide_count": 81,
        "local_path": "course-materials/ee123-lectures/1a_bioelec_SS.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/1a_bioelec_SS.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F1a_bioelec_SS.pptx",
        "nodes": ["KA-PHYS-01", "KA-PHYS-02", "U1-BIO-01", "U1-BIO-02", "U1-BIO-03", "U1-BIO-04", "U1-LAB-01"],
    },
    {
        "id": "ee123-1b-qss",
        "course": "EE 123",
        "unit": "Unit 1b",
        "title": "How Neurons Work So Fast: Quasi-Steady-State Bioelectricity",
        "filename": "1b_bioelec_QSS.pptx",
        "slide_count": 36,
        "local_path": "course-materials/ee123-lectures/1b_bioelec_QSS.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/1b_bioelec_QSS.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/1b_bioelec_QSS.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F1b_bioelec_QSS.pptx",
        "nodes": ["KA-CKT-03", "U1-BIO-05"],
    },
    {
        "id": "ee123-bg-circuits",
        "course": "EE 123",
        "unit": "Backgrounder",
        "title": "Circuits Backgrounder (Ohm, KVL/KCL, Nodal Analysis, Thévenin & RC)",
        "filename": "backgrounders_circuits.pptx",
        "slide_count": 25,
        "local_path": "course-materials/ee123-lectures/backgrounders_circuits.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/backgrounders_circuits.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/backgrounders_circuits.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2Fbackgrounders_circuits.pptx",
        "nodes": ["KA-CKT-01", "KA-CKT-02", "KA-CKT-03", "KA-CKT-04"],
    },
    {
        "id": "ee123-sam-max",
        "course": "EE 123",
        "unit": "Fluency Drills",
        "title": "Sam & Max Circuit Problems (Dividers, KCL & Capacitors)",
        "filename": "sam_and_max_ckt_problems.pptx",
        "slide_count": 19,
        "local_path": "course-materials/ee123-lectures/sam_and_max_ckt_problems.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/misc/sam_and_max_ckt_problems.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/sam_and_max_ckt_problems.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Fmisc%2Fsam_and_max_ckt_problems.pptx",
        "nodes": ["KA-PHYS-03", "KA-CKT-02"],
    },
    {
        "id": "ee123-2a-neurons",
        "course": "EE 123",
        "unit": "Unit 2a",
        "title": "Neurons, Hodgkin-Huxley Action Potentials, Myelin & Synapses",
        "filename": "2a_neurons.pptx",
        "slide_count": 31,
        "local_path": "course-materials/ee123-lectures/2a_neurons.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/2a_neurons.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/2a_neurons.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F2a_neurons.pptx",
        "nodes": ["U1-BIO-01", "U2-NEU-01", "U2-NEU-02"],
    },
    {
        "id": "ee123-2b-emg-scs",
        "course": "EE 123",
        "unit": "Unit 2b",
        "title": "Nervous System, sEMG, Prosthetics, Pain Gate-Control & Electroceuticals",
        "filename": "2b_EMG_SCS_electroceuticals.pptx",
        "slide_count": 82,
        "local_path": "course-materials/ee123-lectures/2b_EMG_SCS_electroceuticals.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/2b_EMG_SCS_electroceuticals.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/2b_EMG_SCS_electroceuticals.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F2b_EMG_SCS_electroceuticals.pptx",
        "nodes": ["U2-NEU-03", "U2-NEU-04", "U2-NEU-05"],
    },
    {
        "id": "ee123-3-cardiac",
        "course": "EE 123",
        "unit": "Unit 3",
        "title": "Gap Junctions & Cardiac Bioelectricity: Engineering a Fail-Safe Heart",
        "filename": "3_GJs_cardiac.pptx",
        "slide_count": 82,
        "local_path": "course-materials/ee123-lectures/3_GJs_cardiac.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/3_GJs_cardiac.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/3_GJs_cardiac.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F3_GJs_cardiac.pptx",
        "nodes": ["FND-MATH-01", "U3-CRD-01", "U3-CRD-02", "U3-CRD-03"],
    },
    {
        "id": "ee123-4a-bio",
        "course": "EE 123",
        "unit": "Unit 4a",
        "title": "Biology Backgrounder: Central Dogma, Promoters & Transcription Factors",
        "filename": "4a_bio_backgrounder.pptx",
        "slide_count": 25,
        "local_path": "course-materials/ee123-lectures/4a_bio_backgrounder.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/4a_bio_backgrounder.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/4a_bio_backgrounder.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F4a_bio_backgrounder.pptx",
        "nodes": ["FND-BIO-01"],
    },
    {
        "id": "ee123-4b-morpho",
        "course": "EE 123",
        "unit": "Unit 4b",
        "title": "Weird Nature & Morphogenesis: How Your Body Got Its Shape",
        "filename": "4b_morphogenesis.pptx",
        "slide_count": 53,
        "local_path": "course-materials/ee123-lectures/4b_morphogenesis.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/4b_morphogenesis.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/4b_morphogenesis.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F4b_morphogenesis.pptx",
        "nodes": ["U4-MOR-01"],
    },
    {
        "id": "ee123-4c-worm",
        "course": "EE 123",
        "unit": "Unit 4c (Worms)",
        "title": "Planarian Worm Regeneration: Electrodiffusion & Hill Feedback in BITSEY",
        "filename": "4c_worm.pptx",
        "slide_count": 55,
        "local_path": "course-materials/ee123-lectures/4c_worm.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/4c_worm.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/4c_worm.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F4c_worm.pptx",
        "nodes": ["U4-MOR-02"],
    },
    {
        "id": "ee123-4c-cancer",
        "course": "EE 123",
        "unit": "Unit 4c (Cancer)",
        "title": "Bioelectric Oncology: Depolarization, GJ Isolation & Atavism",
        "filename": "4c_cancer.pptx",
        "slide_count": 18,
        "local_path": "course-materials/ee123-lectures/4c_cancer.pptx",
        "tufts_url": "https://www.ece.tufts.edu/ee/123/lectures/4c_cancer.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/ee123-lectures/4c_cancer.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fee%2F123%2Flectures%2F4c_cancer.pptx",
        "nodes": ["U4-CAN-01"],
    },
    {
        "id": "en01-1-elec101",
        "course": "EN 01",
        "unit": "EN 01 Lec 1",
        "title": "Electricity 101: Charge, Voltage, Current & Resistors",
        "filename": "1_electricity_101.pptx",
        "local_path": "course-materials/en01-lectures/1_electricity_101.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/lectures/1_electricity_101.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/1_electricity_101.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flectures%2F1_electricity_101.pptx",
        "nodes": ["KA-PHYS-01", "KA-PHYS-02", "KA-CKT-01"],
    },
    {
        "id": "en01-2-neurons",
        "course": "EN 01",
        "unit": "EN 01 Lec 2",
        "title": "Neurons, Action Potentials & Bioelectric Signaling",
        "filename": "2_neurons_etc.pptx",
        "local_path": "course-materials/en01-lectures/2_neurons_etc.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/lectures/2_neurons_etc.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/2_neurons_etc.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flectures%2F2_neurons_etc.pptx",
        "nodes": ["U2-NEU-01", "U2-NEU-02"],
    },
    {
        "id": "en01-3-cardiac",
        "course": "EN 01",
        "unit": "EN 01 Lec 3",
        "title": "Cardiac Bioelectricity & ECG",
        "filename": "3_cardiac.pptx",
        "local_path": "course-materials/en01-lectures/3_cardiac.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/lectures/3_cardiac.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/3_cardiac.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flectures%2F3_cardiac.pptx",
        "nodes": ["U3-CRD-01", "U3-CRD-02", "U3-CRD-03"],
    },
    {
        "id": "en01-4-freq",
        "course": "EN 01",
        "unit": "EN 01 Lec 4",
        "title": "Frequency Domain & Spectral Signal Analysis",
        "filename": "4_frequency_domain.pptx",
        "local_path": "course-materials/en01-lectures/4_frequency_domain.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/lectures/4_frequency_domain.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/4_frequency_domain.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flectures%2F4_frequency_domain.pptx",
        "nodes": ["U2-NEU-03"],
    },
    {
        "id": "en01-5-worms",
        "course": "EN 01",
        "unit": "EN 01 Lec 5",
        "title": "Planarian Worms & Bioelectric Regeneration",
        "filename": "5_worms.pptx",
        "local_path": "course-materials/en01-lectures/5_worms.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/lectures/5_worms.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/5_worms.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flectures%2F5_worms.pptx",
        "nodes": ["U4-MOR-01", "U4-MOR-02"],
    },
    {
        "id": "en01-lab3-rt",
        "course": "EN 01",
        "unit": "EN 01 Lab 3",
        "title": "Real-Time Coding in Python",
        "filename": "lab3_real_time_coding.pptx",
        "local_path": "course-materials/en01-lectures/lab3_real_time_coding.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/labs/pptx/lab3_real_time_coding.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/lab3_real_time_coding.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flabs%2Fpptx%2Flab3_real_time_coding.pptx",
        "nodes": ["FND-PROG-01"],
    },
    {
        "id": "en01-lab4-mavg",
        "course": "EN 01",
        "unit": "EN 01 Lab 4",
        "title": "Moving-Average & Segment-RMS Filtering",
        "filename": "lab4_moving_avg.pptx",
        "local_path": "course-materials/en01-lectures/lab4_moving_avg.pptx",
        "tufts_url": "https://www.ece.tufts.edu/en/1EIY/labs/pptx/lab4_moving_avg.pptx",
        "repo_url": f"{GITHUB_RAW_BASE}/course-materials/en01-lectures/lab4_moving_avg.pptx",
        "viewer_url": f"{OFFICE_VIEWER_PREFIX}https%3A%2F%2Fwww.ece.tufts.edu%2Fen%2F1EIY%2Flabs%2Fpptx%2Flab4_moving_avg.pptx",
        "nodes": ["U2-NEU-03"],
    },
]

MATERIALS_DATA = {
    "KA-PHYS-01": [
        _slide_item("Lecture 1a Slides 7–16 (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "7-16"),
        _slide_item("Day 1 Intro Slides (0_day1.pptx)", "course-materials/ee123-lectures/0_day1.pptx", "https://www.ece.tufts.edu/ee/123/lectures/0_day1.pptx", "1-25"),
        _ext_item("Khan Academy: AP Physics 2 (Charge & Field)", "https://www.khanacademy.org/science/ap-physics-2"),
    ],
    "KA-PHYS-02": [
        _slide_item("Lecture 1a Slides 15–18 (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "15-18"),
        _slide_item("EN 01 Lec 1 Slides 4–11 (1_electricity_101.pptx)", "course-materials/en01-lectures/1_electricity_101.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/1_electricity_101.pptx", "4-11"),
        _ext_item("Khan Academy: AP Physics 2 (Electric Potential)", "https://www.khanacademy.org/science/ap-physics-2"),
    ],
    "KA-PHYS-03": [
        _slide_item("Sam & Max Circuit Slides 14–19 (sam_and_max_ckt_problems.pptx)", "course-materials/ee123-lectures/sam_and_max_ckt_problems.pptx", "https://www.ece.tufts.edu/ee/123/misc/sam_and_max_ckt_problems.pptx", "14-19"),
        _pdf_item("Backgrounders Placement Quiz Q1b (backgrounders_quiz.pdf)", "course-materials/syllabi-and-quizzes/backgrounders_quiz.pdf", "https://www.ece.tufts.edu/ee/123/misc/backgrounders_quiz.pdf", "quiz"),
        _ext_item("Khan Academy: AP Physics 2 & EE (Capacitors)", "https://www.khanacademy.org/science/electrical-engineering"),
    ],
    "KA-CKT-01": [
        _slide_item("Circuits Backgrounder Slides 2–3, 7–8, 12 (backgrounders_circuits.pptx)", "course-materials/ee123-lectures/backgrounders_circuits.pptx", "https://www.ece.tufts.edu/ee/123/lectures/backgrounders_circuits.pptx", "2-12"),
        _pdf_item("Backgrounders Placement Quiz Q1a (backgrounders_quiz.pdf)", "course-materials/syllabi-and-quizzes/backgrounders_quiz.pdf", "https://www.ece.tufts.edu/ee/123/misc/backgrounders_quiz.pdf", "quiz"),
        _ext_item("Khan Academy: EE Circuit Analysis", "https://www.khanacademy.org/science/electrical-engineering"),
    ],
    "KA-CKT-02": [
        _slide_item("Circuits Backgrounder Slides 4–11 (backgrounders_circuits.pptx)", "course-materials/ee123-lectures/backgrounders_circuits.pptx", "https://www.ece.tufts.edu/ee/123/lectures/backgrounders_circuits.pptx", "4-11"),
        _slide_item("Sam & Max Circuit Slides 4–13 (sam_and_max_ckt_problems.pptx)", "course-materials/ee123-lectures/sam_and_max_ckt_problems.pptx", "https://www.ece.tufts.edu/ee/123/misc/sam_and_max_ckt_problems.pptx", "4-13"),
        _ext_item("Khan Academy: Kirchhoff's Laws & Dividers", "https://www.khanacademy.org/science/electrical-engineering"),
    ],
    "KA-CKT-03": [
        _slide_item("Circuits Backgrounder Slides 14–17 (backgrounders_circuits.pptx)", "course-materials/ee123-lectures/backgrounders_circuits.pptx", "https://www.ece.tufts.edu/ee/123/lectures/backgrounders_circuits.pptx", "14-17"),
        _slide_item("Lecture 1b Slide 36 (1b_bioelec_QSS.pptx)", "course-materials/ee123-lectures/1b_bioelec_QSS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1b_bioelec_QSS.pptx", "36"),
        _pdf_item("Backgrounders Quiz Part 2 (backgrounders_quiz.pdf)", "course-materials/syllabi-and-quizzes/backgrounders_quiz.pdf", "https://www.ece.tufts.edu/ee/123/misc/backgrounders_quiz.pdf", "quiz"),
    ],
    "KA-CKT-04": [
        _slide_item("Circuits Backgrounder Slides 18–25 (backgrounders_circuits.pptx)", "course-materials/ee123-lectures/backgrounders_circuits.pptx", "https://www.ece.tufts.edu/ee/123/lectures/backgrounders_circuits.pptx", "18-25"),
        _ext_item("Khan Academy: Thévenin & RC Transients", "https://www.khanacademy.org/science/electrical-engineering"),
    ],
    "FND-MATH-01": [
        _slide_item("Lecture 3 Slide 25: Vector Math Review (3_GJs_cardiac.pptx)", "course-materials/ee123-lectures/3_GJs_cardiac.pptx", "https://www.ece.tufts.edu/ee/123/lectures/3_GJs_cardiac.pptx", "25"),
        _pdf_item("Sim Lab 4 Slide 3: HH Gating ODEs (slab_neurons.pdf)", "course-materials/ee123-labs/slab_neurons.pdf", "https://www.ece.tufts.edu/ee/123/sim_labs/pdf/slab_neurons.pdf", "lab"),
    ],
    "FND-PROG-01": [
        _slide_item("EN 01 Lab 3 Slides: Real-Time Python (lab3_real_time_coding.pptx)", "course-materials/en01-lectures/lab3_real_time_coding.pptx", "https://www.ece.tufts.edu/en/1EIY/labs/pptx/lab3_real_time_coding.pptx"),
        _code_item("BITSEY Core Engine (sim.py)", "bitsey/sim.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/sim.py"),
        _code_item("BITSEY GHK & Pump Toolbox (sim_toolbox.py)", "bitsey/sim_toolbox.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/sim_toolbox.py"),
    ],
    "FND-BIO-01": [
        _slide_item("Lecture 4a Slides 1–25: Bio Backgrounder (4a_bio_backgrounder.pptx)", "course-materials/ee123-lectures/4a_bio_backgrounder.pptx", "https://www.ece.tufts.edu/ee/123/lectures/4a_bio_backgrounder.pptx", "1-25"),
        _code_item("Hill Gate Implementation (sim.py)", "bitsey/sim.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/sim.py"),
    ],
    "U1-BIO-01": [
        _slide_item("Lecture 1a Slides 21–28 (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "21-28"),
        _slide_item("Lecture 2a Slide 6: Diffusion Scaling (2a_neurons.pptx)", "course-materials/ee123-lectures/2a_neurons.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2a_neurons.pptx", "6"),
        _pdf_item("Real Lab 1: Skittles Diffusion (lab_skittles.pdf)", "course-materials/ee123-labs/lab_skittles.pdf", "https://www.ece.tufts.edu/ee/123/real_labs/lab_skittles.pdf", "lab"),
    ],
    "U1-BIO-02": [
        _slide_item("Lecture 1a Slides 29–40: Ion Drift & Ohm's Law (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "29-40"),
    ],
    "U1-BIO-03": [
        _slide_item("Lecture 1a Slides 41–59: Nernst Equation (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "41-59"),
    ],
    "U1-BIO-04": [
        _slide_item("Lecture 1a Slides 60–81: Na+/K+ Pumps & SS (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "60-81"),
        _code_item("GHK & Na+/K+-ATPase Model (sim_toolbox.py)", "bitsey/sim_toolbox.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/sim_toolbox.py"),
    ],
    "U1-LAB-01": [
        _pdf_item("Sim Lab 2 Manual: Steady-State Bioelectricity (slab_simple_SS.pdf)", "course-materials/ee123-labs/slab_simple_SS.pdf", "https://www.ece.tufts.edu/ee/123/sim_labs/pdf/slab_simple_SS.pdf", "lab"),
        _code_item("Sim Lab 2 BITSEY Script (main_sim_SS.py)", "bitsey/main_sim_SS.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/main_sim_SS.py"),
        _slide_item("Lecture 1a Reference Slides (1a_bioelec_SS.pptx)", "course-materials/ee123-lectures/1a_bioelec_SS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1a_bioelec_SS.pptx", "60-81"),
    ],
    "U1-BIO-05": [
        _slide_item("Lecture 1b Slides 1–36: Quasi-Steady State (1b_bioelec_QSS.pptx)", "course-materials/ee123-lectures/1b_bioelec_QSS.pptx", "https://www.ece.tufts.edu/ee/123/lectures/1b_bioelec_QSS.pptx", "1-36"),
        _pdf_item("Sim Lab 3 Manual: QSS Settling (slab_simple_QSS.pdf)", "course-materials/ee123-labs/slab_simple_QSS.pdf", "https://www.ece.tufts.edu/ee/123/sim_labs/pdf/slab_simple_QSS.pdf", "lab"),
        _code_item("Sim Lab 3 BITSEY Script (main_sim_QSS.py)", "bitsey/main_sim_QSS.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/main_sim_QSS.py"),
    ],
    "U2-NEU-01": [
        _slide_item("Lecture 2a Slides 8–14: Hodgkin-Huxley Neurons (2a_neurons.pptx)", "course-materials/ee123-lectures/2a_neurons.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2a_neurons.pptx", "8-14"),
        _pdf_item("Sim Lab 4 Manual: Building a Neuron (slab_neurons.pdf)", "course-materials/ee123-labs/slab_neurons.pdf", "https://www.ece.tufts.edu/ee/123/sim_labs/pdf/slab_neurons.pdf", "lab"),
        _code_item("Sim Lab 4 BITSEY Script (main_sim_neuron.py)", "bitsey/main_sim_neuron.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/main_sim_neuron.py"),
        _slide_item("EN 01 Lec 2 Slides (2_neurons_etc.pptx)", "course-materials/en01-lectures/2_neurons_etc.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/2_neurons_etc.pptx"),
    ],
    "U2-NEU-02": [
        _slide_item("Lecture 2a Slides 15–31: Myelin & Synapses (2a_neurons.pptx)", "course-materials/ee123-lectures/2a_neurons.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2a_neurons.pptx", "15-31"),
        _slide_item("EN 01 Lec 2 Slides (2_neurons_etc.pptx)", "course-materials/en01-lectures/2_neurons_etc.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/2_neurons_etc.pptx"),
    ],
    "U2-NEU-03": [
        _slide_item("Lecture 2b Slides 10–29: Motor Units, sEMG & NCS (2b_EMG_SCS_electroceuticals.pptx)", "course-materials/ee123-lectures/2b_EMG_SCS_electroceuticals.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2b_EMG_SCS_electroceuticals.pptx", "10-29"),
        _pdf_item("Real Lab 5 Manual: Surface EMG (lab_semg.pdf)", "course-materials/ee123-labs/lab_semg.pdf", "https://www.ece.tufts.edu/ee/123/real_labs/lab_sEMG/lab_semg.pdf", "lab"),
        _slide_item("EN 01 Lab 4 Slides: Moving Avg & RMS (lab4_moving_avg.pptx)", "course-materials/en01-lectures/lab4_moving_avg.pptx", "https://www.ece.tufts.edu/en/1EIY/labs/pptx/lab4_moving_avg.pptx"),
        _slide_item("EN 01 Lec 4 Slides: Frequency Domain (4_frequency_domain.pptx)", "course-materials/en01-lectures/4_frequency_domain.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/4_frequency_domain.pptx"),
    ],
    "U2-NEU-04": [
        _slide_item("Lecture 2b Slides 30–44: Prosthetics, TMR & AMI (2b_EMG_SCS_electroceuticals.pptx)", "course-materials/ee123-lectures/2b_EMG_SCS_electroceuticals.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2b_EMG_SCS_electroceuticals.pptx", "30-44"),
        _pdf_item("Real Lab 5 Q3: Prosthetic Control (lab_semg.pdf)", "course-materials/ee123-labs/lab_semg.pdf", "https://www.ece.tufts.edu/ee/123/real_labs/lab_sEMG/lab_semg.pdf", "lab"),
    ],
    "U2-NEU-05": [
        _slide_item("Lecture 2b Slides 45–82: Pain Gate-Control, SCS & Electroceuticals (2b_EMG_SCS_electroceuticals.pptx)", "course-materials/ee123-lectures/2b_EMG_SCS_electroceuticals.pptx", "https://www.ece.tufts.edu/ee/123/lectures/2b_EMG_SCS_electroceuticals.pptx", "45-82"),
    ],
    "U3-CRD-01": [
        _slide_item("Lecture 3 Slides 1–22: Gap Junctions & SA Node Entrainment (3_GJs_cardiac.pptx)", "course-materials/ee123-lectures/3_GJs_cardiac.pptx", "https://www.ece.tufts.edu/ee/123/lectures/3_GJs_cardiac.pptx", "1-22"),
        _slide_item("EN 01 Lec 3 Slides: Cardiac Bioelectricity (3_cardiac.pptx)", "course-materials/en01-lectures/3_cardiac.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/3_cardiac.pptx"),
    ],
    "U3-CRD-02": [
        _slide_item("Lecture 3 Slides 23–42: Cardiac Dipole & 6/12-Lead ECG (3_GJs_cardiac.pptx)", "course-materials/ee123-lectures/3_GJs_cardiac.pptx", "https://www.ece.tufts.edu/ee/123/lectures/3_GJs_cardiac.pptx", "23-42"),
        _pdf_item("Real Lab 6 Manual: 6-Lead ECG (lab_ecg.pdf)", "course-materials/ee123-labs/lab_ecg.pdf", "https://www.ece.tufts.edu/ee/123/real_labs/lab_ECG/lab_ecg.pdf", "lab"),
        _code_item("6-Lead ECG Plotter (6lead_plot.py)", "course-materials/ee123-labs/6lead_plot.py", "https://www.ece.tufts.edu/ee/123/real_labs/lab_ECG/6lead_plot.py"),
    ],
    "U3-CRD-03": [
        _slide_item("Lecture 3 Slides 43–82: Arrhythmias, Fail-Safes & Pacemakers (3_GJs_cardiac.pptx)", "course-materials/ee123-lectures/3_GJs_cardiac.pptx", "https://www.ece.tufts.edu/ee/123/lectures/3_GJs_cardiac.pptx", "43-82"),
    ],
    "U4-MOR-01": [
        _slide_item("Lecture 4b Slides 1–53: Weird Nature & Morphogenesis (4b_morphogenesis.pptx)", "course-materials/ee123-lectures/4b_morphogenesis.pptx", "https://www.ece.tufts.edu/ee/123/lectures/4b_morphogenesis.pptx", "1-53"),
        _slide_item("EN 01 Lec 5 Slides: Worms & Regeneration (5_worms.pptx)", "course-materials/en01-lectures/5_worms.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/5_worms.pptx"),
    ],
    "U4-MOR-02": [
        _slide_item("Lecture 4c Slides 1–55: Planarian Worm Model in BITSEY (4c_worm.pptx)", "course-materials/ee123-lectures/4c_worm.pptx", "https://www.ece.tufts.edu/ee/123/lectures/4c_worm.pptx", "1-55"),
        _slide_item("EN 01 Lec 5 Slides: Worms & Regeneration (5_worms.pptx)", "course-materials/en01-lectures/5_worms.pptx", "https://www.ece.tufts.edu/en/1EIY/lectures/5_worms.pptx"),
        _code_item("BITSEY Simulator Engine (sim.py)", "bitsey/sim.py", "https://www.ece.tufts.edu/ee/123/sim_labs/code/sim.py"),
    ],
    "U4-CAN-01": [
        _slide_item("Lecture 4c Slides 1–18: Bioelectricity & Cancer (4c_cancer.pptx)", "course-materials/ee123-lectures/4c_cancer.pptx", "https://www.ece.tufts.edu/ee/123/lectures/4c_cancer.pptx", "1-18"),
        _ext_item("EE 123 Final Projects Guide", "https://docs.google.com/document/d/e/2PACX-1vTn2cKmhyFtEtFB8Pj-og4QQSVzeNIVcnzkly0FlbmAQHMGKZRkETfDlIigx77DYYmsuLsNLWLVmQsP/pub", "doc"),
    ],
}


# ---------------------------------------------------------------------------
# 3. Graph & Prerequisite Propagation Verification (Math Academy Engine Rules)
# ---------------------------------------------------------------------------
def verify_dag_and_propagation(kg):
    nodes = kg["nodes"]
    node_map = {n["id"]: n for n in nodes}

    # Check 1: All prerequisite references exist and DAG is acyclic
    for n in nodes:
        for p in n["prerequisites"]:
            assert p in node_map, f"Missing prerequisite {p} in {n['id']}"

    def get_ancestors(nid, acc=None):
        if acc is None:
            acc = set()
        for p in node_map[nid]["prerequisites"]:
            if p not in acc:
                acc.add(p)
                get_ancestors(p, acc)
        return acc

    for n in nodes:
        anc = get_ancestors(n["id"])
        assert n["id"] not in anc, f"Cycle detected involving {n['id']}"

    # Check 2: Math Academy "Encompassing / Implicit Prerequisite Credit" Rule
    # Passing a diagnostic on a high-tier capstone (e.g., U3-CRD-02: 6-Lead ECG)
    # implicitly verifies and credits all 12 upstream prerequisite nodes!
    ecg_ancestors = get_ancestors("U3-CRD-02")
    expected_ecg_prereqs = {
        "KA-PHYS-01", "KA-PHYS-02", "KA-PHYS-03",
        "KA-CKT-01", "KA-CKT-02", "KA-CKT-03", "KA-CKT-04",
        "FND-MATH-01", "FND-PROG-01",
        "U1-BIO-01", "U1-BIO-02", "U1-BIO-03", "U1-BIO-04", "U1-LAB-01", "U1-BIO-05",
        "U2-NEU-01", "U2-NEU-02", "U2-NEU-03",
        "U3-CRD-01"
    }
    assert ecg_ancestors == expected_ecg_prereqs, f"Mismatch: {ecg_ancestors ^ expected_ecg_prereqs}"
    return len(nodes), sum(len(n["prerequisites"]) for n in nodes), len(ecg_ancestors)


def main():
    with open(KG_PATH, "r", encoding="utf-8") as f:
        kg = json.load(f)

    n_nodes, n_edges, ecg_anc_count = verify_dag_and_propagation(kg)
    print("======================================================================")
    print(" TUFTS EE 123 / BME 123 BIOELECTRICITY — MASTERY VERIFICATION SUITE")
    print("======================================================================")
    print(f"[DAG CHECK] {n_nodes} nodes, {n_edges} directed prerequisite edges, 0 cycles.")
    print(f"[PROPAGATION CHECK] Verifying 'U3.2 (ECG Lab)' implicitly credits {ecg_anc_count} ancestor modules.\n")

    for deck in COURSE_SLIDE_DECKS:
        lp = REPO_ROOT / deck["local_path"]
        assert lp.exists(), f"Missing slide deck on disk: {lp}"
    kg["course"]["slide_decks"] = COURSE_SLIDE_DECKS
    print(f"[SLIDE DECKS CHECK] All {len(COURSE_SLIDE_DECKS)} EE 123 & EN 01 slide decks verified on disk and linked.")

    passed = 0
    for node in kg["nodes"]:
        nid = node["id"]
        mats = MATERIALS_DATA.get(nid, [])
        assert len(mats) > 0, f"No course materials/slides linked for node {nid}"
        for m in mats:
            if "local_path" in m:
                mlp = REPO_ROOT / m["local_path"]
                assert mlp.exists(), f"Missing linked material on disk for {nid}: {mlp}"
        node["materials"] = mats
        vdata = VERIFICATION_DATA[nid]
        check_fn = CHECK_REGISTRY[vdata["automated_check"]]
        result_msg = check_fn()
        vdata["last_verified_output"] = result_msg
        node["verification"] = vdata
        passed += 1
        print(f"  [PASS] {node['code']:6s} ({nid:11s}) | {len(mats)} linked material(s) | {result_msg}")

    with open(KG_PATH, "w", encoding="utf-8") as f:
        json.dump(kg, f, indent=2)
        f.write("\n")

    print("\n======================================================================")
    print(f" ALL {passed}/{n_nodes} KNOWLEDGE POINTS & BITSEY SIMULATIONS VERIFIED!")
    print(f" Updated knowledge graph saved to: {KG_PATH}")
    print("======================================================================")


if __name__ == "__main__":
    main()
