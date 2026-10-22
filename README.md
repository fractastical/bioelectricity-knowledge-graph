# Tufts EE 123 / BME 123: Bioelectricity — Math Academy-Style Knowledge Graph

A complete module extraction, **Math Academy-style Prerequisite Directed Acyclic Graph (DAG)**, interactive web explorer, and `BITSEY` simulation archive for **Prof. Joel Grodstein's Bioelectricity curriculum at Tufts University** ([`EE 123 / BME 123: Bioelectricity`](https://www.ece.tufts.edu/ee/123/) and its companion hardware/DSP freshman seminar [`EN 01: Electricity Inside of You`](https://www.ece.tufts.edu/en/1EIY/)), including all **Khan Academy Physics & Electrical Engineering prerequisite modules** assigned in the course backgrounders.

---

## Repository Structure

```text
bioelectricity-knowledge-graph/
├── knowledge-graph/
│   └── knowledge_graph.json          # 27-node Math Academy-style DAG, with a diagnostic on every node
├── docs/
│   └── index.html                    # Standalone DAG explorer; mastery requires the diagnostic
├── verify_mastery.py                 # DAG, numerical, and live BITSEY checks for all 27 nodes
├── bitsey/                           # BITSEY Python bioelectric network simulator
│   ├── sim.py                        # Core vectorized cell, ion channel, pump & GJ engine
│   ├── sim_toolbox.py                # Goldman-Hodgkin-Katz (GHK) & Na+/K+-ATPase pump models
│   ├── edebug.py                     # Equivalent-circuit analyzer & flux debugging utilities
│   ├── eplot.py                      # V_mem and ion concentration plotting helpers
│   ├── main_sim_SS.py                # Sim Lab 2: Steady-State (SS) single-cell bioelectricity
│   ├── main_sim_QSS.py               # Sim Lab 3: Quasi-Steady-State (QSS) fast settling
│   └── main_sim_neuron.py            # Sim Lab 4: Hodgkin-Huxley (m^3h, n^4) action potentials
├── course-materials/
│   ├── syllabi-and-quizzes/          # EE 123 & EN 01 syllabi + Khan Academy placement quiz PDF
│   ├── ee123-lectures/               # All 12 EE 123 PowerPoint slide decks (Units 0 through 4c)
│   ├── ee123-labs/                   # All 6 lab manuals (PDFs), ECG/sEMG sample data & scripts
│   └── en01-lectures/                # EN 01 Electricity 101, Frequency Domain, Moving Avg & DSP decks
└── requirements.txt
```

---

## Running the Interactive Knowledge Graph Explorer

The explorer is published at <https://fractastical.github.io/bioelectricity-knowledge-graph/>.

```bash
open docs/index.html
```

Features included in `docs/index.html`:
- **Interactive SVG Prerequisite DAG (Tiers 0–6)**: Click any node to highlight its upstream prerequisite chain and downstream unlocked modules.
- **Math Academy Learning Frontier**: A module unlocks only after every prerequisite is verified. Answering its diagnostic correctly records the XP. Revoking a module also clears every module that depends on it.
- **Khan Academy Bridge & Placement Quiz**: Full mapping of Prof. Grodstein's assigned Khan Academy `AP Physics 2` and `Electrical Engineering` units alongside the solutions to `backgrounders_quiz.pdf`.
- **Module Inspector**: Four-choice diagnostic, equations, and the last automated check result for every node.

Verified progress is stored in the browser under `ee123_verified`.

## Automated verification

`verify_mastery.py` checks the prerequisite DAG, the numerical and circuit diagnostics, and live BITSEY runs (steady state, quasi-steady state, Hodgkin–Huxley, and a 5-cell planarian gradient). It writes each result back onto the matching node in `knowledge_graph.json`.

```bash
python3 verify_mastery.py
```

---

## Running the `BITSEY` Bioelectric Simulator Labs

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
cd bitsey

# Sim Lab 2: Steady-State Bioelectricity (100,000 seconds virtual time)
python3 main_sim_SS.py SS1 100000
python3 main_sim_SS.py SS2 100000
python3 main_sim_SS.py SS3 100000

# Sim Lab 3: Quasi-Steady-State (QSS) Fast Capacitive Settling (2 seconds virtual time)
python3 main_sim_QSS.py QSS1 2

# Sim Lab 4: Hodgkin-Huxley Neuron Action Potential (230 seconds virtual time)
python3 main_sim_neuron.py HH 230
```

---

## Curriculum Tracks & Knowledge Graph Nodes (27 Nodes, 795 Total XP)

### 1. `KHAN_PHYS` — Khan Academy: AP Physics 2 (Electrostatics & Fields)
- **`KA.P1` (`KA-PHYS-01`, Tier 0, 15 XP)**: Electric Charge, Valence & Coulomb Forces
- **`KA.P2` (`KA-PHYS-02`, Tier 1, 15 XP)**: Electric Potential, Work & Relative Voltage
- **`KA.P3` (`KA-PHYS-03`, Tier 1, 20 XP)**: Electric Current & Membrane Capacitance

### 2. `KHAN_EE` — Khan Academy: Electrical Engineering & Circuit Analysis
- **`KA.C1` (`KA-CKT-01`, Tier 1, 20 XP)**: Ohm's Law, Conductance & Circuit Elements
- **`KA.C2` (`KA-CKT-02`, Tier 2, 25 XP)**: KVL, KCL, Series/Parallel & Voltage/Current Dividers
- **`KA.C3` (`KA-CKT-03`, Tier 3, 30 XP)**: Nodal Analysis of Multi-Branch Cell Circuits
- **`KA.C4` (`KA-CKT-04`, Tier 3, 25 XP)**: Thévenin/Norton Equivalents & $RC$ Delay Dynamics

### 3. `FOUNDATIONS` — Math, Python (`BITSEY`) & Molecular Biology Foundations
- **`FND.M1` (`FND-MATH-01`, Tier 0, 20 XP)**: Exponential Relaxation ODEs & 2D Vector Projections
- **`FND.CS1` (`FND-PROG-01`, Tier 1, 25 XP)**: Python, NumPy Vectorization & `BITSEY` Architecture
- **`FND.B1` (`FND-BIO-01`, Tier 0, 20 XP)**: Central Dogma, Transcription Factors & Gene Logic

### 4. `UNIT_1` — Cellular Bioelectricity (Steady State & Quasi-Steady State)
- **`U1.1` (`U1-BIO-01`, Tier 1, 25 XP)**: Fick's Law of Diffusion & Skittles Lab (`lab_skittles.pdf`)
- **`U1.2` (`U1-BIO-02`, Tier 2, 25 XP)**: Ion Drift Flux, Mobility & Membrane Ohm's Law (`1a_bioelec_SS.pptx`)
- **`U1.3` (`U1-BIO-03`, Tier 3, 35 XP)**: Semi-Permeable Membranes & The Nernst Equation (`1a_bioelec_SS.pptx`)
- **`U1.4` (`U1-BIO-04`, Tier 4, 35 XP)**: $\text{Na}^+/\text{K}^+$-ATPase Pumps, Detailed Balance vs. Steady State (`1a_bioelec_SS.pptx`)
- **`U1.L2` (`U1-LAB-01`, Tier 5, 40 XP)**: Sim Lab 2: Steady-State Bioelectricity in `BITSEY` (`slab_simple_SS.pdf`)
- **`U1.5` (`U1-BIO-05`, Tier 5, 40 XP)**: Quasi-Steady State (QSS) & Fast Capacitive Settling (`1b_bioelec_QSS.pptx`, `slab_simple_QSS.pdf`, Oral Quiz #1)

### 5. `UNIT_2` — Neurons, sEMG, Prosthetics, Pain & Electroceuticals
- **`U2.1` (`U2-NEU-01`, Tier 6, 45 XP)**: Hodgkin-Huxley Action Potentials & Gating ODEs (`2a_neurons.pptx`, `slab_neurons.pdf`)
- **`U2.2` (`U2-NEU-02`, Tier 6, 35 XP)**: Axonal Propagation, Myelin Capacitance & Synapses (`2a_neurons.pptx`)
- **`U2.3` (`U2-NEU-03`, Tier 6, 40 XP)**: Motor Units, sEMG, Nerve Conduction (NCS) & DSP Filtering (`2b_EMG_SCS_electroceuticals.pptx`, `lab_semg.pdf`, `lab4_moving_avg.pptx`)
- **`U2.4` (`U2-NEU-04`, Tier 6, 30 XP)**: Myoelectric Prosthetics, TMR & AMI Interfaces (`2b_EMG_SCS_electroceuticals.pptx`)
- **`U2.5` (`U2-NEU-05`, Tier 6, 35 XP)**: Gate-Control Theory of Pain, SCS, PNS & Electroceuticals (`2b_EMG_SCS_electroceuticals.pptx`, Oral Quiz #2)

### 6. `UNIT_3` — Gap Junctions, Cardiac Bioelectricity & Fail-Safe Heart
- **`U3.1` (`U3-CRD-01`, Tier 5, 35 XP)**: Cardiac Conduction System, Funny Channels & Gap Junctions (`3_GJs_cardiac.pptx`)
- **`U3.2` (`U3-CRD-02`, Tier 6, 45 XP)**: ECG Physics, Dipole Vector Projections & 6/12-Lead Geometry (`3_GJs_cardiac.pptx`, `lab_ecg.pdf`)
- **`U3.3` (`U3-CRD-03`, Tier 6, 40 XP)**: Cardiac Arrhythmias, Triple-Redundancy Fail-Safes & Pacemakers (`3_GJs_cardiac.pptx`, Oral Quiz #3)

### 7. `UNIT_4` — Developmental Bioelectricity, Morphogenesis & Cancer
- **`U4.1` (`U4-MOR-01`, Tier 5, 35 XP)**: Morphogenesis as a Layered Bioelectric System (`4b_morphogenesis.pptx`)
- **`U4.2` (`U4-MOR-02`, Tier 6, 45 XP)**: Planarian Polarity Modeling: Electrodiffusion, Hill Feedback & GJ Scaling (`4c_worm.pptx`)
- **`U4.3` (`U4-CAN-01`, Tier 6, 35 XP)**: Bioelectric Oncology: Depolarization, GJ Isolation & Atavism (`4c_cancer.pptx`)

---

## Credits & Attribution

- **Course Author & Instructor**: [Prof. Joel Grodstein](https://www.eecs.tufts.edu/~joelg/) (Department of Electrical and Computer Engineering, Tufts University).
- **`BITSEY` Simulator**: Copyright Alexis Pietak, Cecil Curry, and Joel Grodstein (adapted from `BETSE`, Allen Discovery Center at Tufts University).
- **Khan Academy Refresher Modules**: [Khan Academy AP Physics 2](https://www.khanacademy.org/science/ap-physics-2) & [Khan Academy Electrical Engineering](https://www.khanacademy.org/science/electrical-engineering).
