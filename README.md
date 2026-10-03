# Load Identification in Composite Plates Using Distributed Strain Measurements

**Towards the Identification of Moving Loads from Deformation Measures Distributed on Composite Plates**

**Author:** Gabriele Venturi  
**Supervisor:** Assistant Professor Giancarlo Santamato  
**Institution:** University of Pisa

## 📋 Project Overview

This project addresses the **identification of dynamic moving loads** on composite laminated structures using distributed Fiber Bragg Grating (FBG) sensors. The work combines:

- **Finite Element modeling** of composite plates
- **Analytical models** for load propagation
- **Data-driven approaches** for system identification

### Key Features
- 🔧 **8 FBG sensors** embedded in composite laminate (oriented at 0°, 90°, ±45°)
- 📊 **Moving load model** based on pressure wave propagation
- 🔬 **Multiple scenarios** tested (3 different load motions)
- 📈 **Modal & transient analysis** capabilities

---

## 🎯 Objectives

### Phase 1: Forward Problem (FEM Modeling)
1. Develop a dynamic FEM model of the composite plate
2. Characterize modal properties (frequencies, mode shapes)
3. Simulate structural response under moving loads
4. Extract strain measurements at FBG sensor locations

### Phase 2: Inverse Problem (Load Identification)
1. Formulate and solve the inverse problem
2. Recover load parameters from strain measurements:
   - Load position (initial x₀, y₀)
   - Load velocity (u, v)
   - Load amplitude (pmax)
3. Validate with synthetic and experimental data

---

## 🏗️ Project Structure

```
Load_Identificaiton_on_Composite_Plates/
├── Models/                      # CAD geometries (SLDASM, STEP, x_t, IGS)
├── Documentation/               # Project documentation (PDF)
├── Scripts/                     # ANSYS APDL scripts (preprocessing, solving, post-processing)
├── Analytical/                  # Python scripts for ROM problem
└── README.md                    # This file
```

---

## 🔧 ANSYS Workflow (Scripts/)

### Step 1: Geometry & Mesh Setup
```bash
# Open ANSYS APDL 
# Run in Workbench's built-in script editor or from ANSYS console:
/INPUT, geometry_loader.apdl
```
**Function:** Imports the composite laminate geometry, defines element types, material properties, and generates the FEM mesh.

**Key parameters:**
- Element size (adjustable for convergence studies)
- Material properties (orthotropic laminate)
- Boundary conditions (clamped edges)

---

### Step 2: Coordinate System Definition
```bash
/INPUT, SYS_orienter.apdl
```
**Function:** Creates local coordinate systems on each FEM element aligned with ply axes.

**Why it matters:** 
- Composite Materials are orthotropic and need proper orientation
- Supports multi-directional laminates

---

### Step 3: Modal Analysis
```bash
/INPUT, modal_script.apdl
```
**Function:** Computes natural frequencies and mode shapes.

**Output:**
- Mode frequencies (f₁, f₂, f₃, ...)
- Mode shapes (eigenvectors)
- Damping estimates from subsequent free vibration analysis
- Results stored in `.rst` file

**Next step:** Extract frequencies and modal participation factors for further analysis

---

### Step 4: Transient Dynamic Analysis
Choose one based on load scenario:

#### For moving pressure loads:
```bash
/INPUT, script_transient.apdl
```
**or** (for improved convergence):
```bash
/INPUT, script_transient_conv.apdl
```

**Function:** Time-domain integration of the equation of motion under moving load.

**Input parameters (in script):**
- Load trajectory: x(t) = x₀ + ut, y(t) = y₀ + vt
- Load amplitude: pmax
- Impulse duration: T_imp
- Time stepping parameters

**Output:**
- Nodal displacements over time
- Ready for strain extraction at sensor locations

#### For Static analysis:
```bash
/INPUT, script_static.apdl
```

---

### Step 5: Data Extraction

#### Extract General Results
```bash
/INPUT, data_exporter.apdl
```
**Output:** CSV file with:
- Nodal displacements
- Reaction forces
- Global stress/strain components

#### Extract Strain at FBG Sensors 
```bash
/INPUT, def_exporter.apdl
```
**Output:** CSV file with strain time-series at each of 8 FBG locations
- Configured for 8 sensors (easily modifiable)
- Extracts strain components in local coordinates
- These are the **primary input for the inverse problem**

#### Extract Mode Shapes
```bash
/INPUT, mode_exporter.apdl
```
**Output:** Mode shape data for visualization and analytical model calibration

---

### Step 6: Free Vibration Analysis (Optional)
```bash
/INPUT, free_vibration.apdl
```
**Function:** Extract impulse response and damping ratio from a selected mode.

**Use for:**
- Experimental validation (compare with test data)
- Damping estimation
- Modal frequency verification

---

## 🚀 Reduced Order Model (ROM) - Inverse Problem Solver

The `reduced_order_model.py` is the **primary tool** for solving the inverse load identification problem. It's based on a **Rayleigh-Ritz analytical model** combined with **variable projection** and **multi-stage optimization**, offering significant computational advantages over pure FEM approaches.

### Key Advantages
- ⚡ **~100x faster** than full FEM transient analysis
- 🎯 **Accurate** modal properties and strain predictions
- 🔧 **Flexible** for different load scenarios and sensor configurations
- 📊 **Robust** optimization with multi-level refinement

### How It Works
1. **Builds a Rayleigh-Ritz modal model** from laminate properties
   - Classical Laminate Theory (CLT) for stiffness matrix
   - Pre-diagonalizes M, K matrices (one-time cost)
   
2. **Integrates with ODE solvers** (RK45)
   - Decoupled SDOF equations in modal coordinates
   - Much faster than coupled FEM solver
   
3. **Solves the inverse problem**
   - **Grid search** (coarse global exploration)
   - **Nelder-Mead** (local refinement)
   - **L-BFGS-B** (final polishing with bounds)
   - **Variable Projection** for separable parameters (pmax estimation)

4. **Optional calibration**
   - Scales model predictions to match experimental reference
   - Corrects for modeling uncertainties

---

### Quick Start

```bash
cd Analytical/

# Basic usage: identify load from strain measurements
python reduced_order_model.py \
    --csv strain_measurements.csv

# With calibration: scale model predictions to known load
python reduced_order_model.py \
    --csv unknown_load.csv \
    --cal-csv known_load.csv \
    --cal-pmax 15000 \
    --cal-phi -0.25 0.0 5.0 0.0
```

---

### Complete Parameter Reference

```bash
python reduced_order_model.py \
    --csv <strain_data.csv> \
    [--dt <timestep_ms>] \
    [--t-start <seconds>] \
    [--t-end <seconds>] \
    [--n-poly <order>] \
    [--grid-n <points_per_axis>] \
    [--pmax-min <Pa>] \
    [--pmax-max <Pa>] \
    [--lambda-reg <weight>] \
    [--T-imp <impulse_duration>] \
    [--cal-csv <calibration_file>] \
    [--cal-pmax <known_amplitude>] \
    [--cal-phi <xs0 ys0 u v>]
```

#### Input Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--csv` | ⭐ required | Path to strain measurements (8 columns for FBG1-8) |
| `--dt` | auto-estimated | Time step (seconds). If omitted, estimated from peak position |
| `--t-start` | 0 | Start time window (seconds) |
| `--t-end` | max | End time window (seconds) |
| `--n-poly` | 8 | Ritz polynomial order (8=coarse, 10=moderate, 12+=refined) |
| `--grid-n` | 20 | Grid search resolution per axis (higher=finer but slower) |
| `--pmax-min` | 1000 | Minimum expected load amplitude (Pa) |
| `--pmax-max` | 80000 | Maximum expected load amplitude (Pa) |
| `--lambda-reg` | 0.5 | Regularization strength for pmax bounds |
| `--T-imp` | 0.1 | Impulse duration (seconds) |
| `--no-header` | — | Skip first row of CSV |
| `--k_m` | 0 | Membrane stiffness (advanced tuning) |

#### Calibration Parameters
| Parameter | Description |
|-----------|-------------|
| `--cal-csv` | Reference CSV with known load parameters |
| `--cal-pmax` | True load amplitude from reference (Pa) |
| `--cal-phi` | True load parameters: `xs0 ys0 u v` |
| `--cal-dt` | Timestep for calibration file (if different) |

---

### Example Workflows

#### Scenario 1: Simple Load Identification
```bash
# Identify load from FEM synthetic data
python reduced_order_model.py \
    --csv ../Results/FEM_strain_motion1.csv \
    --n-poly 10 \
    --pmax-min 5000 \
    --pmax-max 50000
```

**Output:**
- Identified parameters: xs0, ys0, u, v, pmax
- RMS error per sensor
- Comparison plot: measured vs reconstructed strain
- Figure saved as `inverse_result_v2.png`

---

#### Scenario 2: With Experimental Calibration
```bash
# Calibrate model using known reference case
# Then identify unknown load
python reduced_order_model.py \
    --csv unknown_load_case.csv \
    --cal-csv known_reference.csv \
    --cal-pmax 25000 \
    --cal-phi -0.25 0.00 5.0 0.0 \
    --n-poly 11
```

**Workflow:**
1. Model runs forward on calibration case with known parameters
2. Computes scaling factor `k_cal = pmax_known / pmax_ritz`
3. Applies same `k_cal` to final identified load
4. Corrects systematic modeling bias

---

#### Scenario 3: Refined Analysis (Fine Grid, Higher Order)
```bash
# For high-accuracy identification with more time
python reduced_order_model.py \
    --csv high_res_measurement.csv \
    --n-poly 12 \
    --grid-n 30 \
    --lambda-reg 1.0
```

**Trade-off:**
- Higher `n-poly` → better accuracy but slower forward evaluations
- Higher `grid-n` → better global search but more evaluations (~30⁴ = 810,000 combos)
- `--grid-n 20` is reasonable default

---

### Available Scenarios (Load Motions)

The model automatically infers the load motion from strain data. Three standard scenarios are:

| Scenario | Position | Velocity | Pattern |
|----------|----------|----------|---------|
| **Moto 1** | xs0=-0.25, ys0=0.00 | u=5.0 m/s, v=0 | Pure X motion (centered) |
| **Moto 2** | xs0=-0.25, ys0=1.00 | u=5.0 m/s, v=0 | X motion (offset in Y) |
| **Moto 3** | xs0=-0.50, ys0=-0.50 | u=10.0 m/s, v=10.0 | Diagonal motion (symmetric) |

The algorithm searches over these parameter space automatically.

---

### Understanding the Output

```
=============================================================
  RISULTATI (Results)
=============================================================
  xs0   = -0.2500 m        ← Initial X position (found)
  ys0   =  0.0050 m        ← Initial Y position (found)
  u     =  4.9850 m/s      ← X velocity (found)
  v     =  0.0150 m/s      ← Y velocity (found)
  p_max Ritz      = 14850.25 Pa       ← Model-estimated amplitude
  k_cal           = 1.0100            ← Calibration scale factor
  p_max calibrato = 15000.00 Pa       ← Corrected amplitude
  J*    = 0.0234 / 8.0                ← Error (lower is better)
=============================================================

  Errore RMS per sensore:
    FBG1: 0.045 ue (2.3%)   ← Low error, good fit
    FBG2: 0.052 ue (2.8%)
    ...
    FBG8: 0.041 ue (1.9%)
```

**Interpretation:**
- **J* near 0** = excellent match (model explains data well)
- **RMS < 5%** = acceptable fit
- **RMS > 10%** = investigate mismatch (wrong scenario, noise, or modeling issue)

---

### Troubleshooting

**Q: "dt stimato" seems wrong?**  
A: Provide explicit `--dt` value based on your actual sampling frequency.

```bash
# If you know sampling was at 5 kHz
python reduced_order_model.py --csv data.csv --dt 0.0002
```

**Q: Optimization keeps converging to (0, 0)?**  
A: Widen search bounds or reduce regularization:
```bash
--lambda-reg 0.1  # Softer bounds on pmax
```

**Q: Runtime too slow?**  
A: Use coarser grid or lower polynomial order:
```bash
--grid-n 15 --n-poly 8    # Faster but less accurate
```

**Q: Results show k_cal far from 1.0?**  
A: Model systematically over/under-estimates. Check:
- Laminate properties (E, G, nu)
- Sensor position in Z (z_fbg parameter)
- Load model assumptions (1/R decay, impulse shape)

---

## 📊 Data Flow Diagram

```
OPTION A: Full FEM Workflow
─────────────────────────────

┌─────────────────┐
│  CAD Model      │
│ (Models/*.x_t)  │
└────────┬────────┘
         │
         ▼
┌──────────────────────────────────────┐
│  FEM Preprocessing & Setup (ANSYS)   │
│  1. geometry_loader.apdl             │
│  2. SYS_orienter.apdl                │
└────────┬─────────────────────────────┘
         │
    ┌────┴────────────┐
    ▼                 ▼
┌────────────┐    ┌─────────────────┐
│ Modal      │    │ Transient       │
│ Analysis   │    │ Dynamic FEM     │
└────────────┘    └────────┬────────┘
                           │
                ┌──────────┴──────────┐
                ▼                     ▼
         ┌──────────────┐      ┌─────────────┐
         │ def_exporter │      │ mode_exporter
         └──────┬───────┘      └──────────────┘
                │
                ▼
        ┌──────────────────┐
        │ strain_           │
        │ measurements.csv  │ ← 8 columns (FBG1-8)
        └────────┬─────────┘
                 │
                 ▼
    ┌────────────────────────────────┐
    │  reduced_order_model.py        │
    │  (Rayleigh-Ritz Inverse)       │
    └────────┬───────────────────────┘
             │
             ▼
    ┌──────────────────────────┐
    │ Identified Load Params   │
    │  xs0, ys0, u, v, pmax    │
    │  + error metrics         │
    │  + reconstruction plot   │
    └──────────────────────────┘


OPTION B: Fast Path (Experimental Data + Calibration)
──────────────────────────────────────────────────────

┌──────────────────────┐      ┌─────────────────────┐
│ Experimental Strain  │      │ FEM Reference Case  │
│ (unknown load)       │      │ (known parameters)  │
└──────────┬───────────┘      └──────────┬──────────┘
           │                             │
           │                             ▼
           │                    ┌──────────────────┐
           │                    │ strain_ref.csv + │
           │                    │ --cal-pmax 25000 │
           │                    │ --cal-phi -0.25  │
           │                    └────────┬─────────┘
           │                             │
           └─────────────┬───────────────┘
                         ▼
            ┌────────────────────────────────┐
            │  reduced_order_model.py        │
            │  (with calibration factor k)   │
            └────────┬───────────────────────┘
                     │
                     ▼
            ┌──────────────────────────┐
            │ Calibrated Load Params   │
            │  xs0, ys0, u, v, pmax*   │
            │  (pmax* = pmax/k_cal)    │
            └──────────────────────────┘
```

---

## 🔄 Complete Analysis Workflow

### For Solving an Inverse Load Problem:

**Option A: From FEM Simulation (Full Workflow)**

**1. FEM Setup & Solution (ANSYS)**
```
Open pannello4.wbpj in ANSYS Workbench
Run: geometry_loader.apdl       (load geometry, create mesh)
Run: SYS_orienter.apdl          (setup local coordinates)
Run: modal_script.apdl          (extract natural frequencies)
Run: script_transient.apdl      (run dynamic analysis)
```

**2. Extract Strain Data**
```
Run: def_exporter.apdl          (extract strain at FBG sensors) ⭐ CRITICAL
→ Produces: strain_measurements.csv (8 columns)
```

**3. Solve Inverse Problem (Python)**
```bash
cd Analytical/
python reduced_order_model.py \
    --csv strain_measurements.csv \
    --n-poly 10 \
    --grid-n 20
```
→ Outputs: identified load parameters + visualization

---

**Option B: From Experimental Measurements (Fast Path)**

If you already have experimental strain data:
```bash
cd Analytical/
python reduced_order_model.py \
    --csv experimental_strain.csv \
    --cal-csv fem_reference.csv \
    --cal-pmax 25000 \
    --cal-phi -0.25 0.0 5.0 0.0
```

No need to run full FEM if calibration reference exists!


---

## 🧪 Test Scenarios

Three pre-configured load motions are available:

### Scenario 1: Pure X-Motion
```python
xs0, ys0 = -0.25, 0.00
u, v = 5.0, 0.0
# Load enters from left, crosses center, exits right
```

### Scenario 2: Parallel Offset
```python
xs0, ys0 = -0.25, 1.00
u, v = 5.0, 0.0
# Similar to Scenario 1 but offset in Y
```

### Scenario 3: Diagonal Motion
```python
xs0, ys0 = -0.50, -0.50
u, v = 10.0, 10.0
# Load enters diagonally, moves across plate
```

Modify `--moto` parameter in Python scripts to switch scenarios.

---

## 📁 Key Files

### ANSYS Scripts (Scripts/)
| File | Lines | Purpose |
|------|-------|---------|
| `geometry_loader.apdl` | 208 | Load CAD, generate mesh |
| `SYS_orienter.apdl` | 190 | Define local coordinates for sensors |
| `modal_script.apdl` | 13 | Modal analysis |
| `script_transient.apdl` | 51 | Dynamic time-domain solver |
| `script_transient_conv.apdl` | 43 | Solver with convergence control |
| `def_exporter.apdl` | 75 | **Extract strain at FBG sensors** |
| `data_exporter.apdl` | 42 | Export general FEM data |
| `mode_exporter.apdl` | 50 | Extract eigenvectors |
| `free_vibration.apdl` | 66 | Analyze free vibration response |

### Python Scripts (Analytical/)
| File | Purpose | Status |
|------|---------|--------|
| **`reduced_order_model.py`** | **Main inverse solver using Rayleigh-Ritz ROM** | ✅ **Active** |
|  | • Builds analytical modal model from CLT | |
|  | • Grid search + Nelder-Mead + L-BFGS-B optimization | |
|  | • Variable projection for pmax estimation | |
|  | • Optional calibration with known reference cases | |
|  | • ~100x faster than FEM | |

**Removed (no longer maintained):**
- `direct_problem.py` — Replaced by ROM
- `direct_problem_w_inversion.py` — Deprecated v1
- `direct_problem_w_inversion_v2.py` — Deprecated v2
- `inverse_fem_problem_v2.py` — Deprecated FEM-based approach

### Model Files
| File | Format | Use |
|------|--------|-----|
| `Pannello_assieme.SLDASM` | Solidworks | Design modification |
| `Pannello_assieme.STEP` | STEP | CAD exchange |
| `pannello.x_t` | Parasolid | ANSYS import (preferred) |
| `pannello4.wbpj` | ANSYS Workbench | Main FEM project file |

---

## 📊 Results Organization

### Results/Transient/
- `strain_data.csv` — Time-series strain at 8 FBG sensors
- `*.png` — Visualization of transient response

### Results/Static/
- Static analysis results (if run)

### Results/Modal/
- Mode shapes and frequencies
- Modal participation factors

---

## 🛠️ Dependencies & Requirements

### ANSYS
- ANSYS Workbench 2021 R2 or later
- ANSYS Parametric Design Language (APDL)
- Composite Material library

### Python
```bash
pip install numpy scipy matplotlib
```

### Supported Platforms
- Windows 10/11
- Linux (WSL2 tested)
- macOS (requires ANSYS cross-platform license)

---

## 📈 Analysis Capabilities

✅ **Implemented**
- **Full FEM workflow:** Preprocessing, modal analysis, transient dynamics (ANSYS APDL)
- **Fast ROM solver:** Rayleigh-Ritz analytical model with multi-stage optimization
- **Load identification:** Recover position, velocity, and amplitude from strain measurements
- **Advanced optimization:** Grid search + Nelder-Mead + L-BFGS-B pipeline
- **Calibration:** Scale model predictions to match known reference cases
- **Variable separation:** Efficient handling of linear (pmax) vs. nonlinear (kinematics) parameters
- **Flexible input:** Auto-detects time step, handles multiple scenarios

🔜 **Future Extensions**
- Real experimental FBG data validation
- Non-linear composite behavior (damage, plasticity)
- Noise robustness analysis & regularization techniques
- Bayesian uncertainty quantification
- Deep learning for scenario classification

---

## 🔗 Related Work & References

**Project Based On:**
- Classical Laminate Theory (CLT)
- Kirchhoff-Love plate theory
- Finite Element Method (FEM)
- Inverse problem regularization (Tikhonov)
- Ritz method for analytical models

**Applications:**
- Structural health monitoring (SHM)
- Damage detection in composite structures
- Load monitoring in real-time systems
- Impact localization and identification

---

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@thesis{venturi2024loads,
  author = {Venturi, Gabriele},
  title = {Verso l'Identificazione di Carichi in Movimento da Misure di 
           Deformazione Distribuite su Piastre Composite},
  school = {University of Pisa},
  year = {2024},
  advisor = {Santamato, Giancarlo}
}
```

---

## 📧 Contact & Support

**Author:** Gabriele Venturi  
**Supervisor:** Assoc. Prof. Giancarlo Santamato  
**Department:** Department of Civil and Industrial Engineering, University of Pisa

---

## 📜 License

This project is provided for academic and research purposes. Please see individual files for specific licensing terms.

---

## ✅ Checklist for New Users

- [ ] Read this README completely
- [ ] Review `Documentation/Project_Documentation.pdf` for theory
- [ ] Open `pannello4.wbpj` in ANSYS Workbench
- [ ] Run preprocessing scripts (geometry_loader → SYS_orienter)
- [ ] Run modal analysis to verify model
- [ ] Run transient analysis for one scenario
- [ ] Extract strain data using def_exporter.apdl
- [ ] Run `direct_problem.py` for validation
- [ ] Experiment with inverse problem scripts
- [ ] Validate results with analytical models

---

**Version:** 2.0  
**Last Updated:** 2026-10-03  
**Status:** Refined with Reduced Order Model (ROM) solver, Ready for GitHub publication

### What's New in v2.0
- ✨ **Primary solver:** Replaced multiple FEM-based Python scripts with single, optimized `reduced_order_model.py`
- ⚡ **Performance:** ~100x speedup using Rayleigh-Ritz analytical model
- 🎯 **Robustness:** Multi-stage optimization (grid search → Nelder-Mead → L-BFGS-B)
- 📋 **Calibration:** Integrated reference-based calibration for systematic bias correction
- 📚 **Documentation:** Comprehensive parameter guide, examples, and troubleshooting
- 🗑️ **Cleanup:** Removed deprecated Python inversion scripts with incorrect results
