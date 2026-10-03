# Load Identification in Composite Plates Using Distributed Strain Measurements

**Towards the Identification of Moving Loads from Deformation Measures Distributed on Composite Plates**

**Author:** Gabriele Venturi  
**Supervisor:** Assistant Professor Giancarlo Santamato  
**Institution:** University of Pisa

## 📋 Project Overview

This project addresses the **identification of dynamic moving loads** on composite laminated structures using distributed Fiber Bragg Grating (FBG) sensors. The work combines:

- **Finite Element modeling** of composite plates
- **Analytical models** for load propagation
- **Inverse problem solving** to recover load parameters from sensor measurements
- **Data-driven approaches** for system identification

### Key Features
- 🔧 **8 FBG sensors** embedded in composite laminate (oriented at 0°, 90°, ±45°)
- 📊 **Moving load model** based on pressure wave propagation
- 🎯 **Forward-Inverse framework** for load identification
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
CAD_FEM/
├── Models/                      # CAD geometries (SLDASM, STEP, x_t, IGS)
├── Documentation/               # Project documentation (PDF)
├── Results/                     # Extracted results and visualizations
│   ├── Transient/              # Time-domain analysis results
│   ├── Static/                 # Static analysis results
│   └── Modal/                  # Modal analysis results
├── Scripts/                     # ANSYS APDL scripts (preprocessing, solving, post-processing)
├── Simulations/
│   ├── Sim_1, Sim_2, Static, Static2/  # FEM workbench setups
│   └── Analytical/             # Python scripts for inverse problem
├── Archive/                     # Archived/old versions
└── README.md                    # This file
```

---

## 🔧 ANSYS Workflow (Scripts/)

### Step 1: Geometry & Mesh Setup
```bash
# Open ANSYS Workbench and load the project
pannello4.wbpj

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
**Function:** Creates local coordinate systems on each FEM element aligned with FBG sensor axes.

**Why it matters:** 
- FBG sensors measure strain along their fiber axis
- Local coordinates allow proper strain extraction at sensor locations
- Supports multi-directional sensors (0°, 90°, ±45°)

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

#### For quasi-static analysis:
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

#### Extract Strain at FBG Sensors ⭐ **CRITICAL**
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

## 🐍 Inverse Problem & Analysis (Simulations/Analytical/)

After extracting strain data from FEM, use Python scripts to solve the inverse problem:

### Quick Start
```bash
cd Simulations/Analytical/

# 1. Forward model validation
python direct_problem.py \
    --csv ../../Results/data_extracted.csv \
    --moto 1 \
    --k_m 2400 \
    --n-poly 11
```

### Available Scenarios

| Scenario | Position (m) | Velocity (m/s) | Description |
|----------|--------------|----------------|-------------|
| Moto 1   | (-0.25, 0.00) | (5.0, 0.0)    | Load moving in X direction |
| Moto 2   | (-0.25, 1.00) | (5.0, 0.0)    | Load moving parallel but offset |
| Moto 3   | (-0.50, -0.50) | (10.0, 10.0) | Diagonal motion |

---

### Phase 1: Forward Model
```bash
python direct_problem.py \
    --csv <FEM_data.csv> \
    --moto 1 \
    --k_m 2400 \
    --n-poly 11 \
    --grid-n 5
```

**Function:** 
- Loads FEM data extracted from ANSYS
- Builds analytical 1D Ritz model of the plate
- Simulates structural response using modal superposition
- Predicts strain at FBG sensors
- Compares with reference analytical solution

**Parameters:**
- `--csv`: Path to exported FEM data (from `data_exporter.apdl`)
- `--moto`: Load scenario (1, 2, or 3)
- `--k_m`: Membrane stiffness (adjust for material properties)
- `--n-poly`: Polynomial order for Ritz basis (11 recommended)
- `--lambda-reg`: Regularization parameter for inverse (default 0.5)

**Output:**
- Plots comparing FEM vs analytical response
- Error metrics
- Extracted frequencies and mode contributions

---

### Phase 2: Inverse Problem with Calibration
```bash
python direct_problem_w_inversion_v2.py \
    --csv <FEM_data_unknown_load.csv> \
    --moto 2 \
    --cal-csv <FEM_data_known_load.csv> \
    --cal-pmax 15000 \
    --cal-phi -0.25 1.0 5.0 0.0 \
    --lambda-reg 0.5
```

**Function:**
- Calibrates model on known load case
- Identifies unknown load from sensor measurements
- Recovers: position (x₀, y₀), velocity (u, v), amplitude (pmax)

**Parameters:**
- `--cal-csv`: Reference measurement for calibration
- `--cal-pmax`: Known load amplitude (Pa)
- `--cal-phi`: Known load parameters [xs0 ys0 u v]
- `--lambda-reg`: Tikhonov regularization strength

**Output:**
- Identified load parameters
- Residual errors
- Convergence plots

---

### Phase 3: Full Inverse Problem (FEM-Based)
```bash
python inverse_fem_problem_v2.py \
    --csv <strain_measurements.csv> \
    --n-poly 11 \
    --pmax-min 5000 \
    --pmax-max 80000
```

**Function:**
- Direct inverse problem without forward model
- Uses strain measurements to identify load
- More robust for noisy data

**Parameters:**
- `--pmax-min/max`: Search bounds for load amplitude

---

### Analytical Validation
Compare FEM results with 1D analytical solutions:

```bash
# Concentrated load
python 1D_concentrated_load.py

# Distributed load
python 1D_distributed_load.py

# Fourier harmonic load
python 1D_Fourier_load.py
```

These generate reference solutions for validation.

---

## 📊 Data Flow Diagram

```
┌─────────────────┐
│  CAD Model      │
│ (Models/*.x_t)  │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  FEM Preprocessing (ANSYS)          │
│  1. geometry_loader.apdl            │
│  2. SYS_orienter.apdl               │
└────────┬────────────────────────────┘
         │
    ┌────┴────┐
    ▼         ▼
┌────────┐  ┌─────────────────┐
│ Modal  │  │ Transient FEM   │
│Analysis│  │ (Dynamic Loads) │
└────┬───┘  └────────┬────────┘
     │               │
     └───────┬───────┘
             ▼
     ┌──────────────────┐
     │ Data Extraction  │
     │ def_exporter     │ ← STRAIN DATA
     └────────┬─────────┘
              │
              ▼
    ┌─────────────────────────────────┐
    │  Python Inverse Problem         │
    │  1. direct_problem.py           │
    │  2. direct_problem_w_inv_v2.py  │
    │  3. inverse_fem_problem_v2.py   │
    └────────┬────────────────────────┘
             │
             ▼
    ┌──────────────────────┐
    │ Identified Load      │
    │ Parameters: (x₀, y₀) │
    │           (u, v)     │
    │           (pmax)     │
    └──────────────────────┘
```

---

## 🔄 Complete Analysis Workflow

### For a New Load Case:

**1. FEM Setup (ANSYS)**
```
Open pannello4.wbpj in ANSYS Workbench
Run: geometry_loader.apdl
Run: SYS_orienter.apdl
```

**2. Solve Dynamics**
```
Run: modal_script.apdl          (to characterize system)
Run: script_transient.apdl      (with your load scenario)
```

**3. Extract Results**
```
Run: data_exporter.apdl         (general data)
Run: def_exporter.apdl          (strain at sensors) ⭐ IMPORTANT
Run: mode_exporter.apdl         (mode shapes)
```

**4. Reduced Order Model**
```bash
cd Analytical/
python reduced_order_model.py \
    --moto 1 \
    --lambda-reg 0.5
```


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

### Python Scripts (Simulations/Analytical/)
| File | Lines | Purpose |
|------|-------|---------|
| `direct_problem.py` | 317 | Forward model validation |
| `direct_problem_w_inversion.py` | 359 | Forward + inverse with calibration |
| `direct_problem_w_inversion_v2.py` | 528 | Robust version v2 |
| `inverse_fem_problem_v2.py` | 518 | Full FEM-based inverse problem |
| `1D_concentrated_load.py` | 101 | Analytical solution (point load) |
| `1D_distributed_load.py` | 156 | Analytical solution (distributed load) |
| `condition_number.py` | 105 | Numerical stability check |

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
- Forward FEM modeling with moving loads
- Modal analysis and frequency extraction
- Transient dynamic response
- Strain field computation at arbitrary points
- Inverse problem solving (load identification)
- Data-driven model calibration
- Analytical validation

🔜 **Future Extensions**
- Experimental validation with real FBG data
- Non-linear material models
- Composite damage modeling
- Uncertainty quantification
- Machine learning-based load classification

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

**Version:** 1.0  
**Last Updated:** 2024-10-03  
**Status:** Ready for GitHub publication
