# 2D MASWavesPy

## Overview

2D MASWavesPy extends the functionality of the original **MASWavesPy** Python package developed by Magnús Snorri Bjarnason (Faculty of Civil and Environmental Engineering, University of Iceland).

The original package consists of the modules:

- `wavefield`
- `dispersion`
- `combination`
- `inversion`

with helper modules:

- `dataset`
- `select_dc`

The `wavefield` module provides methods to import recorded shot gathers as `RecordMC` objects. The phase shift method (1) is used to transform each shot gather into the frequency-phase velocity domain. The `dataset` module can be used to import a set of shot gathers in the form of a `Dataset` object through a .csv file. 

The `dispersion` module, along with the supplementary `select_dc` module, provides methods for visualization of the phase velocity spectrum and dispersion curve (DC) identification using a GUI (Graphical User Interface). An `ElementDC` object stores the frequency-phase velocity domain representation of a given `RecordMC` and the corresponding DC (referred to as an elementary DC). 

The `combination` module provides methods to combine elementary DCs obtained from multiple shot gathers into a composite DC (2) (a `CombineDCs` object) and to assess and view the spread in the dispersion data, either as a function of frequency or wavelength. A `Dataset` object can contain multiple pairs of `RecordMC` and `ElementDC` objects (one pair for each shot gather) and provides routines for initializing a `CombineDCs` for the set of records or a particular subset of records. 

The `inversion` module provides methods to evaluate the shear wave velocity profile of the tested site. The inversion methods, along with routines for post-processing of the inversion results, are defined on an `InvertDC` object that is initialized using an experimental DC. The fast delta matrix algorithm (3) is used for forward computations and a Monte-Carlo global search algorithm (4) for searching the solution space for the optimal set of model parameters. 

A more comprehensive description is provided in (5). 

## 2D Implementation

This 2D implementation adds support for the typical CMP (Common Midpoint) acquisition geometry used in standard 2D MASW surveys and provides:

- CMP sorting
- Dispersion curve calculation
- Interactive dispersion curve picking
- Bayesian inversion
- 2D shear-wave velocity section generation

Processing parameters are controlled through a YAML configuration file.

---

## Installation
The requirements for the packages are:

- obspy (1.5.0)
- cython (3.2.9)
- pandas (3.0.5)
- pyyaml (6.0.3)

It is best to create an environment for the code e.g.

```bash
mamba create -n 2Dmaswavespy python=3.11
```

Activate the environment
```bash
mamba activate 2Dmaswavespy
```

The easiest installation is using the provided environment.yml file.

```bash
mamba env create -f environment.yml
```

or if you don't use mamba

```bash
conda env create -f environment.yml
```

## Manual Installation

If you want to install packages manually you can use

```bash
pip install obspy
```

```bash
pip install cython
```

```bash
pip install pandas
```

```bash
pip install pyyaml
```

and compile the cython components of maswavespy

```bash
pip install -e .
```

which uses setup.py. Note that in the process seom packages might be
reinstalled.

You might also want to install ipython (if using notebooks)

```bash
pip install ipython
```

You might need to add the path to the maswavespy repository to your
Python path

```bash
export PYTHONPATH=("$PYTHONPATH:/path/to/your/2Dmaseavespy")
```

or similar for your shell. 


---
# Running the code

## Required Input Files

### 1. Configuration File

Example:

```text
MASWavesPy_config.yaml
```

The configuration file should define:

- Site name
- Profile name
- SEG2 inventory file
- Initial velocity model
- Processing options

---

### 2. SEG2 Inventory

A text file containing paths to the original shot gather files, one file per line.

Example:

```text
Example_Data/3006.dat
Example_Data/3007.dat
Example_Data/3008.dat
```

Example filename:

```text
Line_3000.csv
```

---

### 3. Initial Velocity Model

Example:

```text
ElCuchillo_initial.csv
```

Required format:

| Layer | Thickness h (m) | Vs (m/s) | Density (kg/m³) | Saturation | Vp (m/s) | ν |
|---------|---------|---------|---------|---------|---------|---------|
| 1 | 1 | 150 | 1850 | sat | 400 | 0.3 |
| 2 | 2 | 175 | 1900 | sat | 550 | 0.3 |
| 3 | 8 | 200 | 1950 | sat | 900 | 0.3 |
| 4 (half-space) | - | 750 | 1950 | sat | 1500 | 0.3 |

---

## Workflow Control

Processing stages can be enabled or disabled individually.

Example:

```yaml
cmp_sorting:
  enabled: true

processing:
  enabled: true

plotting:
  create_section: true
```

---

## Running the Code

```bash
python MASWavesPy_processing.py MASWavespy_config.yaml
```

---

# CMP Sorting

CMP spacing is defined as a multiple of the geophone spacing.

Example:

```yaml
cmp_bin_factor: 2
```

This creates CMP bins that are twice the geophone spacing.

Output directory:

```yaml
output_dir: CMP_Gathers
```

For standard CMP processing:

```yaml
forward_only: true
```

This uses only forward shots and maintains consistency with SeisImager2D.

### Notes

- Geometry files are generated automatically for each CMP.
- These geometry files are later used when reading the dataset in CMP mode.
- If CMP sorting is skipped, a valid `CMP_inventory.csv` must already exist.

---

# Dispersion Curve Picking

Each CMP must be processed individually.

## Step 1: Review the Wavefield

The CMP dataset is displayed as a record section with increasing source-receiver distance.

Information such as:

- Source location
- Receiver location
- Offset range

is printed to the command line.

After reviewing data quality, close the wavefield window.

---

## Step 2: Compute the Dispersion Image

The software calculates the CMP dispersion image and opens the interactive picking GUI.

### Picking Options

You can define phase velocities using:

#### Method A

Select a frequency range interactively using a mouse-drawn box.

#### Method B

Pick individual phase velocity measurements manually.

#### Method C

Provide frequency indices manually using a comma-separated list:

```text
10-15,16,18,...
```

Selected points become labelled and highlighted.

When finished:

1. Click **Stop**
2. Click **Save Dispersion Curve**
3. Click **Close**

---

## Step 3: Initial Model Comparison

The picked dispersion curve is plotted together with the theoretical dispersion curve generated from the initial velocity model.

The initial misfit is printed to the command line.

A poor initial fit is acceptable provided inversion search bounds are sufficiently wide.

---

# Bayesian Inversion

Closing the comparison window starts the inversion.

Progress is displayed in the terminal.

Several diagnostic windows are shown during processing.

## Sampled Models

Displays:

- All sampled models
- Associated misfits
- Corresponding Vs structures

Higher-quality models are shown using brighter colours.

Review the parameter space coverage and close the window.

---

## Accepted Vs Profiles

Displays:

- All accepted velocity profiles
- Rejected models in grey
- Potential multimodal behaviour in the posterior distribution

Review and close the window.

---

## Median Velocity Model

Displays:

- Median Vs profile
- Velocity uncertainty
- Layer-depth uncertainty
- Theoretical dispersion curve of the median model

Review uncertainty estimates and close the figure.

---

## Best-Fitting Models

Displays:

- Ten lowest-misfit Vs profiles
- Their corresponding dispersion curves

Additional statistical summaries are printed to the terminal, including:

- Median velocities at selected depths
- Lowest-misfit model properties

After closing, processing automatically continues to the next CMP.

---

# Final 2D Section

After all CMPs have been processed, the final 2D velocity model is created.

The output contains:

## Velocity Models

- Blocky inversion model
- Gaussian-smoothed velocity model

## Uncertainty Panels

- Velocity uncertainty
- Interface-depth uncertainty

---

## Saved Output

All inversion products are saved to the `Results` directory:

- Picked dispersion curves
- Inverted models
- Pickle files

The final 2D section is additionally exported as:

```text
Results/*.pdf
```

Intermediate figures are***not automatically saved** and sh*uld be saved manually if required.*
---

# Acknowledgements

This pro*ect builds upon **MASWavesPy**:

h*tps://github.com/Mazvel/maswavespy*
## References

1. *Park, C.B., Miller, R.D., Xia, J.* (1998).  
   *Imging dispersion curves of surface *aves on multi*channel record.*  
   SEG Technica* Program Expanded Abstracts.  
   *ttps://doi.org/10.1190/1.1820161

2. Olafsdottir, E.A., Bessason, B.,*Erlingsson, S. (2018).  
   *Combination of dispersion curves from MA*W measurements*.  
   Soil Dynamics and Earthquake Engineering, 113, *73-487.  
   https://doi.org/10.10*6/j.soildyn.2018.05.025

3. Buchen P.W., Ben-Hador, R. (1996).  
   *Free-mode surface-wave computation*.*  
   Geophysical Journal Intern*tional, 124(3), 869-887.  
   http*://doi.org/10.1111/j.1365-246X.199*.tb05642.x

4. Olafsdottir, E.A., Erlingsson, S., Bessason, B. (2020)*  
   *Open-Source MASW Inversion *ool Aimed at Shear Wave Velocity Profiling for Soil Site Explorations**  
   Geosciences, 10(8), 322.  
*  https://doi.org/10.3390/geoscien*es10080322

5. Olafsdottir, E.A., *essason, B., Erlingsson, S., Kayni*, A.M. (2024).  
   *A Tool for Pr*cessing and Inversion of MASW Data and a Study of Inter-Session Variability of MASW.*  
   Geotechnical Testing Journal, 47(5), 1006-1025. *
   https://doi.org/10.1520/GTJ202*0380
