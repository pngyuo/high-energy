# Jet-Induced Identified Hadron Productions in pp Collisions

This repository contains the analysis code and data for the study "Investigating jet induced identified hadron productions from the relative transverse activity classifier in pp collision at LHC" (arXiv:2507.23306).

## 📄 Paper Information
- **Title**: Investigating jet induced identified hadron productions from the relative transverse activity classifier in pp collision at LHC
- **Authors**: Yuhao Peng, Yan Wu, Xinye Peng, Zhongbao Yin, and Liang Zheng
- **arXiv**: 2507.23306
- **Journal**: Preprint
- **Collision System**: pp collisions at $\sqrt{s}=13$ TeV
- **Model**: AMPT with PYTHIA8 initial conditions

## 🎯 Research Overview
This study systematically investigates jet-associated identified hadron (π, K, p) productions using the relative transverse activity classifier $R_T$ in pp collisions. The analysis focuses on:

- Transverse momentum ($p_T$) spectra in toward and transverse regions
- Particle ratios (K/π, p/π) as functions of $R_T$ and $p_T$
- In-jet hadron production through toward-transverse subtraction
- Effects of partonic and hadronic final-state interactions
- Average transverse momentum $\langle p_T\rangle$ analysis

## 📁 Repository Structure

### 📊 Data Analysis Folders

#### `Pt/`
- Contains analysis of transverse momentum ($p_T$) spectra
- Includes toward and transverse region spectra for pions, kaons, and protons
- $R_T$-dependent spectra analysis
- Comparison with ALICE experimental data

#### `ration_Pt/`
- Particle ratio analysis as function of $p_T$
- K/π and p/π ratios in toward and transverse regions
- $R_T$ dependence studies
- Final-state interaction effects

#### `averge_Pt_fit/`
- Average transverse momentum $\langle p_T\rangle$ calculations
- Analysis across different topological regions
- $R_T$ dependence of $\langle p_T\rangle$
- Fitting procedures and results

#### `injet_dndy/`
- In-jet particle yield calculations
- Yield subtraction: $N^{\text{In-Jet}} = N^{\text{Toward}} - N^{\text{Transverse}}$
- $p_T$-differential in-jet yields
- Event activity dependence

#### `injet_ration_Pt/`
- In-jet particle ratios (K/π, p/π)
- $p_T$ dependence in different $R_T$ classes
- Crossing behavior analysis at intermediate $p_T$
- Final-state interaction effects on in-jet ratios

#### `paper_complete/`
- Complete paper manuscript and figures
- Reference management
- Supplementary materials

## 🔧 Methodology

### Event Topology Classification
- **Toward region**: $|\Delta\varphi| < 60^\circ$
- **Transverse region**: $60^\circ \leq |\Delta\varphi| < 120^\circ$
- **Away region**: $|\Delta\varphi| \geq 120^\circ$

### Relative Transverse Activity
$$R_T = \frac{N_T}{\langle N_T \rangle}$$

### In-Jet Production Definition
$$\frac{d^2N^{\text{In-Jet}}}{dp_Tdy} = \frac{d^2N^{\text{Toward}}}{dp_Tdy} - \frac{d^2N^{\text{Transverse}}}{dp_Tdy}$$

### AMPT Model Configurations
- **0 mb w/o ART**: No final-state interactions
- **0.15 mb w/o ART**: Partonic interactions only ($\sigma=0.15$ mb)
- **0.15 mb w/ ART**: Full partonic + hadronic interactions

## 🛠️ Dependencies and Requirements

### Required Software
- ROOT data analysis framework
- AMPT model with PYTHIA8 initial conditions
- Python 3.x with scientific computing libraries
- LaTeX (for paper compilation)

### Key Libraries
- NumPy, SciPy, Matplotlib
- PyROOT for ROOT integration
- Statistical analysis tools

## 📈 Key Findings

1. **$R_T$ Dependence**: Particle ratios show strong $R_T$ dependence driven by final-state interactions
2. **Crossing Behavior**: In-jet p/π ratios exhibit crossing at intermediate $p_T$ between low and high $R_T$ events
3. **Hadronic Effects**: Pion wind effect significantly enhances kaon and proton $\langle p_T\rangle$
4. **Jet Modification**: Evidence of jet medium interactions in high multiplicity pp events

## 🚀 Usage

### Running Analysis
1. Configure AMPT model parameters in respective scripts
2. Execute analysis scripts in each folder for specific observables
3. Generate plots and comparison with experimental data
4. Compile results for paper figures

### Data Processing
- Event selection based on $R_T$ classification
- Topological region separation
- Particle identification (π, K, p)
- Statistical error propagation

## 📋 Output Files

Each analysis folder generates:
- Data files with analysis results
- Plot files (PDF, PNG)
- Statistical analysis outputs
- Comparison tables with experimental data

## 📚 References

- ALICE Collaboration, JHEP 06 (2023) 027
- AMPT model: Lin et al., Phys. Rev. C 72 (2005) 064901
- PYTHIA8: Sjöstrand et al., Comput. Phys. Commun. 191 (2015) 159

## 📄 License

This project contains research code and data. Please cite the original paper when using these results.