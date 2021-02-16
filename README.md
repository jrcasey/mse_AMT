# Microbial Simulation Environment

The Microbial Simulation Environment (MSE) is a MATLAB toolbox for the simulation of genome-scale metabolic models (GEMs) in the natural environment. It was designed to relax the static parameterizations of cell physiology typical of standard flux balance analysis (FBA) simulations, by allowing for cellular acclimation processes (nutrient transport, photoacclimation, macromolecular allocation) which are important under environmental conditions which depart from 'optimal' laboratory growth conditions. 
  
The impact of physiological acclimation on fitness can be substantial. Below is a histogram of the relative increase in growth rates for a collection of 69 strains of a cyanobacterial genus *Prochlorococcus* along a meridional section of the North and South Atlantic Oceans:

![FitnessGains](https://github.com/jrcasey/mse_AMT/blob/master/docs/images/FitnessGains.jpg?raw=true)

### MSE Pipeline
This distribution (mse_AMT) includes a complete dataset for you to try out, since there is quite a bit of fairly specialized experimental data that are needed. The MSE pipeline is implemented in several sequential steps

1. A GEM is loaded
2. Environmental parameters are parsed. These currently include temperature, light (spectral downwelling irradiance), and nutrients (ammonia, nitrite, nitrate, phosphate). 
3. Specification of the GEM. Assigns initial physiological parameters from experimental data collected during exponential growth in 'optimal' laboratory growth conditions (i.e., nutrient replete, optimal temperature, optimal light). These data include cell size, macromolecular compositions, metabolomics, DNA and RNA base stoichiometry, maximum photosynthetic rates, transporter abundances from quantitative proteomics, growth and non-growth associated maintenance ATP requirements, and pigment concentrations.
4. Specification of constraints. Constraints include upper and lower bounds on cell size, macromolecular compositions, and pigment concentrations and pigment ratios. Upper bounds on cell size set upper bounds on nutrient transporter abundances, based on available membrane surface area.   
5. Nutrient transport optimization (OptTrans). Adjusts membrane transporter abundance and cell size to maximize growth rate. Based on the model of [Casey and Follows, 2020](https://jrcasey.github.io/assets/docs/CaseyFollows2020.pdf).
6. Photoacclimation (PhysOpt) - optimally adjusts pigment content to maximize growth rate, based on the *in situ* downwelling irradiance spectrum and the absorption spectra of each light-harvesting or photoprotective pigment, subject to constraints on their concentrations and photosystem ratios. Macromolecular allocation (also part of PhysOpt). optimally adjusts the macromolecular composition of biomass to maximize growth rate. The relative abundance of protein, lipid, carbohydrate, RNA, cell wall, pigments, and several intracellular dissolved metabolite pools are allowed to vary within some experimentally defined bounds, subject to the substrate uptake rate constraints determined in the nutrient transport acclimation step.  
7. $l_{1}$-norm FBA. Solves the FBA problem (max $c^{T}v$, s.t. $Sv=0$, $v^{lb}\leq v < v^{ub}$), stores the function value $f$, and then minimizes the euclidian of the flux vector $v$ subject to $Sv=0$, $c^{T}v = f$ and $v^{lb}\leq v < v^{ub}$.
8. Adjusts growth rate and fluxes for temperature effects. A parameterization of growth rate dependence on temperature is applied, specifying the optimal growth temperature predicted by a machine learning algorithm called [TOME](https://github.com/EngqvistLab/tome_cool) which requires only the translated genome sequence.

Additionally, there are quite a number of post-run clean-up steps, sanity checks and various analyses included to help you make sense of the output.


### Installation
Dependencies: [RAVEN 2.0](https://github.com/SysBioChalmers/RAVEN/wiki)\*, [Mosek 9](https://www.mosek.com/downloads/). 

1. Clone the repository to a local directory (using your SSH key): `git clone git@github.com:jrcasey/mse_AMT.git`
2. Add to MATLAB path: `addpath(genpath(path/to/mse_AMT/));`
3. Two options to run MSE: Locally: open `mse_AMT/toolbox/wrappers/AMT_Wrapper.m`, specify options, gridding, and filenames. Choose a station, depth, and strain ID and run. On a cluster: Specify options, gridding, and filenames in `mse_AMT/toolbox/wrappers/AMT_Wrapper.m`. Specify job loop parameters for SLURM array in `Job_loop.sh` Run SLURM array using `sh Job_loop.sh`

*\*You will need to sign up for a free academic license for Mosek. Not painless, unfortunately. If you plan to run mse_AMT on a cluster, you will need to install on the server, place the license in the Mosek directory, and make sure the solver is working by running the included test `mse_AMT/test/mosekdiag.m`*


