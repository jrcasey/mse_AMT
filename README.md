# Microbial Simulation Environment

The Microbial Simulation Environment (MSE) is a MATLAB toolbox for the simulation of genome-scale metabolic models (GEMs) in the natural environment. It was designed to relax the static parameterizations of cell physiology typical of standard flux balance analysis (FBA) simulations, by allowing for cellular acclimation processes. Acclimation steps currently include:

1. Nutrient transport - optimally adjusts membrane transporter abundance and cell size to maximize growth rate. Based on the model of [Casey and Follows, 2020](https://jrcasey.github.io/assets/docs/CaseyFollows2020.pdf).
2. Photoacclimation - optimally adjusts pigment content to maximize growth rate, based on the *in situ* downwelling irradiance spectrum and the absorbtion spectra of each light-harvesting or photoprotective pigment, subject to boundary constraints on their biomass-specific abundances and ratios. 
3. Macromolecular acclimation - optimally adjusts the macromolecular composition of biomass to maximize growth rate. The relative abundance of protein, lipid, carbohydrate, RNA, cell wall, pigments, and several intracellular dissolved metabolite pools are allowed to vary within some experimentally defined bounds, subject to the substrate uptake rate constraints determined in the nutrient transport acclimation step. 

Step 1. is implemented first, and is followed by steps 2. and 3. which are solved in a single optimization. 