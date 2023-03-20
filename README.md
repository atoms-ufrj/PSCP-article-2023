# 

This repository contains additional information about the article "Revisiting the Pseudo-supercritical path Method: An Improved Formulation for the Alchemical Calculation of Solid-Liquid Coexistence". 

The melting point calculation of the propenal system is considered for exemplification. It is attached free energy calculation results and some scripts applied in molecular dynamics simulations and data post-processing.

Molecular dynamics simulations were performed with LAMMPS package (version 2021). A customization of the "change_box" command was developed for the affine transformations of Stage 2.

The multistate Bennett acceptance ratio (MBAR) method was applied in the data post-processing with the PYMBAR v. 3.0.5 package, which can be friendy manipulated in our analyses with the mics package - https://github.com/craabreu/mics.
