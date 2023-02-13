# Difference-Density-Natural-Orbitals
This program carries out the Density Difference Natural Orbital analysis
Computed the promotion number, excitation number, and the DDNO's for
an excited state calculation

-A. J. Bovill, 2022

The Promotion number was proposed by Martin H. Gordon et al. to find the
number of electrons that have been promoted to an excited state.
"Analysis of Electronc Transitons as the Difference of Electron Attachment
and Detachment Densities"
DOI: 10.1021/j100039a012

The Excitation number was an improvement over the promotion number and was 
proposed by Peter M. W. Gill et al. to represent the
number of electrons in the excited state that lie 
in the unoccupied space of the ground state.
"Excitation Number: Characterizing Mutiply Excited States"
DOI: 10.1021/acs.jctc.7B00963

The Density Difference Natural Orbital anaylsis in this program is an
extension of previous work in our lab involving Natural Ionization
Orbitals (NIO') and applied to vertical excitations instead of
ionizations.
"Natural ionization orbitals for interpreting electron detachment
processes"
DOI: 10.1063/1.5941738
Previous NIO python code can be found here:
https://github.com/hphratchian/nio 

This program uses our lab' Merced Quantum chemistry program for working
with matrix files produced from gaussian calculations. The OOP format
hides away mathematical and works with a number of quantum chemistry codes
for easy mathematical interpretation
MQC code can be found here:
https://github.com/MQCPack/mqcPack
