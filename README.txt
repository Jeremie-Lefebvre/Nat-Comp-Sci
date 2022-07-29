README: Talidou et al. NAT COMP SCI 2022

Results in the manuscript were generated using the custom C++ programs contained in this folder. 

Each code has the same structure, data structures and generates similarly organized data files. Parameter values used are defined in each program independently. 

To generate network simulations with previously learned conduction velocities, each program imports CV data for different network dimensions 
(Omega=1, 10, 100, 1000, called CVS1.txt, CVS10.txt, CVS100.txt and CVS1000.txt, respectively) as well as associated axonal tract lenghts (LENGTHS1.txt ,LENGTHS10.txt ,LENGTHS100.txt and LENGTHS1000.txt). These files must be placed in the same folder as the code is executed. 
To generate simulations for one particular dimension, the parameter Omega in the codes must be changed to a desired value, and the files imported should match that chosen value (e.g. if Omega=10, then import CVS10.txt and LENGTHS10.txt).
Using these, once can then import stabilized connectivities and reproduce the data presented in Figures 4, 5 and 6.

Alternatively, one can use the programs and explore ADM learning by using baseline conduction velocities.  For these, axonal tract lengths are imported the same way as mentionned previously, but CVs are not (all set to baseline, and changing afterwards).
These can be used to reproduce Figures 2 and 3. 

In all cases, the codes will output .txt datafiles that need to be plotted using a third party graphics software.

Summary of codes use and output: 

ADM - Figure 2.cpp : outputs the dynamics of our network model as ADM plasticity is turned on and started from baseline. Imports axonal lengths, based on the choice of network spatial dimension (i.e. Omega). ADM learning is ON.
Outputs the mean firing rate, spike train, time evolution of sample CVs, connectivity, conditional firing rates, stability metric, and delay/CV distribution before and after ADM plasticity. These quantities are defined in our manuscript.

ADM - Figure 3.cpp : ouputs the mean firing rate, its variance, firing rate correlations, as well as theoretical predictions across independent trials for different network spatial dimension(i.e. Omega). Imports axonal conduction velocity (CVS1-1000.txt and LENGTHS1-1000.txt). ADM learning is OFF. 
Output data files related to the network dynamics, activity, evolution of CVs etc are those for the last trial performed. 

ADM - Figure 4.cpp (various versions): Network mean firing rate and its variance as a function of changing stimulation parameters and for various waveforms, across independent trials for network spatial dimension Omega 10mm. Imports axonal conduction velocity (CVS1-1000.txt and LENGTHS1-1000.txt). ADM learning is OFF. 

ADM - Figure 5.cpp (various versions): Outputs the mutual information and its variance, across independent trials for network spatial dimension Omega 10mm. Imports axonal conduction velocity (CVS1-1000.txt and LENGTHS1-1000.txt). ADM learning is OFF. 

ADM - Supplementary Figure 2: Outputs the theoretical predictions, namely the mean firing rate, mean firing rate variance, as a function of conduction velocity. Omega ranges from 1 -1000mm. ADM learning is OFF. No data is imported. 

ADM - Supplementary Figure 3:  outputs the dynamics of our network model as ADM plasticity is turned on and started from baseline, while the network is set in a balanced state. Imports axonal lengths, based on the choice of network spatial dimension (i.e. Omega). ADM learning is ON.
Outputs the mean firing rate, spike train, time evolution of sample CVs, connectivity, conditional firing rates, stability metric, and delay/CV distribution before and after ADM plasticity. These quantities are defined in our manuscript

Data files (.xlsx) contains all the data needed to reproduce all the figures (1-5, Supplementary 2 and 3) from the paper. Name of file specifies which figure the data belongs to, while figures panels are specified within each file as different tabs. 
 
