# MiMiC-ANT-testbench

Legged robots are valuable platforms to mimic animals. Multiple studies compare robots locomotion, navigation and other implemented behaviours, to performances of real animals (Ijspeert et al. 2014). Nevertheless, there is still a big challenge to define what are the performances of a legged robot in terms of locomotion efficiency and endurance. Various measurements have been conducted on walking animals. Typically, insect studies by mean of test benches are composed of motion capture, force measures (Dallmann et al. 2016), and charge transportation (Zill et al. 2004). Fewer studies dealt with the estimation of insect's energetic cost of locomotion (Lighton et al. 2016). On the other side, robot’s characteristics are represented mostly by their velocity, cost of transport, energy consumption, or structure robustness. For now, regarding the available data and methods, there is still no defined procedure grouping and estimating with a single experimental setup all the meaningful characteristics of a robotic leg. Such a setup would allow researchers to compare any legs with a same procedure.   

In this context, we decided to build a new robotic leg test bench, combining the measurements techniques applied on animals and robots in a one and a same test bench. Our test bench allows to establish a robotic leg’s performance table, required mainly, for robot performances assessment, structure and locomotion optimization. Moreover, we used the resulting values to compare a classical hexapod robotic leg with an ant to better understand the gap between a robot and its biological model.

https://hal.archives-ouvertes.fr/hal-03464001v1  


1. System requirements:
The complete set of procedures was tested on:
Windows 10
Minimum 16 Go RAM is recommended
Recommended free disk storage of 20 Go (depending the database size)
Intel i7 2.3 GHz 8 cores
Matlab version R2018a

2. Files list  
> System folder  
Contains the MiMiC-Ant test bench script controlling the sensors, and recording processes.  
This folder includes a custom database management script.  

> PostProcess folder
Groups all the script used to process the raw data, and plot energetic profile curves.  
Postprocess is composed of 3 main parts: full robot calculations, simulation, and test bench data processing.  

PostProcess script need additional tools from: https://github.com/IlyaBrod/vz_tools  

> pictures  
Test bench pictures and result curves.  

> data  
Some data used for the calculations.  
Full datasets used during the current study are available from the corresponding author on reasonable request.  


For more informations, see the readme file in each folder.  



