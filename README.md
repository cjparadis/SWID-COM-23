SWID-COM-23 is a computer program that numerically (finite difference) simulates solute transport (advection, mechanical dispersion, linear-equilibrium sorption, and first-order decay) during a single-well injection-drift (SWID) test and is fully executable from the command line via user inputs and automated outputs.

Below is a list and a brief description of the ten files that come with the program. Note that both the computation (GNU Fortran) and visualization (GNUPLOT) software are open source, i.e., SWID-COM-23 and the software used to build it are free.

FILES:

#1
Name: read_me.text
Description: readme file that you are reading right now :)

#2
Name: swid_test_v4.f95
Description: source code written in Fortran 95

#3
Name: swid_test_v4.exe	
Description: compiled source code via GNU Fortran (GCC) 12.2.0

#4
Name: swid_out_s1.conc
Description: default output file for concentration versus radial distance about the well at the end of stress period 1 (injection phase)

#5
Name: swid_out_s2.conc	
Description: default output file for concentration versus time at the well during stress period 2 (drift phase)

#6
Name: swid_out_s1.disc	
Description: default output file for mass budget for during stress period 1 (injection phase)

#7
Name: swid_out_s2.disc	
Description: default output file for mass budget for during stress period 2 (drift phase)

#8
Name: swid_plot_conc.gnu	
Description: source code written for GNUPLOT 5.2.8 for visualization of the default output files for concentration

#9
Name: swid_plot_disc.gnu	
Description: source code written for GNUPLOT 5.2.8 for visualization of the default output files for mass budget

#10
Name: swid_comp.txt
Description: detailed information on the computational methods (numerical finite-difference model) of the program

HOW TO RUN THE PROGRAM (WINDOWS):

1) Download the folder named: SWID-COM-23
2) Open the Windows command prompt
3) Change the directory to the drive and path of the downloaded folder 
4) Execute the program by typing "swid_test_v4.exe" in the command prompt followed by hitting enter
5) Enter the input parameters for the injection phase, default parameters are used if you simply enter a single comma
6) Enter the input parameters for the drift phase, default parameters are used if you simply enter a single comma for each of the three simulations
7) Enter "Yes" or "No" for simulations of the injection phase from an analytical model, default is "Yes" via typing a single comma
8) Enter "Yes" or "No" for a mass budget calculation, default is "Yes" via typing a single comma
9) Examine both the command line output and the output files (#4 through #7 above) for the results of the simulations
10) Optional: run files #8 and #9 via GNUPLOT to graphaically visualize the simulations

QUESTIONS, COMMENTS, & BUG REPORTS:

Write to Dr. C. J. Paradis at paradisc@uwm.edu with a subject heading as follows: SWID-COM-23
