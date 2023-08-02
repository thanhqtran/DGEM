'Dynamic General Equilibrium Modeling: Computational Methods and Applications'
(2nd edition, 2009)

By Burkhard Heer and and Alfred Maußner --- JULIA Code for Chapter 9



The files in Julia_code_Chapter9_heer_maussner_dsge.rar are:
========================================================

readme_ch9_julia_code.txt: this file

AK60_value_main.jl: main program file, equivalent to the GAUSS program 'RCh91v.g', as described in Chapter 9.1.2
AK60_value_procs.jl: procedures for main program file "AK_value_main.jl"
=> computes the steady state in the Auerbach-Kotlikoff-60-Period OLG model in Section 9.1.2 in Heer/Maussner, 2009, 2nd ed.

Instructions:
=============

In order to run the Julia programs, you simply need to store the files in the same directory and run "AK60_value_main.jl".

In order to run the program, you should have installed JULIA and the packages (as given in the first lines of AK60_value_main.jl)

A good introduction into JULIA and a description of the installation of JULIA can be found on
the website of Jesse Perla, Thomas Sargent and John Starchuski

https://julia.quantecon.org/

If you have any questions with respect to the code, send an email to: Burkhard.Heer@wiwi.uni-augsburg.de

Burkhard Heer

last update: June 4, 2021