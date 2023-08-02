'Dynamic General Equilibrium Modeling: Computational Methods and Applications'
(3nd edition, scheduled for fall 2022)

By Burkhard Heer and and Alfred Maußner --- Python Code for Chapter 10
____________________________________________________________________________

Comments & Questions: Burkhard.Heer@wiwi.uni-augsburg.de


In this rar-file, you find the solution to 

1. a multi-dimensional value function iteration problem with individual uncertainty (Chapter 10.1)
2. computation of the dynamics of a large-scale OLG model with the Krusell-Smith algorithm. (Chapter 10.2.2)

The material (1) is completely new and intended for a new third edition of the book. Its purpose is to solve
a multi-dimensional heterogeneous agent overlapping-generations model and the endogenous distribution of assets and cumulated contributions
to the pension system.


Chapter 10.1 (revised):
=======================

1. Model with individual uncertainty 

as described by

slides_book_heer_maussner_Chapter_10_1.pdf (attached) and preliminary
draft of Chapter 10.1.1-10.1.3 (attached)

Burkhard Heer also provides a youtube presentation:

https://youtube/D2tl5imrrNk

2. Model with contributions-based pensions as described in the preliminary
draft of Chapter 10.1.4 (attached)

Comments welcome: Burkhard.Heer@wiwi.uni-augsburg.de

You need to store the PYTHON code 'AK70_stochastic_income.py' 
and 'AK70_prog_pensions.py' and the input file with the age-dependent
survival probabilities and age-efficiency profile
'survival_probs.xlsx' in the same directory.

Chapter 10.2.2 (revised):
=========================

1. Krusell-Smith Algorithm: 

Model with individual and aggregate uncertainty as described by

slides_book_heer_maussner_Chapter_10_2_2.pdf (attached) and preliminary
draft of Chapter 10.2. (attached)


You need to store the PYTHON codes 'OLG_krusell_smith_ss.py' 
and 'OLG_krusell_smith_dyn.py' and the input file with the age-dependent
survival probabilities and age-efficiency profile
'survival_probs.xlsx' in the same directory.

First run 'OLG_krusell_smith_ss.py' to compute the non-stochastic steady state 
(input into OLG_krusell_smith_dyn). Runtime: approx 10 minutes

Second run 'OLG_krusell_smith_dyn.py' to compute the dynamics. 
RUNTIME: approx 2 DAYS

=> you may want to speed up the code using multithreading, Numba or parallelization



The files in Python_Chapter10_heer_maussner_dsge.rar are:
========================================================

* readme_ch10_python.txt: this file

* AK70_stochastic_income.py: main PYTHON program file to compute model 10.1.1

* AK70_prog_pensions.py: main PYTHON program file to compute model 10.1.4

* survival_probs.xlsx: Excel input file for the PYTHON program

* OLG_krusell_smith_ss.py

* OLG_krusell_smith_dyn.py


* slides_book_heer_maussner_Chapter_10_1.pdf: model 10.1.1 description

* slides_book_heer_maussner_Chapter_10_2_2.pdf: model 10.2.2 description


* Ch101.pdf: Preliminary draft of Chapter 10.1 for 3rd edition Heer/Maussner

* Ch102.pdf: Preliminary draft of Chapter 10.1 for 3rd edition Heer/Maussner


Instructions:
=============

In order to run the Python programs, you simply need to store the file and the Excel file in the same directory and run it. You also 
need to install libraries numpy and scipy. You may also need to install quantecon and 

Run time: approximately 52 hours (model 10.1.1) and 6 days (model with endogenous pensions, Chapter 10.1.4) 
and 2 days (model Chapter 10.2.2)

We have documented the PYTHON code AK70_stochastic_income.py in detail (line by line discussion) at the following online manuscript:

https://assets.uni-augsburg.de/media/filer_public/12/1f/121f73c8-7007-47f6-b930-7493ff3180d1/script_dge_python_chapter1013may2021.html

A good introduction into Python and a description of the installation of Python can be found on
the website of Thomas Sargent and John Starchuski

https://python-programming.quantecon.org/intro.html

In case of questions, please send me an email: Burkhard.Heer@wiwi.uni-augsburg.de

last update: March 7, 2022