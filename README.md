# Simulating the common factor model
This repository contains R code that can be used to simulate the
common factor model. It is based on an approach described by Tucker,
Koopman, and Linn [[1]](#1). This approach was subsequently used by
MacCallum, Widaman, Zhang, and Hong [[2]](#2) to identify study design
features that determined the minimum *N* needed to obtain stable
factors.

The code was based on/translated from SAS code that used PROC IML to
simulate factors using Tucker et al.'s approach.[[3]](#3). To use it
set the following variables:

nvars: Number of items/variables

nfactors: Number of factors

d_frac: Proportion of variables that are dichotomous

commun_type: 1 = low, 2 = wide, 3 = high (see MacCallum et al. for
definitions)

replicat: Number of samples to draw from the population

numobs: Number of observations (e.g., participants) per sample

## To do list
1. Turn into an R function of some sort
2. Allow users to input their factor loading pattern and communalities

## References
<a id="1">[1]</a>
Tucker L, Koopman RF, Linn RL. 1969. Evaluation of factor analytic research procedures by means of 
simulated correlation matrices. Psychometrika 34:421-459. http://dx.doi.org/10.1007/BF02290601

<a id="2">[2]</a>
MacCallum RC, Widaman KF, Zhang S, Hong S. 1999. Sample size in factor
analysis. Psychol Meth
1:84-99. https://psycnet.apa.org/doi/10.1037/1082-989X.4.1.84

<a id="3">[3]</a>
Coughlin K, Kromrey J, Hibbard S. 2013 Using predetermined factor
structures to simulate a variety of data conditions. Southeast SAS
	Users Group.
	https://www.researchgate.net/publication/282132032_Using_Predetermined_Factor_Structures_to_Simulate_a_Variety_of_Data_Conditions
