# Prediction of Cross-Fitness for Adaptive Evolution to Different Environmental Conditions: Consequence of Phenotypic Dimensional Reduction

Source Code of T.U.Sato, C.Furusawa and K.Kaneko "Prediction of Cross-Fitness for Adaptive Evolution to Different Environmental Conditions: Consequence of Phenotypic Dimensional Reduction" [to be updated] and T.U.Sato, C.Furusawa and K.Kaneko "Prediction of Cross-Fitness for Adaptive Evolution to Different Environmental Conditions: Consequence of Phenotypic Dimensional Reduction" arXiv, 2022, https://arxiv.org/abs/2209.07756.

All code is written in c++ and was verified to be compiled with the Intel c++ compiler (icpc) version 2021.3.0 with the options "-std=c++11 -qopenmp" on CentOS.
Explanation for source codes is the follwing;

## Evo_lowdim_.cpp

A program for evolutionary simulation to introduce low-dimensional structure into gene regulatory networks (using in the Sec.III in the paper).

(ex.) The case with $\tilde{P}=1,\alpha=0.45$ and random_seed_genotype : 0.
```bash
icpc -std=c++11 -qopenmp Evo_lowdim_.cpp
a.out 1 0.45 0 
```

## Calc_VE_.cpp (Calc_VG_.cpp)

A program to obtain phenotypic changes due to environmental stress (genotypic mutations) in the evolved gene regulatory network obtained with "Evo_lowdim_.cpp" (using in the Sec.III B in the paper).

(ex.) The case to calculate phenotypic changes due to environmental stresses with the cells of the 2nd fittest cells. 
```bash
icpc -std=c++11 Calc_VE_.cpp
a.out 1 0.45 0 1
```


## Evo_antibiotics_.cpp

A program to simulate adaptive evolution to environmental stress from the gene regulatory network obtained with "Evo_lowdim_.cpp" (using in the Sec.IV in the paper).

(ex.) The case with $\tilde{P}=1, \alpha=0.45$, random_seed_genotype : 0 and random_seed_antibiotics : 1.
```bash
icpc -std=c++11 -qopenmp Evo_antibiotics_.cpp
a.out 1 0.45 0 1
```

## Calc_base_.cpp

A program to calculate the fitness with environmental stress in the parental network in "Evo_antibiotics_.cpp" (using in the  Sec.IV in the paper).

(ex.) The case with $\tilde{P}=1, \alpha=0.45$, random_seed_genotype : 0 and random_seed_antibiotics : 1.
```bash
icpc -std=c++11 calc_base_.cpp
a.out 1 0.45 0 1
```

## Calc_cross_.cpp

A program to calculate the cross-fitness (using in the Sec.IV in the paper).

(ex.) The case with $\tilde{P}=1, \alpha=0.45$, random_seed_genotype : 0 and random_seed_antibiotics in evolved cell : 1 and random_seed_antibiotics to calculate cross-fitness : 2.
```bash
icpc -std=c++11 -qopenmp Evo_antibiotics_.cpp
a.out 1 0.45 0 1 2
```
