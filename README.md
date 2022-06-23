# ecFactory: A multi-step method for prediction of metabolic engineering gene targets

![fig1](https://user-images.githubusercontent.com/26483972/175298382-c3ceb172-bf59-4eb7-85e6-aa0ec60928c9.jpg)

The ecFactory method is a series of sequential steps for identification of metabolic engineering gene targets. These targets show which genes should be subject to overexpression, modulated expression (knock-down) or deletion (knock-out), with the objective of increasing production of a given metabolite. This method was developed by combining the principles of the FSEOF algorithm (flux scanning with enforced objective function) together with the features of GECKO enzyme-constrained metabolic models (ecModels), which incorporate enzymes as part of genome-scale metabolic networks.

## Required Software
* A functional Matlab installation (MATLAB 7.3 or higher). 
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

## Installation

Clone this repository into an accesible directory in your computer. No further steps are needed.

## Tutorial

A case study for prediction of metabolic engineering targets for increased production of 2-phenylethanol in *S. cerevisiae* cells using [ecYeastGEM](https://github.com/SysBioChalmers/yeast-GEM) and the **ecFactory** method is explained in detail in a MATLAB live script. To run this example, open the [live script](https://github.com/SysBioChalmers/ecFactory-case-studies/blob/main/code/find_gene_targets.mlx) in MATLAB and run it! with this, you will see the outputs of the method scripts in real time. 

All the relevant outputs of the method are stored in the `results folder` in this repository.

Last update: 2022-06-23

This repository is administered by [Iván Domenzain](https://github.com/IVANDOMENZAIN), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

