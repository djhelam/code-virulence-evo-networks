# Landscape structure drives eco-evolution in host-parasite systems.

## Code for:

Jhelam N. Deshpande, Vasilis Dakos, Oliver Kaltz and Emanuel A. Fronhofer: Landscape structure drives eco-evolution in host-parasite systems.

## Abstract: 

As all biological and many artificial systems, hosts and their parasites are most often spatially structured. Such spatial structure becomes evident when one considers, for example, cities connected by commuters or habitat islands connected by dispersing individuals. Within this spatial context, parasites may change through time due to rapid evolutionary processes, including mutation and selection, which strongly impacts how epidemics unfold. As a consequence, host-parasite systems can only be understood in the light of evolution and in their spatial context. However, past work has largely ignored relevant spatial complexity. How parasites evolve in realistically complex landscapes remains therefore unclear, hampering the translation of theoretical predictions to real ecological systems. Therefore, we here develop an eco-evolutionary metapopulation model of host-parasite dynamics in which hosts and parasites disperse through realistically complex spatial networks. Parasite virulence, a key parasite life-history trait that impacts host fitness (here: reproduction), is able to evolve. Our model therefore captures feedback loops between host demography and parasite evolution in space. In order to gain a general understanding of parasite eco-evolution in space, we analyse our model for spatial networks that represent characteristic terrestrial (represented by random-geometric graphs; RGG) and riverine aquatic (represented by optimal channel networks; OCN) landscapes. We find that parasites may evolve higher virulence in characteristically riverine compared to terrestrial networks. More generally, we explore how a key host trait, dispersal, modulates virulence evolution within these networks. Consistent with previous work, we show that parasite virulence evolution is driven by kin selection, because dispersal and landscape structure both impact patterns of relatedness. Our model yields readily testable predictions, including that terrestrial parasites should be more virulent than aquatic parasites at low dispersal rates and vice versa as dispersal increases. These differences in evolved virulence directly lead to differences in system stability, with more virulent parasites more often leading to their own extinction. Thus, in this study we highlight the role of landscape structure in driving eco-evolutionary dynamics of parasites.

## Description:

There are two folders:

1. cpp_codes: this contains the C++ code host_parasite_network.cpp for the individual-based model of host-parasite evolution in networks, an input file input.txt corresponding to one of the scenarios explored and a Makefile.

2. cpp_codes_shuffled_simulations: this contains the C++ code to run simulations in which parasite genotypes are re-shuffled in each generation to break kin structure. This is done to test the role of kin selction in dribing virulence evolution.

3.generate_adjacency_matrices: this contains R-scripts used to generate all adjacency matrices for the networks: RGG, OCN, circle, grid, modular, spiky. 


