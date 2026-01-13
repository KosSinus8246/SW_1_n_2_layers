Linear Rotating Shallow Water codes in Fortran 90 !

this repository contains the 1-layer and 2-layers versions of the LRSW model. 3 different simulations are possible with a gaussian 2D perturbation :

1) Nothing : we let the perturbation propagate itself.
2) Island : we put at the center of a grid a square that simulates an island.
3) Detroit : we put 2 wall face to face to simulate a detroit (e.g. Gibraltar)

Integrator : Leap-Frog with centered finite differences

BC's : Rigid wall for u,v and 0-Neumann for eta
