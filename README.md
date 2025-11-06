# Project-paper
Code files necessary for project paper in orbital mechanics

Notes on particle class:
init_conds_cart returns first element of array with cartesian coordinates,
this works for 1 particle, but if later extend to several particles, this must
change

Scipy solver:
Does not calculate beta, only position and velocity
t0 not an argument in class, always starts at t=0 s