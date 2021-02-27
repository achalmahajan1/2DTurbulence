# 2DTurbulence in a periodic box
This pseudo-spectral code simulates the temporal evolution of vorticity contours in 2D homogeneous isotropic turbulence using 1024 spatial grid points in both x and y direction. 
The Navier Stokes equations are written in terms of vorticity/streamfunction with no stretching term. The linear term is solved implicitly whereas the non linear (convection) term is solved explicitly.

The time evolution is performed using integrated RKW3 IMEX scheme and spatially discretized using Pseudospectral methods.


## Check out the vorticity dynamics in this video
https://www.youtube.com/watch?v=3a-gHK3Mc_k&ab_channel=AchalMahajan

If you like this code, feel free to use it though it's not the best (or optimized) version that you can find on the internet, but will highly recommend if you are just trying to scratch the surface.
