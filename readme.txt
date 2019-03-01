The content of the zip file with supplemental code is as follows:

./finite_differences/cont_elast.cpp
C++ source file for finite difference simulation of an elastic shell harboring a rectangular contractile domain.

./analytics/el_anal_xy_LxLy_comments.mws
Maple file for deriving Fourier coefficients from model parameters.

./analytics/el_anal.cpp
C++ source file that uses the Fourier coefficients which can be found from ./analytics/el_anal_xy_LxLy_comments.mws to construct the deformation field by summing the Fourier series.

./finite_elements/ScriptFreeFem_LameFixedFree.edp
A script to be used with FreeFem++ for simulating an elastic shell harboring a rectangular contractile domain using finite element analysis.

./ablation/string_tri.cpp
A C++ source file to simulate ablation experiments. Must be compiled together with ./ablation/outpt.h. The ./ablation/meshfiles folder contains files serving input (initial conditions) for this simulation.