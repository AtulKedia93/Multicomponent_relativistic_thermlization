# README #

This is a repostitory by University of Notre Dame Ph.D. student Atul Kedia to share the codes written on
MATLAB for simulating thermalization of nuclei present during big bang nucleosynthesis as described in :

Phys. Rev. D 101, 123532 - https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.123532, https://arxiv.org/abs/1911.07334

Phys. Rev. E 103, 032101 - https://journals.aps.org/pre/abstract/10.1103/PhysRevE.103.032101, https://arxiv.org/abs/2004.13186

Here you will find a lot of codes and data that look very similar to each other. As a general rule, please don't
use codes/data from here that are older than Dec-2019 (apart from exceptions like 'randpdf.m' mentioned below).
Those are legacy codes used to make a more current version which have a similar name. There might be some
repetition in newer codes as well, but for all, in general please use the most recent version.

The most important files for the thermalization are 

3D : code_e3D_vfv_w_random_w_xsection.m, code_e3D_vfv_w_random.m
     e3D_Rel_FDM_KT01_inject_GeV.m (for high energy injected nuclei)

2D : code_e2D_vfv_w_random_w_xsection.m, 

1D : code_e1D_vfv.m
     e1D_Rel_FDM_KT1_inject_GeV.m (for high energy injected nuclei)

Some of the 3D and 2D files typically need a module called "randpdf.m"
(Adam Nieslony (2021). Random numbers from a user defined distribution (https://www.mathworks.com/matlabcentral/fileexchange/26003-random-numbers-from-a-user-defined-distribution), MATLAB Central File Exchange.)

Total size of the repository is approx. : 110 MB with the main files listed above taking < 50 kB.

Atul - 10/March/2021