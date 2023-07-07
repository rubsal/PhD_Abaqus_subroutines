** MATERIALS
** 
*Material, name=Durable
*Density
 1.12e-1
*User Material, constants=30
**     muA,       lamLA,       muB,      kapB,    alpha,    dgam0,          m,      tau0,
        6.3,       1.22,        262,      2500.,      0.1,    1e-4,       6.7,      15.3,  
**    tauss,      h            chi,      maxit,   tol         Gf,    delta_tc,    sigc_m,         
      13.,       300.,        0.01,     100,      1e-4,      0.5,        1e-4,      1300,  
**    sigc_std, sigc_min,  sigc_max,    seed,     sflag,     xmin,        xmax,      ymin,
        22,        40,        1000,        1,         1,        0,          5,         0,                  
**    ymax,      zmin,        zmax,       nx,        ny,      nz,           
        30,         0,          5,        24,        150,       24,            
*Depvar, DELETE=32
 32
1, FiB_11,
2, FiB_22,
3, FiB_33,
4, FiB_12,
5, FiB_23,
6, FiB_31,
7, FiB_21,
8, FiB_32,
9, FiB_13,
10, NiB_11,
11, NiB_22,
12, NiB_33,
13, NiB_12,
14, NiB_23,
15, NiB_31,
16, NiB_21,
17, NiB_32,
18, NiB_13,
19, iter
20, res
21, tauhat
22, tauOld
23, dgamb
24, Nsub
25, sigc
26, epseq
27, epseq_c
28, epseq_u
29, delta_t_sig
30, damage
31, conFlag
32, failFlag