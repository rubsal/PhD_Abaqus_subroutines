# Polymer model VUMAT: hyper-viscoelasticity and stochastic stress-based fracture criterion
This repository contains the Abaqus/Explicit user subroutines for a polymer model developed as part of my PhD.
Example input files are included for:

1. A tensile test with axisymmetric solid elements
2. A quasi-static compression test of an octet lattice structure

To run the jobs, type the following command:

'abaqus double inter job=INPUT_FILE user=VUMAT.f'

where 'INPUT_FILE' is the name of the abaqus input file without the '.inp' extension. 

If you want to run the jobs with multiple cpus, add the following line to the command above:

'cpus=CPU_NUM mp_mode=threads'

where 'CPU_NUM' is the number of cpus you want to use and 'mp_mode=threads' signalises a threads-based parallel execution.