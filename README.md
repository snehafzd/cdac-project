# **DFT study of MB(mythlene blue) dye on Graphene-ZnO Nano hybrid**

This project is designed to integrate the use of HPC in the world of material science. MB dye known to be very dangerous industry by product for animals and humans. it is not only cause water pollution but it is hard to degrade MB dye because of pi bonds. 

For this project, we inspired from [ZnO-graphene nanohybrids for photocatalytic degradation of methylene blue dye](https://www.sciencedirect.com/science/article/abs/pii/S0925963525008489). Mr. rao achieved up to 93 % degradation of a 5 ppm MB solution
within 130 min, this indicate that graphene-zno is highly efficient in solving the water pollution due to MB dye.

In our project we aimed to built reaction path and simulation of the G-ZnO using DFT in quantum espresso tool. For this,we first built the molecule and .in file of  graphene, zno(wurtzite), mb by using ASE module in pyhton using parallel programming (multiprocessing).Further, we feed those input file to quantum espresso, which was preinstalled in HPC(param utkarsh provided by CDAC-BANGLORE) using slurm and linux pipeline, we used #SBATCH -array to parallelize our work for different set of molecules in slurm. We were aiming to calculate the SCF, Relax and Band calculation for MB, Graphene, ZnO(wurtzite), G-ZnO and MB on G-ZnO. 

 
### Technology Used:<br>
DFT,Quantum Espresso, Pyhton(ASE, multiprocessing, numpy,py3Dmol),linux, Slurm, HPC






Though our efforts does not provided satisfied results. But it was worth a try to implement and execute the tools like quantum espresso, slurm , linux pipeline and python module like ASE, py3Dmol etc. for the first time.
