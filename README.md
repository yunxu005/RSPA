# RSPA
A reference implementation of the unified surface energy model to fracture, as described in:  "A unified surface energy model for predicting micro-mechanics of heterogeneous composites from elasticity to fracture" 

To run the code and the notebooks the recommended steps are the following:

1. Download and install mef90 at https://github.com/bourdin/mef90.
 
2. Run Fig03a.m to get Figure 3a.

3. Run Fig03b.m to get Figure 3b.

4. Run Fig04a.m to get Figure 4a.

5. Run Fig04b.m to get Figure 4b.

6. For Figure06

  1) Download directory Figure06
  
  2) Run compute_parameter.m to generate parameters in Tension.yaml

  3) mpirun -np 1 vDefP -prefix Figure06_input -options_file_yaml Tension.yaml

  4) Change the loading velocity to obtain Figure06_0*_out.gen under different velocities
  
  5) Visulization the result file Figure06_0*_out.gen

7. For Figure07
  
  1) Download directory Figure07
  
  2) Run compute_parameter.m to generate parameters in Tension.yaml

  3) mpirun -np 1 vDefP -prefix Figure07_input -options_file_yaml Tension.yaml

  4) Change the loading velocity to obtain Figure07_0*_out.gen under different velocities
  
  5) Visulization the result file Figure07_0*_out.gen

8. For Figure08 
 
  1) Download directory Figure08
  
  2) mpirun -np 1 vDefP -prefix Figure08_input -options_file_yaml Shear.yaml
   
  3) Visulization the result file Figure08_out.gen
