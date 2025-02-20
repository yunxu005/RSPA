# RSPA
A unified surface energy model for predicting micro-mechanics of heterogeneous composites from elasticity to fracture

A reference implementation of the unified surface energy model to fracture, as described in:  "A unified surface energy model for predicting micro-mechanics of heterogeneous composites from elasticity to fracture" 

To run the code and the notebooks the recommended steps are the following:

1.Download and install mef90 at https://github.com/bourdin/mef90.
 
2. run Fig03a.m to get the Figure3a.

3. run Fig03b.m to get the Figure3b.

4. run Fig04a.m to get the Figure4a.

5. run Fig04b.m to get the Figure4b.

6.for Figure06
  1) download directory Figure06
  2) run compute_parameter.m to generate 
  3) visulization the result file Figure06_0*_out.gen by Vislt.

7.for Figure07
   1)download directory Figure07
   2)run compute_parameter.m to generate 
   3)visulization the result file Figure07_0*_out.gen by Vislt.

8.for figure 8 
   1）download directory Figure08
   2）mpirun -np 1 vDefP -prefix Figure08_input -options_file_yaml Shear.yaml
   3）visulization the result file Figure08_out.gen by Vislt.
