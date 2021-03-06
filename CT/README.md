## Credits
In our implement, we use the following packages developed by other researchers:

the "TFOCS" package by Becker et al. (https://github.com/cvxr/TFOCS);

the "AIR Tools II" package by Hansen and Jorgensen (https://github.com/jakobsj/AIRToolsII).

## Instructions
1. Unzip the two zip packages here in this folder. 

2. Open Matlab

3. Run set_up.m to load the data and the optimal function value

4. There are 6 algorithms to test:

   run PDHG.m to test the vanilla PDHG;
   
   run DP_PDHG.m to test the diagonal-preconditioned PDHG;
   
   run iPrePDHG_BCD.m to test our proposed "inexact preconditioning" algorithm;
   
   run ADMM_quasiexact to test PrePDHG with nearly exact solving of subproblems;
   
   run APDHG to test accelerated PDHG;
   
   run ALADMM to test accelrated (primal) linearized ADMM.
   
If an algorithm requires more than 500s to finish, we will add annotation in the beginning of its corresponding file.
