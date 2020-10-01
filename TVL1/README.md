## Credits
In our implement, we use the following packages developed by other researchers:

the "TFOCS" package by Becker et al. (https://github.com/cvxr/TFOCS);

the "toolbox_signal" and "toolbox_general" packages by Peyre (https://github.com/gpeyre/numerical-tours/tree/master/matlab).

The input image in the TV-L1 experiment comes from http://www.hlevkin.com/TestImages/man.bmp.


## Instruction

1. Unzip the three zip files here in this folder.

2. Open Matlab

3. Run set_up.m

4. There are 7 algorithms to test:

     run PDHG.m to test the vanilla PDHG;

     run DP_PDHG.m to test the diagonal-preconditioned PDHG;

     run iPrePDHG_BCD.m to test iPrePDHG with BCD as the inner loops;

     run iPrePDHG_FISTA.m to test iPrePDHG with FISTA as the inner loops;

     run ADMM_quasiexact to test PrePDHG with nearly exact solving of subproblems;

     run APDHG to test accelerated PDHG;

     run ALADMM to test accelrated (primal) linearized ADMM.

  If an algorithm requires more than 500s to finish, we will add annotation in the beginning of its corresponding file.
