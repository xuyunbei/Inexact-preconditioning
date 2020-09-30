%%
clear all
M=256;
N=256;
MN = M*N;

% Addpath
addpath('AIRToolsII-master/testprobs')
addpath('functions')

% load a test problem 
[A,b,x_true,theta,~,R,d] = fancurvedtomo(N, 0:10:359);

[rA,cA]=size(A);
a = b;
% matrix B is the discrete gradient operator
[B1,B2] = generate_B_Neumann(M,N);
B=[B1;B2];
mu=1;
%normAB2=normest([A;B]'*[A;B]);
%normA2=normest(A'*A);
load normAB2;
load normA2;


%% Compute f_end using CVX

%cvx_solver gurobi
%cvx_begin 
    %variable x(MN)
    %minimize 1/2*sum((A*x-b).^2)+mu*norm(B*x,1)   
%cvx_end 
%optval = cvx_optval;
%solution_cvx=x;

load optval;


