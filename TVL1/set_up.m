% TV-L1
clear all
M=1024;
N=1024;
MN = M*N;

% Addpath
addpath('toolbox_signal') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('toolbox_general') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('functions')
addpath('TFOCS-master')

% load an image
name = 'man';
f0 = load_image(name);
%f0(:,:) = f00(:,:,1);
f0 = rescale(crop(f0,N));
% display it
clf;
imageplot(f0);
load g0
% diaplay it
clf;
imageplot(g0); 
clf;
% Matrix
[B1,B2] = generate_B_Neumann(M,N);
B=[B1;B2];
b = zeros(MN,1);
g=reshape(g0,[MN,1]);
lambda=1;

%% Compute f_end using CVX
load x_cvx;



