% weighted gradient operator
% Convex Region-Based Image Segmentation
%%
clear all;
addpath('toolbox_signal') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('toolbox_general') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('functions')
addpath('TFOCS-master')
%%
name = '1211_25-Most-Beautiful-Blue-Flowers_174945887.jpg_1';
m = 660;
n = 720;
%I = rescale(crop(load_image(name),n));
I = rescale(load_image(name));
I = I(1:m, 1:n, 1:3);
clf;
imageplot(I);


%% Compute w0 and w1. Compute and display the segmentation
c0 = [0;0;1];
c1 = [0;1;0];
compute_w = @(I,c)sum( (I - repmat(reshape(c,[1 1 3]), [m n 1])).^2, 3);
w0 = compute_w(I,c0);
w1 = compute_w(I,c1);
Omega = w0<w1;
display_segmentation = @(u)repmat(reshape(c0,[1 1 3]), [m n 1]) .* repmat(u>.5, [1 1 3]) + ...
        repmat(reshape(c1,[1 1 3]), [m n 1]) .* repmat(u<.5, [1 1 3]);
clf; 
imageplot(display_segmentation(Omega));

w = w0-w1;

beta=10;
Diag=zeros(m,n,2);
Diag(:,:,1)=exp(-beta*(abs(I([2:m,1],:,1)-I(:,:,1))+...
    abs(I([2:m,1],:,2)-I(:,:,2))+abs(I([2:m,1],:,3)-I(:,:,3))));
%Diag(:,:,1)=n*n*Diag(:,:,1)./sum(sum(Diag(:,:,1)));
Diag(:,:,2)=exp(-beta*(abs(I(:,[2:n,1],1)-I(:,:,1))+...
    abs(I(:,[2:n,1],2)-I(:,:,2))+abs(I(:,[2:n,1],3)-I(:,:,3))));
%Diag(:,:,2)=n*n*Diag(:,:,2)./sum(sum(Diag(:,:,2)));

lambda = 2;
alpha = 1/lambda;

clf;
%% obtain the optimal solution
[B1,B2]=generate_B_Neumann(m,n);
B=[B2;B1];
d1=reshape(Diag(:,:,1),[m*n,1]);
d2=reshape(Diag(:,:,2),[m*n,1]);
w_r=reshape(w,[m*n,1]);



load optval; %optimal function value

 %cvx_solver gurobi 
 %cvx_begin 
     %variable x(m*n)
     %minimize norm(d1.*(B2*x),1)+norm(d2.*(B1*x),1)+alpha*sum(x.*w_r)
     %subject to
              %0<=x<=1
 %cvx_end 
 %optval = cvx_optval;
 %solution_cvx=x;
 
 %f=reshape(x,m,n);imageplot(f);
 %optval=norm(d1.*(B2*x),1)+norm(d2.*(B1*x),1)+alpha*sum(x.*w_r)


%%
Diag1=diag(sparse(d1));
Diag2=diag(sparse(d2));

load C;
load Cnorm;
% C=[Diag1*B2;Diag2*B1];
% Cnorm=normest(C*C');

%%  Convex Discrete Formulation
options.bound = 'sym';
Grad = @(x)grad(x,options);
Div = @(x) div(x,options);
Div_w=@(x)Div(Diag.*x);

ProxG = @(f,tau)max(0,min(1,f));

Amplitude = @(u)sqrt(sum(u.^2,3));
F=@(u)sum(sum(sum(abs(u))));
ProxFS = @(u,sigma)max(-1,min(u,1));
Proj= @(u) max(-1,min(u,1));








