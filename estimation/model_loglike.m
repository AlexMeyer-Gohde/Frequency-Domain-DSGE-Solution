function [loglik, varargout]=model_loglike(para)
% Prepares for the block Levinson algorithm to compute the
%               loglikelihood by generating the
%               sequence of autocovariance matrices of a model
%
% ALGORITHMS
%   Meyer-Gohde, Alexander (2024). "Solving and Analyzing DSGE Models in the Frequency Domain" 
%
% Authors: Alexander Meyer-Gohde, 08/2024
%
% Copyright (C) 2024 Alexander Meyer-Gohde
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Based on file model_loglike.m:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%model_loglike.m
%
%This program:
%               prepares the block Levinson algorithm for computing the
%               loglikelihood by calling on the model to get the recusive
%               solution and then, using spectral methods, generate the
%               sequence of autocovariance matrices
%
%Input:         parameters: a vector of parameters to be passed to the
%                   model
%               selection_vector: a vector containing an integer at the
%                   i'th entry that denotes which endogenous variable in
%                   the model file corresponds to the i'th observable
%               Data: vectorized data set
%               T: number of observations
%               K: number of observables
%               model_name: a function that executes the model file
%
%Output:        loglik: the value of the log-likelihood function with the
%                   given parameters and data set 
%
%THIS VERSION: 1.0 March 27, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tau = para(1);
kappa = para(2);
psi1 = para(3);
psi2 = para(4);
rA = para(5);
piA = para(6);
gammaQ = para(7);
rho_R = para(8);
rho_g = para(9);
rho_z = para(10);
sigma_R = para(11);
sigma_g = para(12);
sigma_z = para(13);

Data = load('us.txt');
[output] = model_solution_AMG(para);
[T K]=size(Data);
Data=Data-[gammaQ  piA piA+rA+4*gammaQ];
Data=Data';
Data=Data(:);
selection_vector=[4 2 3];
selection_scale=[1 4 4];

H = diag([(0.20*0.579923)^2  (0.20*1.470832)^2  (0.20*2.237937)^2]);


% selection_matrix_0=zeros(K,n_y);
% selection_matrix_1=zeros(K,n_y);
% selection_matrix_0(eq_y,y_t) =  1;
% selection_matrix_1(eq_y,y_t) =  -1; 
% selection_matrix_0(eq_y, z_t) = 1;
% selection_matrix_0(eq_pi,pi_t) =  4;
% selection_matrix_0(eq_ffr,R_t) = 4;

existence_uniqueness=output.existence_uniqueness;

A=output.A;C=output.C;B=output.B; delta=output.BlkTol;
D=output.D;
D_z=output.D_z;
D_z_vec=output.D_z_vec;
Y_0=output.Y_0;
Omega=output.Omega;

[A_row,A_col] = size(A);[C_row,C_col] = size(C);[B_row,B_col] = size(B);
[D_row,D_col] = size(D);
if A_row~=A_col || C_row~=C_col || B_row~=B_col || A_row~=C_col || A_row~=B_col
    display('Error')
    display('A, B, and C must be square matrices of the same dimension')
    output=[];
    return
end
n_y=A_row;
n_w=D_col;

selection_matrix=zeros(K,n_y);
selection_matrix((selection_vector-1).*K+[1:K])=selection_scale;

if isnumeric(existence_uniqueness)==1
if sum(abs(output.eigs)>1)~=n_y
%loglik = 10^50;
loglik = -10^50;
else




grid_size=max(2^10,2^(nextpow2(2*K)+1));
%grid_size=2^16;
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;
s_Y=zeros(K^2,grid_size);
E_minus=exp( -1i*grid_points);
D_zs=1./E_minus.*D_z_vec(E_minus);
D_zst=kron(eye(n_w,n_w),D)*D_zs;
A_Y_0=A*Y_0;
 for n = 1 : grid_size
e_minus=E_minus(n);
temp1=(1/e_minus*A+B+e_minus*C);
temp2=D_zst(:,n);
temp2=reshape(temp2, [n_y,n_w]);
temp2=(1/e_minus*A_Y_0-temp2);
temp=temp1\temp2;
%temp=(1/e_minus*A+B+e_minus*C)\(1/e_minus*A*X_0-D*reshape(D_zs(:,n), [n_w,n_w]));
temp=selection_matrix*temp;
temp=temp*Omega*temp';
temp=temp+H;
s_Y(:,n)=temp(:);
 end
s_Y=ifft(transpose(s_Y)/(2*pi),'symmetric')*2*pi;
s_Y=reshape(transpose(s_Y), [K K*grid_size]);

s_Y=[s_Y(1:K,1:K*T) s_Y(1:K,end-K*T+1:end)];

loglik = block_levinson_loglik(Data, s_Y, 1); 
%loglik=-loglik;
end
else
%loglik = 10^50;
loglik = -10^50;
end
