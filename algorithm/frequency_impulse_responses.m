function [output] = frequency_impulse_responses(output)
% Calculates irfs using and ifft for the model in output
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

% selection_matrix=zeros(K,n_y);
% selection_matrix((selection_vector-1).*K+[1:K])=selection_scale;

R = chol(Omega);

if isnumeric(existence_uniqueness)==1





grid_size=max(2^10,2^(nextpow2(2*n_y)+1));

   s_Y=zeros(n_y*n_w,grid_size);
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;
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
temp=temp*R';
s_Y(:,n)=temp(:);
 end
 s_Y=real(ifft(transpose(s_Y)))';
 s_Y=reshape(s_Y, [n_y  n_w grid_size]);
 s_Z=D_zs.*R(:);%./E_minus;
 s_Z=real(ifft( transpose(s_Z)))';
 s_Z=reshape(s_Z, [n_w n_w grid_size]);
 output.irf=zeros(n_w+n_y,n_w,grid_size);
 output.irf(1:n_y,:,:)=s_Y;
 output.irf(n_y+1:end,:,:)=s_Z;
 output.irf = permute(output.irf,[1 3 2]);
end