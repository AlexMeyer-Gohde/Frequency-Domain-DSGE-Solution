function [input] = covs_vec_new(input)
% Calculates covariances and variance decomposition for the model in input
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



CORRELATION_HORIZON=40;


A=input.A;C=input.C;B=input.B; delta=input.BlkTol;
D=input.D;
D_z=input.D_z;
D_z_vec=input.D_z_vec;
Y_0=input.Y_0;
Omega=input.Omega;

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



grid_size=max(2^10,2^(nextpow2(2*n_y)+1));
%grid_size=2^16;
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;
s_Y=zeros((n_y+n_w)^2,grid_size);
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
temp=[temp;D_z(e_minus)/e_minus];
%temp=(1/e_minus*A+B+e_minus*C)\(1/e_minus*A*X_0-D*reshape(D_zs(:,n), [n_w,n_w]));
temp=temp*Omega*temp';
s_Y(:,n)=temp(:);
 end
s_Y=ifft(transpose(s_Y)/(2*pi),'symmetric')*2*pi;
%s_Y=reshape(transpose(s_Y), [n_y n_y*grid_size]);


for n = 1 : grid_size
s_temp = reshape(s_Y(n,:),(1)*n_y+(1)*n_w,(1)*n_y+(1)*n_w);
s_YY{n}= s_temp;
end
for j=1:length(s_YY)+1
if j-length(s_YY)/2<=0
Gamma{j}=s_YY{j+length(s_YY)/2};
else
Gamma{j}=s_YY{j-length(s_YY)/2};
end
end

input.Gamma=Gamma;
input.cov=Gamma{length(s_YY)/2+1};

for i=1:(1)*n_y+(1)*n_w
factor(:,:,i)=((Gamma{length(s_YY)/2+1}(i,i).^(1/2))*(diag(Gamma{length(s_YY)/2+1}).^(1/2)));
end
for j=1:length(Gamma)
for i=1:(1)*n_y+(1)*n_w
Gamma_xcorr{i}(:,j)=Gamma{j}(:,i)./factor(:,:,i);
end
end

% Gamma_xxcorr=zeros((1)*n_y+(0)*n_w,(1)*n_y+(0)*n_w,length(Gamma));
% for j=1:length(Gamma)
% for i=1:(1)*n_y+(0)*n_w
% Gamma_xxcorr(:,j,i)=Gamma{j}(:,i)./factor(:,:,i);
% end
% end

for j=1:2*CORRELATION_HORIZON+1
for i=1:(1)*n_y+(1)*n_w
    xcorr{j}(:,i)=Gamma_xcorr{i}(:,grid_size/2+1-CORRELATION_HORIZON+j-1);
end
end
input.xcorr=xcorr;

R = chol(Omega);
for k=1:n_w
    select=zeros(n_w,n_w);
    select(k,k)=1;
    s_Y=zeros((n_y+n_w)^2,grid_size);
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
temp=[temp;reshape(D_zs(:,n)/e_minus, [n_w,n_w])];
%temp=(1/e_minus*A+B+e_minus*C)\(1/e_minus*A*X_0-D*reshape(D_zs(:,n), [n_w,n_w]));
temp=temp*R'*select*R*temp';
s_Y(:,n)=temp(:);
 end
s_Y=ifft(transpose(s_Y)/(2*pi),'symmetric')*2*pi;
%s_Y=reshape(transpose(s_Y), [n_y n_y*grid_size]);


for n = 1 : grid_size
s_temp = reshape(s_Y(n,:),(1)*n_y+(1)*n_w,(1)*n_y+(1)*n_w);
s_YY{n}= s_temp;
end
for j=1:length(s_YY)+1
if j-length(s_YY)/2<=0
Gamma{j}=s_YY{j+length(s_YY)/2};
else
Gamma{j}=s_YY{j-length(s_YY)/2};
end
end
cov_temp=Gamma{length(s_YY)/2+1};
decomp(:,k)=diag(cov_temp)./diag(input.cov);
end
input.decomp=decomp;
