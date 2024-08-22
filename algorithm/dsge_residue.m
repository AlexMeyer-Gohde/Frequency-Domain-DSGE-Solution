function [input] = dsge_residue(input)
% The main algorithm. Solves the model in input using residues
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Checks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=input.A;B=input.B;C=input.C; delta=input.BlkTol; D=input.D;
[A_row,A_col] = size(A);[C_row,C_col] = size(C);[B_row,B_col] = size(B);
[D_row,D_col] = size(D);
if A_row~=A_col || B_row~=B_col || C_row~=C_col ||  A_row~=B_col || A_row~=C_col
    display('Error')
    display('A, B, and C must be square matrices of the same dimension')
    output=[];
    return
end
n_y=A_row;
n_w=D_col;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the companion linearization
A=[eye(n_y,n_y),zeros(n_y,n_y);zeros(n_y,n_y),A];
B=[zeros(n_y,n_y),eye(n_y,n_y);-C,-B];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generalized Schur
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S,T,Q_star,M] = qz(A,B,'complex');
eigs=ordeig(T,S);
[eigs eigs_index]=sort(eigs,'ComparisonMethod','abs');
if abs(eigs(n_y+1))<1
    %disp('Too few unstable roots')
    input.existence_uniqueness=2;
elseif abs(eigs(n_y))>1
    %disp('Too few unstable roots')
    input.existence_uniqueness=3;
else
    input.existence_uniqueness=1;
end
[~,ord] = sort(eigs_index);
ord = max(ord)-ord+1;
[S,T,Q_star,M] = ordqz(S,T,Q_star,M,ord);
input.eigs=ordeig(T,S);
% Determine reordering of lower half of Schur form into block form.
lower_half_ordering = blocking(S(n_y+1:end,n_y+1:end),T(n_y+1:end,n_y+1:end),delta);
[ord, ind] = swapping(lower_half_ordering);  % Gives the blocking.
ord = max(ord)-ord+1;        % Since ORDSCHUR puts highest index top left.
final_ord=[2*n_y:-1:n_y+1, ord];
[S,T,Q_star,M] = ordqz(S,T,Q_star,M,final_ord);
S_u=S(n_y+1:end,n_y+1:end);
T_u=T(n_y+1:end,n_y+1:end);
eigs_u=ordeig(T_u,S_u);
D_tilde=Q_star*[zeros(n_y,n_y);eye(n_y,n_y)]*D;
D_u=D_tilde(n_y+1:end,:);
U_0=NaN(n_y,n_w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Derivatives of z W(z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
highest_d_derivative=max(cellfun(@length, ind));
if sum(abs(eigs_u(ind{end}))>1e16)>0
    highest_d_derivative=highest_d_derivative+1;
end
%matlabFunction(sym('z')*sym(input.D_z))
if highest_d_derivative>1
    D_z{1}=input.D_z;
    if isfield(input,"D_z_i") && length(input.D_z_i)>=highest_d_derivative-1
        D_z(2:1+length(input.D_z_i))=input.D_z_i;
    else
        D_z_start=sym(D_z);
        for j=2:highest_d_derivative
            D_z_start=diff(D_z_start);
            D_z_i{j-1}=matlabFunction(D_z_start);
        end
        D_z(2:1+length(D_z_i))=D_z_i;
    end
else
    D_z=input.D_z;
end


U_u=NaN(n_y,n_w,highest_d_derivative,length(ind)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Main algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=length(ind):-1:1 %Go through the blocks, starting at the last
    k=length(ind{m}); %Size of the block m
    % Set matrices for block m
    S_m=S_u(ind{m},ind{m});
    T_m=T_u(ind{m},ind{m});
    S_hat_m=triu(S_m,1);
    T_hat_m=triu(T_m,1);
    S_m_plus=S_u(ind{m},max(ind{m})+1:end);
    T_m_plus=T_u(ind{m},max(ind{m})+1:end);
    % Calculate residues
    if sum(abs(eigs_u(ind{m}))>1e16)>0 %"Infinite eigenvalues"
        Theta=S_hat_m/T_m;
        %T_D=T_m\D_u(ind{m},:);
        U_0(ind{m},:)=(1/k)*D_u(ind{m},:)*D_z{k+1}(0);
        for jj=k-1:-1:1
            U_0(ind{m},:)=(1/jj)*(Theta*U_0(ind{m},:)+D_u(ind{m},:)*D_z{jj+1}(0));
        end
        U_0(ind{m},:)=T_m\U_0(ind{m},:);
    else % Finite eigenvalues
        mu=1/eigs_u(ind{m});mu=mu(1);
        Theta=(S_hat_m-mu*T_hat_m)/T_m;
        T_D=T_m\D_u(ind{m},:);
        U_0(ind{m},:)=zeros(length(ind{m}),n_w);
        for jj=k-1:-1:1
            U_0(ind{m},:)=(1/jj)*Theta*(U_0(ind{m},:)+T_D*D_z{jj+1}(mu)+(S_m_plus-mu*T_m_plus)*U_u(max(ind{m})+1:end,:,jj,m)-jj*T_m_plus*U_u(max(ind{m})+1:end,:,jj-1,m));
        end
        U_0(ind{m},:)=S_m\(U_0(ind{m},:)+D_u(ind{m},:)*D_z{1}(mu)+(S_m_plus-mu*T_m_plus)*U_u(max(ind{m})+1:end,:,1,m)-S_u(ind{m},max(ind{m})+1:end)*U_0(max(ind{m})+1:end,:));
    end
% Calculate for all blocks mm lower than m (come afterwards in recursion)
% U_m(mu_mm) and necessary derivatives i for block mm of size length(ind{mm})
    for mm=m-1:-1:1
        mu=1/eigs_u(ind{mm});mu=mu(1);
        LHS=S_m-mu*T_m;
        for i=1:length(ind{mm}) %derivatives needed for eigenvalue block mm
            if i==1
                RHS=-D_u(ind{m},:)*D_z{i}(mu)-(S_m_plus-mu*T_m_plus)*U_u(max(ind{m})+1:end,:,i,mm)+S_m_plus*U_0(max(ind{m})+1:end,:)+S_m*U_0(ind{m},:);
                U_u(ind{m},:,1,mm)=LHS\RHS;
            else
                RHS=-D_u(ind{m},:)*D_z{i}(mu)-(S_m_plus-mu*T_m_plus)*U_u(max(ind{m})+1:end,:,i,mm)+(i-1)*T_m_plus*U_u(max(ind{m})+1:end,:,i-1,mm)+(i-1)*T_m*U_u(ind{m},:,i-1,mm);
                U_u(ind{m},:,i,mm)=LHS\RHS;
            end
        end
    end
end
% Translate U_0 into Y_0
M_star=M';
M_star_22=M_star(n_y+1:end,n_y+1:end);
Y_0=M_star_22\U_0;
if max(max(abs(imag(Y_0))))<1e12
    input.Y_0=real(Y_0);
else
    disp('Imaginary components - check')
    disp(max(max(abs(imag(Y_0)))))
    input.Y_0=Y_0;
end
end

