function [output] = model_solution_AMG(para)
% Returns the solution of a DSGE model
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
% Based on file model_solution.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide:
% Solve the Small-Scale DSGE model in the textbook 
% INPUT
% para:  structural parameters
%
% OUTPUT
% para:  structural parameters
%
% DATE: 2/13/2016
% Written by Minsu Chang


%=========================================================================
%                     Paramaters
%=========================================================================

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

bet = 1/(1+rA/400);


%=========================================================================
%                        SOLVE DSGE MODEL
%=========================================================================




% /* Variable indices */

%=========================================================================
%                       DEFINE OBJECTS
%=========================================================================

% Equation indices

eq_1    = 1;  %** (2.1) on \hat{y}(t) **/
eq_2    = 2;  %** (2.1) on \hat{pi}(t) **/
eq_3    = 3;  %** (2.1) on \hat{R}(t) **/
eq_4    = 4;  %** obs  \Delta\hat{y}(t) **/
eq_5    = 5;  %** obs  E_t shock combination

% Variable indices 

y_t   = 1;
pi_t   = 2;
R_t   = 3;
obs_delta_y_t=4;
Ezg_t1=5;

%Variable names
input.ENDOGENOUS_VARIABLE_NAMES={'Output Gap'
'Inflation'
'Nominal Interest Rate'
'Output Growth'
'Exp Shocks'};

       
%Shock indices (eps)

z_t = 1;
g_t = 2;
R_sh = 3;

input.EXOGENOUS_VARIABLE_NAMES={'Productivity'
'Government Expenditures'
'Monetary Policy'};

%SUMMARY
n_y = 5;
n_w = 3;


% /** initialize matrices **/
input.A = zeros(n_y,n_y);
input.B = zeros(n_y,n_y);
input.C = zeros(n_y,n_y);
input.D = zeros(n_y,n_w);


input.Omega=diag([sigma_z sigma_g sigma_R ].^2);
input.BlkTol=1e-10;
input.PERIOD=4;


%=========================================================================
%                EQUILIBRIUM CONDITIONS: CANONICAL SYSTEM
%=========================================================================

%=========================================================================
%         1. 
%=========================================================================

input.B(eq_1,y_t)   =  1;
input.B(eq_1,R_t)   =  1/tau;
input.D(eq_1,g_t)  = -1;
input.A(eq_1, y_t) = -1;
input.A(eq_1, pi_t) = -1/tau;
input.A(eq_1,Ezg_t1)  = -1;

%=========================================================================
%         2. 
%=========================================================================

input.B(eq_2,y_t)   =  -kappa;
input.B(eq_2,pi_t)   = 1;
input.D(eq_2,g_t)   =  kappa;
input.A(eq_2, pi_t) = -bet;

%=========================================================================
%         3. 
%=========================================================================

input.B(eq_3,y_t)   = -(1-rho_R)*psi2;
input.B(eq_3,pi_t)  = -(1-rho_R)*psi1;
input.B(eq_3,R_t)  = 1;
input.D(eq_3,g_t) = (1-rho_R)*psi2;
input.C(eq_3,R_t) = -rho_R;
input.D(eq_3,R_sh) = -1;

%=========================================================================
%         4. 
%=========================================================================

input.B(eq_4,obs_delta_y_t)   = 1;
input.B(eq_4,y_t)   = -1;
input.C(eq_4,y_t)   = 1;
input.B(eq_4,z_t)   = -1;

%=========================================================================
%         5. 
%=========================================================================

input.B(eq_5,Ezg_t1)   = 1;
input.D(eq_5,z_t)   = -1/tau;
input.D(eq_5,g_t)   = 1;

  
%=========================================================================
%         Exogenous Process in z 
%=========================================================================

% input.D_z_0=@(z)[1./(1-rho_z.*z) 0 0; 0 1./(1-rho_g.*z) 0; 0 0 1];
input.D_z=@(z)[z./(1-rho_z.*z) 0 0; 0 z./(1-rho_g.*z) 0; 0 0 z];
input.D_z_i{1}=@(z)[1./(rho_z.*z - 1).^2 0 0; 0 1./(rho_g.*z - 1).^2 0; 0 0 1];
input.D_z_i{2}=@(z)[-(2*rho_z)./(rho_z.*z - 1).^3 0 0; 0 -(2*rho_g)./(rho_g.*z - 1).^3  0; 0 0 0];
input.D_z_i{3}=@(z)[(6*rho_z.^2)./(rho_z.*z - 1).^4 0 0; 0 (6*rho_g.^2)./(rho_g.*z - 1).^4 0; 0 0 0];
      
% f_string = func2str(input.D_z);
% f_string = strrep(f_string, ',', ';');
% eval(['input.D_z_vec_new=',f_string,';']);
% 
f_string = func2str(input.D_z);
f_string = strrep(f_string, ',', ';');
f_string = strrep(f_string, ';', '+0.*z;');
f_string = strrep(f_string, ']', '+0.*z]');
eval(['input.D_z_vec=',f_string,';']);

% f_string = func2str(input.D_z_0);
% f_string = strrep(f_string, ',', ';');
% f_string = strrep(f_string, ';', '+0.*z;');
% f_string = strrep(f_string, ']', '+0.*z]');
% eval(['input.D_z_0_vec=',f_string,';']);

% 
% highest_d_derivative=4;
% input.D_z=matlabFunction(sym('z')*sym(input.D_z));
% D_z_start=sym(input.D_z);
% for j=2:highest_d_derivative
%     D_z_start=diff(D_z_start);
%     D_z_i{j-1}=matlabFunction(D_z_start);
% end
% input.D_z_i=D_z_i;


%=========================================================================
%           QZ(generalized Schur) decomposition by GENSYS
%=========================================================================

[output] = dsge_residue(input);
output.para=para;
end
