% Optimizes an approximation of the mode and constructs the candidate for Metropolis Haastings
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
% Based on file Candidate.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide:
%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%        Constructing the Candidate Density for MH Algorithm
%
%
%
% Author: Minsu Chang        minsuc@sas.upenn.edu
% Modified from Candidate.m by Luigi Bocola
%==========================================================================


clear
clc
close all
delete *.asv

tic

l = path;

% path('../Herbst_schofheide/DSGE Estimation/Mfiles',path);
% path('../Herbst_schofheide/DSGE Estimation/Optimization Routines',path);
%path('../Herbst_schofheide/DSGE Estimation/LRE',path);
%path('../Herbst_schofheide/DSGE Estimation/Matfiles',path);
path('Matfiles',path);
path('../algorithm',path);
path('../Optimization Routines',path);
path('../estimation',path);

disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: THE CANDIDATE DENSITY      ');
disp('                                                                  ');

% load('MH_candidate.mat')
% param=mode;
%load('mhdraws.mat')
%param=mean(Thetasim(40001:end,:));
param = [2.09 0.6530 2.00 0.65 0.34 3.16 0.51 1.5 4 2.5 0.19 0.65 0.5]; % starting values
param = [4.09 0.6530 2.00 0.65 0.76 3.08 0.55 1.5 10 10 0.19 0.65 0.5]; % starting values
param=[3      0.5       1.6      0.72 ...
       1.8       4.5     0.5      2 ...
      -4       -4       .68307729       .47189761 ...
     0.22179423 ];
param = [3 -0.75 1.6 0.75 1.8 4.5 0.5 0 0 0 1 5 5]; % starting values
param = [4 0.4 1.3 0.3 1.2 3.1 0.52 0.77 -0.7 -0.1 0.4 1 0.22]; for ii=[2 8]; param(ii)=log(param(ii))-log(1-param(ii)); end; for ii=[9 10]; param(ii)=log(1+param(ii))-log(1-param(ii)); end% starting values
%param = [3 0.4 1.3 0.3 1.2 3.1 0.52 0.77 0 0  1 1 1]; for ii=[2 8]; param(ii)=log(param(ii))-log(1-param(ii)); end; for ii=[9 10]; param(ii)=log(1+param(ii))-log(1-param(ii)); end% starting values
%param=[2.91251920000000	0.483084225219626	1.50576900000000	0.0163358840000000	1.21450340000000	2.94751200000000	0.544500510000000	0.646670225384118	-0.9	0.000237525865532939	2.20499140000000	1.00944560000000	0.263870990000000];
%for ii=[2 8]; param(ii)=log(param(ii))-log(1-param(ii)); end; for ii=[9 10]; param(ii)=log(1+param(ii))-log(1-param(ii)); end
param=[4.2637  0.3794  1.2620  0.1933  1.3145  3.2525  0.4550  0.7652 -0.7043 -0.2029  0.4113  0.8103  0.6362]; 
param=[4.0621  0.3358  1.2875  0.2024  1.3746  3.4902  0.4592  0.7659 -0.6940 -0.1901  0.4590  0.8312  0.4437];
param=[3.5  0.4  1.2  0.34  1.2  3.1  0.5  0.7659 -0.8 -0.901  2.7  1.2  2];
for ii=[2 8]; param(ii)=log(param(ii))-log(1-param(ii)); end; for ii=[9 10]; param(ii)=log(1+param(ii))-log(1-param(ii)); end% starting values


[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('objectiveconstr',param',eye(13)*1e-8,[] ,10^(-8),200);
%options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',1e10,'MaxIterations',1e5);
%x = fminunc('objectiveconstr',param',options);



Theta=x';

Theta(2) = exp(Theta(2))/(1+exp(Theta(2)));
Theta(8) = exp(Theta(8))/(1+exp(Theta(8)));

Theta(9) = 1-(2)/(1+exp(Theta(9)));
Theta(10) = 1-(2)/(1+exp(Theta(10)));



%disp('   tau    kappa    psi1    psi2    rA    piA    gammaQ    rho_R    rho_g    rho_z    sigma_R    sigma_g    sigma_z   ');
%disp(num2str(Theta))
%disp('                                                                  ');                 
  

mode  = Theta; 
Sigma = nhess(@objectiveunconstr,mode');
Sigma = inv(Sigma);


save Matfiles/MH_candidate Sigma mode 

path(l);

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;



