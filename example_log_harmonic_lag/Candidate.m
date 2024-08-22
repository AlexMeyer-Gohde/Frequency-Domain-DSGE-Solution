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
param = [2.09 0.6530 2.00 0.65 0.76 3.08 0.55 1.5 -1 -1 0.19 0.65 0.5]; % starting values
param = [3	0.65	1.6	0.3	1.7	3.5	0.5	5	5	5	0.06	1.9	0.2];
param = [2.82 0.78 1.8 0.63 0.42 3.3 0.52 0.77 0.98 0.88 0.22 0.71 0.31]; for ii=[2 8 9 10]; param(ii)=log(param(ii))-log(1-param(ii)); end% starting values

[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('objectiveconstr',param',eye(13)*1e-2,[] ,10^(-8),200);
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



