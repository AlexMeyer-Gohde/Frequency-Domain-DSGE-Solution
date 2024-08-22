% Display some posterior results
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
% Based on file disp_results.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide
[Nsim,Npam] = size(Thetasim);

[yy, yy05, yy95]=moment(Thetasim(:,1:13));



%% description
    
sum_vec = [yy' yy05' yy95'];
vartype     = {'\tau','\kappa','\psi_1','\psi_2','r^{(A)}',...
               '\pi^{(A)}','\gamma^{(Q)}',...
               '\rho_{r}','\rho_{g}', '\rho_{z}', ...
               '\sigma_{r}','\sigma_{g}', '\sigma_{z}'};
      
col1=[yy; yy05];
col1=col1(:);
col2=[NaN(size(yy)); yy95];
col2=col2(:);

final_table=[NaN(size(col1)),col1,col2 ];

marginal = marginal_density(Thetasim,logposterior);