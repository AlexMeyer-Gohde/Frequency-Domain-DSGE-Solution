function marginal = marginal_density(simulations,lposterior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%marginal_density.m
%
%This program:
%               Calculates the marginal density from the posteriors using
%               Geweke's modified harmonic mean estimator
%Input:         simulations: a matrix containing the draws from the posterior
%               lposterior: a vector containing the log posterior densities
%                   associated with the draws
%Ouput:         marginal: the marginal data density
%
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%  This file is exclusively based on the .m files: marginal_density from dynare_v3 and dynare_4_0_4.
%
%  The primary difference is that this file creates a function on the desktop with all the
%  constants already calculated, circumventing the repetitious calculations
%  in the dynare source code - the savings in computational time is
%  essentially irrelevant, as the likelihood calculations are much more computationally intensive.
%
% Accordingly:
%   Copyright (C) 2003-2008 Dynare Team
%
%   This file is part of Dynare.
%
%   Dynare is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   Dynare is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
% A copy of the GNU General Public License associated with Dynare can be
% found in the directory /dynare_gpl_fdl.  
% If not, see <http://www.gnu.org/licenses/>%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[N,npara]=size(simulations);

lpost_mode = max(lposterior);

MU = mean(simulations)';
SIGMA = zeros(npara,npara);
for i=1:N;
    SIGMA = SIGMA + (simulations(i,:)'-MU)*(simulations(i,:)'-MU)';
end;
SIGMA = SIGMA/N;

DetSIGMA = det(SIGMA);
InvSIGMA = inv(SIGMA);
check_coverage = 1;
increase = 1;
attempts=0;
while check_coverage>0.01 && attempts<30
    marginal = [];
    for p = 0.1:0.1:0.9;
        critval = chi2inv(p,npara);
        tmp = 0;
        for i = 1:N;
            deviation  = (simulations(i,:)-MU')*InvSIGMA*((simulations(i,:)-MU'))';
            if deviation <= critval;
                lftheta = -log(p)-(npara*log(2*pi)+log(DetSIGMA)+deviation)/2;
                tmp = tmp + exp(lftheta - lposterior(i)+lpost_mode);
            end;    
        end;
        marginal = cat(1,marginal,[p,lpost_mode-log(tmp/N)]); 
    end;    
    check_coverage=abs((marginal(9,2)-marginal(1,2))/marginal(9,2));
    attempts=attempts+1;
    increase = 1.2*increase;
    invSIGMA = inv(SIGMA*increase);      
    detSIGMA = det(SIGMA*increase);
end
marginal=mean(marginal(:,2));