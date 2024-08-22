function [objective] = objectiveconstr(Theta)
% Defines the constrained version of the objective for the posterior
% maximization
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
% Based on file objectiveconstr.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide
Theta(2) = exp(Theta(2))/(1+exp(Theta(2)));
Theta(8) = exp(Theta(8))/(1+exp(Theta(8)));

Theta(9) = 1-(2)/(1+exp(Theta(9)));
Theta(10) = 1-(2)/(1+exp(Theta(10)));


prio = prior(Theta);

    
liki = model_loglike(Theta);
objective = -(liki+prio);


end

