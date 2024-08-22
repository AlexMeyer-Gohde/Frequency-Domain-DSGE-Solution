function [objective] = objectiveunconstr(Theta)
% Defines the unconstrained version of the objective for the posterior
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
% Based on file objectiveunconstr.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide    
prio = prior(Theta);

if prio==-Inf
    objective = -1000000000000000;
else

 liki = model_loglike(Theta);

objective = -(liki+prio);

end

