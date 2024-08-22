function [y,y05,y95] = moment(IRF)
% Calculates the mean, 5% and 95% pointwise credible set boundaries
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
% Based on file moment.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide
[nsim,H] = size(IRF);
y = (mean(IRF));
y05 = zeros(1,H);
y95 = zeros(1,H);
for i=1:H
 
  d= sort(IRF(:,i)); 
  y05(:,i) = d(round(0.05*nsim),:);
  y95(:,i) = d(round(.95*nsim),:);
end
end
