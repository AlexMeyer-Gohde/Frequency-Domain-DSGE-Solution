% Plot posteriors
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
% Based on file MetropolisHastings.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide


%=========================================================================
%                  FIGURE 1: RECURSIVE AVERAGES 
%=========================================================================

pnames = strvcat('\tau','\kappa', '\psi_{1}','\psi_{2}','r^{(A)}',...
    '\pi^{(A)}','\gamma^{(Q)}','\rho_{r}','\rho_{g}', '\rho_{z}', '\sigma_{r}', '\sigma_{g}', '\sigma_{z}');

figure('Position',[20,20,900,600],'Name',...
    'Recursive Averages','Color','w')



for i=1:(Npam)
        if i>8
        loca=i+2;
    else
        loca=i;
    end
subplot(ceil((Npam)/3),3,loca), plot(rmean(:,i),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
title(pnames(i,:),'FontSize',12,'FontWeight','bold');    
end



%=========================================================================
%                  FIGURE 2: POSTERIOR MARGINAL DENSITIES 
%=========================================================================

pnames = strvcat('\tau','\kappa', '\psi_{1}','\psi_{2}','r^{(A)}',...
    '\pi^{(A)}','\gamma^{(Q)}','\rho_{r}','\rho_{g}', '\rho_{z}', '\sigma_{r}', '\sigma_{g}', '\sigma_{z}');

figure('Position',[20,20,900,600],'Name',...
    'Posterior Marginal Densities','Color','w')


for i=1:(Npam)
    if i>8
        loca=i+2;
    else
        loca=i;
    end
    xmin = min(Thetasim(:,i));
    xmax = max(Thetasim(:,i));
    grid = linspace(xmin,xmax,100);
    u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
subplot(ceil((Npam)/3),3,loca), plot(grid,ksdensity(Thetasim(:,i)),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
    'Color','black','LineWidth',2.5 ), hold off
axis([xmin xmax 0 u]);
title(pnames(i,:),'FontSize',12,'FontWeight','bold');    
end





