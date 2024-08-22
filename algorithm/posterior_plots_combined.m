% Plots posterior irfs for all four specifications combined
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

ENDOGENOUS_VARIABLE_NAMES={'Output Gap'
'Inflation'
'Nominal Interest Rate'
'Output Growth'
'Exp Shocks'
'Productivity'
'Government Expenditures'
'Monetary Policy'};

EXOGENOUS_VARIABLE_NAMES={'Productivity'
'Government Expenditures'
'Monetary Policy'};

n_y=8;
n_w=3;

scale=1;
horizon=40;
x_lim_low = 1; 
x_lim_up  = horizon;


RGB_wide = [.85 .85 .85]; % color for the outer area
RGB_edge_wide = [.75 .75 .75]; % color for the edge of the outer area
RGB_middle = [.7 .7 .7]; % color for the midlle area
RGB_edge_middle = [.6 .6 .6]; % color for edge of the middle area
RGB_close = [.6 .6 .6]; % color for the inner area
RGB_edge_close = [.5 .5 .5]; % color for the edge of the inner area

tt=(x_lim_low:1:x_lim_up);
o =   zeros(1,size(tt,2)); 
patch_time = [tt flipdim(tt,2)]'; % runs through time and back




exo_select=[1 2 3];
variable_select=[1:4, 6:8];


CC=[0	                0.447000000000000	0.741000000000000
    0.915294117647059	0.281568627450980	0.287843137254902
    0.281568627450980	0.715294117647059	0.287843137254902
    1	                0.598431372549020	0.200000000000000];

for i=1:length(exo_select)
    j=exo_select(i);
    for k=1:length(variable_select)
        figure;
        for ll=1:4
        %        if min(scale*combined_bayesian_irf_05((j-1)*(M_.endo_nbr)+variable_select(k),1:horizon)')~=max(scale*combined_bayesian_irf_95((j-1)*(M_.endo_nbr)+variable_select(k),1:horizon)')
        
        %Median and zero line:
        plot(tt, scale*combined_bayesian_irf_50(variable_select(k),1:horizon,j,ll),'color', CC(ll,:), 'Linewidth', 1.5, 'HandleVisibility','off'); hold on
        plot(tt, o, 'k', 'Linewidth', 1.5, 'HandleVisibility','off');
        Hpatch  = [patch_time [scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j,ll)';flipdim(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j,ll),2)']];
        patch(Hpatch(:,1),Hpatch(:,2),RGB_wide,'edgecolor',RGB_edge_wide, 'HandleVisibility','off'); hold on;

        Hpatch2  = [patch_time [scale*combined_bayesian_irf_16(variable_select(k),1:horizon,j,ll)';flipdim(scale*combined_bayesian_irf_84(variable_select(k),1:horizon,j,ll),2)']];
        patch(Hpatch2(:,1),Hpatch2(:,2),RGB_middle,'edgecolor',RGB_edge_middle, 'HandleVisibility','off');
        end
        %Median and zero line:
        for ll=1:4
        plot(tt, scale*combined_bayesian_irf_50(variable_select(k),1:horizon,j,ll), 'color', CC(ll,:), 'Linewidth', 2);
        end
        plot(tt, o, 'k', 'Linewidth', 1.5, 'HandleVisibility','off');
        % Recover vertical axis:
        % line([x_lim_low x_lim_low],[scale*combined_bayesian_irf_05(variable_select(k),1,j) scale*combined_bayesian_irf_95(variable_select(k),1,j)], 'Color','k', 'HandleVisibility','off');
        % if scale>0
        %     y_lim_low=((sign(min(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)'))-1)*1.1/(-2)+(sign(min(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)'))+1)*0.9/2)*min(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)');
        %     y_lim_up=((sign(max(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)'))-1)*0.9/(-2)+(sign(max(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)'))+1)*1.1/2)*max(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)');
        % else
        %     y_lim_low=((sign(min(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)'))-1)*1.1/(-2)+(sign(min(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)'))+1)*0.9/2)*min(scale*combined_bayesian_irf_95(variable_select(k),1:horizon,j)');
        %     y_lim_up=((sign(max(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)'))-1)*0.9/(-2)+(sign(max(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)'))+1)*1.1/2)*max(scale*combined_bayesian_irf_05(variable_select(k),1:horizon,j)');
        % end
        % axis([x_lim_low x_lim_up y_lim_low y_lim_up ])
        eval(sprintf('title(''Impulse Responses of %s to a Shock in %s'')',ENDOGENOUS_VARIABLE_NAMES{variable_select(k)},EXOGENOUS_VARIABLE_NAMES{exo_select(i)}))
        grid on;
        set(gca,'FontSize',8); %,'FontWeight','bold');
       % ylabel(ylabel_string,'FontSize',12);
        xlabel('Quarters','FontSize',12);
        %        end
        legend('AR (1)', 'MA(1)', 'Log Lag', 'Log Harmonic Lag')
    end
end