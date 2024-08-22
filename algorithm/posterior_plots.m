% Plots posterior irfs 
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

A=output.A;C=output.C;B=output.B; delta=output.BlkTol;
D=output.D;
D_z=output.D_z;
D_z_vec=output.D_z_vec;
Y_0=output.Y_0;
Omega=output.Omega;

[A_row,A_col] = size(A);[C_row,C_col] = size(C);[B_row,B_col] = size(B);
[D_row,D_col] = size(D);
if A_row~=A_col || C_row~=C_col || B_row~=B_col || A_row~=C_col || A_row~=B_col
    display('Error')
    display('A, B, and C must be square matrices of the same dimension')
    output=[];
    return
end
n_y=A_row;
n_w=D_col;

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



bayesian.irf_05= prctile(bayesian.irf, 5,4);
bayesian.irf_16= prctile(bayesian.irf, 16,4);
bayesian.irf_50= prctile(bayesian.irf, 50,4);
bayesian.irf_84= prctile(bayesian.irf, 84,4);
bayesian.irf_95= prctile(bayesian.irf, 95,4);

exo_select=[1 2 3];
variable_select=[1:4, 6:8];


for i=1:length(exo_select)
    j=exo_select(i);
    for k=1:length(variable_select)
        %        if min(scale*bayesian.irf_05((j-1)*(M_.endo_nbr)+variable_select(k),1:horizon)')~=max(scale*bayesian.irf_95((j-1)*(M_.endo_nbr)+variable_select(k),1:horizon)')
        figure;
        %Median and zero line:
        plot(tt, scale*bayesian.irf_50(variable_select(k),1:horizon,j), 'b', 'Linewidth', 1.5, 'HandleVisibility','off'); hold on
        plot(tt, o, 'k', 'Linewidth', 1.5, 'HandleVisibility','off');
        Hpatch  = [patch_time [scale*bayesian.irf_05(variable_select(k),1:horizon,j)';flipdim(scale*bayesian.irf_95(variable_select(k),1:horizon,j),2)']];
        patch(Hpatch(:,1),Hpatch(:,2),RGB_wide,'edgecolor',RGB_edge_wide, 'HandleVisibility','off'); hold on;

        Hpatch2  = [patch_time [scale*bayesian.irf_16(variable_select(k),1:horizon,j)';flipdim(scale*bayesian.irf_84(variable_select(k),1:horizon,j),2)']];
        patch(Hpatch2(:,1),Hpatch2(:,2),RGB_middle,'edgecolor',RGB_edge_middle, 'HandleVisibility','off');

        %Median and zero line:
        plot(tt, scale*bayesian.irf_50(variable_select(k),1:horizon,j), 'b', 'Linewidth', 1.5);
        plot(tt, o, 'k', 'Linewidth', 1.5, 'HandleVisibility','off');
        % Recover vertical axis:
        % line([x_lim_low x_lim_low],[scale*bayesian.irf_05(variable_select(k),1,j) scale*bayesian.irf_95(variable_select(k),1,j)], 'Color','k', 'HandleVisibility','off');
        % if scale>0
        %     y_lim_low=((sign(min(scale*bayesian.irf_05(variable_select(k),1:horizon,j)'))-1)*1.1/(-2)+(sign(min(scale*bayesian.irf_05(variable_select(k),1:horizon,j)'))+1)*0.9/2)*min(scale*bayesian.irf_05(variable_select(k),1:horizon,j)');
        %     y_lim_up=((sign(max(scale*bayesian.irf_95(variable_select(k),1:horizon,j)'))-1)*0.9/(-2)+(sign(max(scale*bayesian.irf_95(variable_select(k),1:horizon,j)'))+1)*1.1/2)*max(scale*bayesian.irf_95(variable_select(k),1:horizon,j)');
        % else
        %     y_lim_low=((sign(min(scale*bayesian.irf_95(variable_select(k),1:horizon,j)'))-1)*1.1/(-2)+(sign(min(scale*bayesian.irf_95(variable_select(k),1:horizon,j)'))+1)*0.9/2)*min(scale*bayesian.irf_95(variable_select(k),1:horizon,j)');
        %     y_lim_up=((sign(max(scale*bayesian.irf_05(variable_select(k),1:horizon,j)'))-1)*0.9/(-2)+(sign(max(scale*bayesian.irf_05(variable_select(k),1:horizon,j)'))+1)*1.1/2)*max(scale*bayesian.irf_05(variable_select(k),1:horizon,j)');
        % end
        % axis([x_lim_low x_lim_up y_lim_low y_lim_up ])
        %eval(sprintf('title(''Impulse Responses of %s to a Shock in %s'')',M_.endo_names(variable_select(k),:),M_.exo_names(j,:)))
        grid on;
        set(gca,'FontSize',14); %,'FontWeight','bold');
       % ylabel(ylabel_string,'FontSize',12);
        xlabel('Quarters','FontSize',12);
        %        end
    end
end


