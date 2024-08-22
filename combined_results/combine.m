% Combine the posteriors from the different specifications
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


clear all
load('../example_ar/bayesian_imp_short.mat','bayesian')
i=1
combined_bayesian_corr_05(:,:,:,i)= bayesian.corr_05;
combined_bayesian_corr_16(:,:,:,i)= bayesian.corr_16;
combined_bayesian_corr_50(:,:,:,i)= bayesian.corr_50;
combined_bayesian_corr_84(:,:,:,i)= bayesian.corr_84;
combined_bayesian_corr_95(:,:,:,i)= bayesian.corr_95;
combined_bayesian_irf_05(:,:,:,i)= bayesian.irf_05;
combined_bayesian_irf_16(:,:,:,i)= bayesian.irf_16;
combined_bayesian_irf_50(:,:,:,i)= bayesian.irf_50;
combined_bayesian_irf_84(:,:,:,i)= bayesian.irf_84;
combined_bayesian_irf_95(:,:,:,i)= bayesian.irf_95;
bayesian_decomp_05(:,:,i)=bayesian.decomp_05;
bayesian_decomp_16(:,:,i)=bayesian.decomp_16;
bayesian_decomp_50(:,:,i)=bayesian.decomp_50;
bayesian_decomp_84(:,:,i)=bayesian.decomp_84;
bayesian_decomp_95(:,:,i)=bayesian.decomp_95;
load('../example_ma/bayesian_imp_short.mat','bayesian')
i=2
combined_bayesian_corr_05(:,:,:,i)= bayesian.corr_05;
combined_bayesian_corr_16(:,:,:,i)= bayesian.corr_16;
combined_bayesian_corr_50(:,:,:,i)= bayesian.corr_50;
combined_bayesian_corr_84(:,:,:,i)= bayesian.corr_84;
combined_bayesian_corr_95(:,:,:,i)= bayesian.corr_95;
combined_bayesian_irf_05(:,:,:,i)= bayesian.irf_05;
combined_bayesian_irf_16(:,:,:,i)= bayesian.irf_16;
combined_bayesian_irf_50(:,:,:,i)= bayesian.irf_50;
combined_bayesian_irf_84(:,:,:,i)= bayesian.irf_84;
combined_bayesian_irf_95(:,:,:,i)= bayesian.irf_95;
bayesian_decomp_05(:,:,i)=bayesian.decomp_05;
bayesian_decomp_16(:,:,i)=bayesian.decomp_16;
bayesian_decomp_50(:,:,i)=bayesian.decomp_50;
bayesian_decomp_84(:,:,i)=bayesian.decomp_84;
bayesian_decomp_95(:,:,i)=bayesian.decomp_95;
load('../example_log_lag/bayesian_imp_short.mat','bayesian')
i=3
combined_bayesian_corr_05(:,:,:,i)= bayesian.corr_05;
combined_bayesian_corr_16(:,:,:,i)= bayesian.corr_16;
combined_bayesian_corr_50(:,:,:,i)= bayesian.corr_50;
combined_bayesian_corr_84(:,:,:,i)= bayesian.corr_84;
combined_bayesian_corr_95(:,:,:,i)= bayesian.corr_95;
combined_bayesian_irf_05(:,:,:,i)= bayesian.irf_05;
combined_bayesian_irf_16(:,:,:,i)= bayesian.irf_16;
combined_bayesian_irf_50(:,:,:,i)= bayesian.irf_50;
combined_bayesian_irf_84(:,:,:,i)= bayesian.irf_84;
combined_bayesian_irf_95(:,:,:,i)= bayesian.irf_95;
bayesian_decomp_05(:,:,i)=bayesian.decomp_05;
bayesian_decomp_16(:,:,i)=bayesian.decomp_16;
bayesian_decomp_50(:,:,i)=bayesian.decomp_50;
bayesian_decomp_84(:,:,i)=bayesian.decomp_84;
bayesian_decomp_95(:,:,i)=bayesian.decomp_95;
load('../example_log_harmonic_lag/bayesian_imp_short.mat','bayesian')
i=4
combined_bayesian_corr_05(:,:,:,i)= bayesian.corr_05;
combined_bayesian_corr_16(:,:,:,i)= bayesian.corr_16;
combined_bayesian_corr_50(:,:,:,i)= bayesian.corr_50;
combined_bayesian_corr_84(:,:,:,i)= bayesian.corr_84;
combined_bayesian_corr_95(:,:,:,i)= bayesian.corr_95;
combined_bayesian_irf_05(:,:,:,i)= bayesian.irf_05;
combined_bayesian_irf_16(:,:,:,i)= bayesian.irf_16;
combined_bayesian_irf_50(:,:,:,i)= bayesian.irf_50;
combined_bayesian_irf_84(:,:,:,i)= bayesian.irf_84;
combined_bayesian_irf_95(:,:,:,i)= bayesian.irf_95;
bayesian_decomp_05(:,:,i)=bayesian.decomp_05;
bayesian_decomp_16(:,:,i)=bayesian.decomp_16;
bayesian_decomp_50(:,:,i)=bayesian.decomp_50;
bayesian_decomp_84(:,:,i)=bayesian.decomp_84;
bayesian_decomp_95(:,:,i)=bayesian.decomp_95;
save combined_results
