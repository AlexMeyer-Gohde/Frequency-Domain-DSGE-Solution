% Calculate posterior statistics like irfs using every 20th draw
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


load('mhdraws.mat');
posterior_sample=length(Thetasim):-20:1;
%posterior_sample=posterior_sample(1:1000);

for j=1:length(posterior_sample)
    clear output
    posterior_selection=posterior_sample(j);
    para=Thetasim(posterior_selection,:);
    [output] = model_solution_AMG(para);
    [output] = frequency_impulse_responses(output);
    bayesian.irf(:,:,:,j)=output.irf;
    [output] = covs_vec_new(output);
    bayesian.cov(:,:,j)=output.cov;
    bayesian.decomp(:,:,j)=output.decomp;
    for i=1:40
    bayesian.corr(:,:,i,j)=output.xcorr{40+i};
    end
end

bayesian.cov_05= prctile(bayesian.cov, 5,3);
bayesian.cov_16= prctile(bayesian.cov, 16,3);
bayesian.cov_50= prctile(bayesian.cov, 50,3);
bayesian.cov_84= prctile(bayesian.cov, 84,3);
bayesian.cov_95= prctile(bayesian.cov, 95,3);

bayesian.corr_05= prctile(bayesian.corr, 5,4);
bayesian.corr_16= prctile(bayesian.corr, 16,4);
bayesian.corr_50= prctile(bayesian.corr, 50,4);
bayesian.corr_84= prctile(bayesian.corr, 84,4);
bayesian.corr_95= prctile(bayesian.corr, 95,4);

bayesian.decomp_05= prctile(bayesian.decomp, 5,3);
bayesian.decomp_16= prctile(bayesian.decomp, 16,3);
bayesian.decomp_50= prctile(bayesian.decomp, 50,3);
bayesian.decomp_84= prctile(bayesian.decomp, 84,3);
bayesian.decomp_95= prctile(bayesian.decomp, 95,3);

bayesian.irf_05= prctile(bayesian.irf, 5,4);
bayesian.irf_16= prctile(bayesian.irf, 16,4);
bayesian.irf_50= prctile(bayesian.irf, 50,4);
bayesian.irf_84= prctile(bayesian.irf, 84,4);
bayesian.irf_95= prctile(bayesian.irf, 95,4);