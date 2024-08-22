% Produce the table with the posterior variance decompositions
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

[n_y n_w n_k]=size(bayesian_decomp_50);

for j=1:4%n_y
    for k=1:n_k
        for i=1:n_w
        final_decom_table((j-1)*(n_k*2)+2*k-1,2*i-1)=bayesian_decomp_50(j,i,k);
        final_decom_table((j-1)*(n_k*2)+2*k-1,2*i)=NaN;
        final_decom_table((j-1)*(n_k*2)+2*k,2*i-1)=bayesian_decomp_05(j,i,k);
        final_decom_table((j-1)*(n_k*2)+2*k,2*i)=bayesian_decomp_95(j,i,k);
        end
    end
end
num2str(final_decom_table,' %.4f')