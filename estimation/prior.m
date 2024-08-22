function [prior] = prior_text(Theta)
% Construct the prior densities
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
% Based on file priors.m
% Associated with "Bayesian Estimation of DSGE Models" Herbst and
% Schorfheide
%
%  The priors considered are:
%  1) tau is GAMMA with mean 2 and st.d. 0.50
%  2) kappa is UNIFORM with [0,1]
%  3) psi1 is GAMMA with mean 1.50 and st.d 0.25
%  4) psi2 is GAMMA with mean 0.50 and st.d.0.25
%  5) rA is GAMMA with mean 0.50 and st.d. 0.50
%  6) piA is GAMMA with mean 7.00 and variance 2.00
%  7) gammaQ is Normal with mean 0.40 and variance 0.20
%  8) rho_R is UNIFORM with [0,1)
%  9) rho_g is UNIFORM with [0,1)
% 10) rho_z is UNIFORM with [0,1)
% 11) sigma_R is InvGAMMA with mean 0.40 and st.d. 4.00
% 12) sigma_g is InvGAMMA with mean 1.00 and st.d. 4.00
% 13) sigma_z is InvGAMMA with mean 0.50 and st.d. 4.00

%Prior from uniform pdf
a = [0, 0, -1, -1];
b = [1, 1, 1, 1];

P2 = 1/(b(1)-a(1));
P8 = 1/(b(2)-a(2));
P9 = 1/(b(3)-a(3));
P10 = 1/(b(4)-a(4));


%Prior from gamma pdf
para1 = [2, 1.5, 0.5, 2, 7];
para2 = [0.5, 0.25, 0.25, 0.5, 2];

b = para2.^2./para1;
a = para1./b;
P1 = gampdf(Theta(1),a(1),b(1));
P3 = gampdf(Theta(3),a(2),b(2));
P4 = gampdf(Theta(4),a(3),b(3));
P5 = gampdf(Theta(5),a(4),b(4));
P6 = gampdf(Theta(6),a(5),b(5));

% Prior from normal pdf
P7 = normpdf(Theta(7),0.40,0.20);

% Prior from inverse Gamma
P11 = exp(lnpdfig(Theta(11),0.40,4));
P12 = exp(lnpdfig(Theta(12),1.00,4));
P13 = exp(lnpdfig(Theta(13),0.50,4));


f = P1*P2*P3*P4*P5*P6*P7*P8*P9*P10*P11*P12*P13;

%if f <= 10^(-13)
%   prior = -Inf;
%else
   prior =log(f);
%end

end
function y = lnpdfig(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);
end




