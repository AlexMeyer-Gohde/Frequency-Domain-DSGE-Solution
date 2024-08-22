% Generates a sample from the posterior using Metropolis Haastings
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
% Schorfheide:
%==========================================================================
%                       DSGE MODEL ESTIMATION 
%                   Metropolis-Hastings Algorithm
%                  (Small NK model in the textbook)
%
% Author: Minsu Chang        minsuc@sas.upenn.edu
% Last modified: 2/24/2016
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================
clc
clear 
close all
delete *.asv

tic

l = path;

% path('../Herbst_schofheide/DSGE Estimation/Mfiles',path);
% path('../Herbst_schofheide/DSGE Estimation/Optimization Routines',path);
%path('../Herbst_schofheide/DSGE Estimation/LRE',path);
%path('../Herbst_schofheide/DSGE Estimation/Matfiles',path);
path('Matfiles',path);
path('../algorithm',path);
path('../Optimization Routines',path);
path('../estimation',path);

disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: METROPOLIS-HASTINGS        ');
disp('                                                                  ');

%=========================================================================
%                  METROPOLIS-HASTINGS ALGORITHM 
% (Report the Acceptance Rate and Recursive Averages Every 500 draws) 
%=========================================================================

load MH_candidate

Nsim          = input('How Many Posterior Draws?:  ');
disp('                                                                  ');
disp('                                                                  ');

Sigma = dirty_hh(Sigma);
c             = 0.2;
c0            = 0.2;
Nburn         = int32(0.50*Nsim)+2;
Thetasim      = zeros(Nsim,13);

% Initialize by taking a draw from a distribution centered at mode
go_on = 0;
while go_on == 0
   Thetac = mvnrnd(mode',c0*Sigma);
    go_on = (Thetac(11)>=0)*(Thetac(12)>=0)*(Thetac(13)>=0)*(Thetac(8)>=0)*(Thetac(9)>=-1)*(Thetac(10)>=-1)*(Thetac(8)<=1)*(Thetac(9)<=1)*(Thetac(10)<=1)*(Thetac(2)<=1);  % bounds
end
Thetasim(1,:) = Thetac;

accept        = 0;
obj           = model_loglike(Thetasim(1,:)) + prior(Thetasim(1,:));
counter       = 0;
logposterior  = obj*ones(Nsim,1);

for i=1:Nsim
    
    Thetac = mvnrnd(Thetasim(i,:),c*Sigma);
    CheckBounds = (Thetac(11)>=0)*(Thetac(12)>=0)*(Thetac(13)>=0)*(Thetac(8)>=0)*(Thetac(9)>=-1)*(Thetac(10)>=-1)*(Thetac(8)<=1)*(Thetac(9)<=1)*(Thetac(10)<=1)*(Thetac(2)<=1);  % bounds

    if CheckBounds == 1 
    
       prioc = prior(Thetac);    
       likic = model_loglike(Thetac);
       objc  = prioc+likic;       
       
       if objc == -Inf
       
          Thetasim(i+1,:) = Thetasim(i,:);
          logposterior(i+1) = obj;
          
       else % objc > -Inf

          alpha = min(1,exp(objc-obj));
          u = rand(1);

          if u<=alpha
             Thetasim(i+1,:)   = Thetac;
             accept            = accept+1;
             obj               = objc;
             logposterior(i+1) = objc;
          else
             Thetasim(i+1,:)   = Thetasim(i,:);
             logposterior(i+1) = obj;
          end
          
       end % if objc == -Inf
       
    else % CheckBounds NE 1
  
       Thetasim(i+1,:) = Thetasim(i,:);
       logposterior(i+1) = obj;
       
    end  % if CheckBounds == 1

    acceptancerate     = accept/i;
    counter            = counter + 1;

    if counter==500
       disp('                                                                  ');
       disp(['                               DRAW NUMBER:', num2str(i)]         );
       disp('                                                                  ');
       disp('                                                                  ');    
       disp(['                           ACCEPTANCE RATE:', num2str(acceptancerate)]);
       disp('                                                                  ');
       disp('                                                                  ');    
       disp('                            RECURSIVE AVERAGES                    ');
       disp('                                                                  ');
       disp('   Tau       Kappa      Psi1       Psi2        rA        piA       gammaQ       rho_R       rho_g       rho_z       sigma_R       sigma_g       sigma_z   ');
       disp(num2str(mean(Thetasim(1:i,:))));  
       disp('                                                                  ');
       counter = 0;
    end % if counter==500
    
end %for i=1:Nsim

Thetasim    = Thetasim(Nburn:end,:);

logposterior= logposterior(Nburn:end);

save Matfiles/mhdraws Thetasim logposterior   % Save posterior draws

[Nsim,Npam] = size(Thetasim);

[yy, yy05, yy95]=moment(Thetasim(:,1:13));



%% description
    
sum_vec = [yy' yy05' yy95'];
vartype     = {'\tau','\kappa','\psi_1','\psi_2','r^{(A)}',...
               '\pi^{(A)}','\gamma^{(Q)}',...
               '\rho_{r}','\rho_{g}', '\rho_{z}', ...
               '\sigma_{r}','\sigma_{g}', '\sigma_{z}'};
      
disp('=========================================================================');
disp(' Variable Name                       Mean         5%        95%         ');
disp('=========================================================================');
for hh=1:length(vartype);
    fprintf('%-30s %10.4f %10.4f %10.4f\n',vartype{hh},sum_vec(hh,1),...
        sum_vec(hh,2),sum_vec(hh,3));    
end
disp('========================================================================='); 


%=========================================================================
%                  FIGURE 1: RECURSIVE AVERAGES 
%=========================================================================

pnames = strvcat('\tau','\kappa', '\psi_{1}','\psi_{2}','r^{(A)}',...
    '\pi^{(A)}','\gamma^{(Q)}','\rho_{r}','\rho_{g}', '\rho_{z}', '\sigma_{r}', '\sigma_{g}', '\sigma_{z}');

figure('Position',[20,20,900,600],'Name',...
    'Recursive Averages','Color','w')

rmean = zeros(Nsim,Npam);

for i=1:Nsim
    rmean(i,:) = mean(Thetasim(1:i,:),1);
end

for i=1:(Npam-1)
    
subplot((Npam-1)/3,3,i), plot(rmean(:,i),'LineStyle','-','Color','b',...
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


for i=1:(Npam-1)
    xmin = min(Thetasim(:,i));
    xmax = max(Thetasim(:,i));
    grid = linspace(xmin,xmax,100);
    u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
subplot((Npam-1)/3,3,i), plot(grid,ksdensity(Thetasim(:,i)),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
    'Color','black','LineWidth',2.5 ), hold off
axis([xmin xmax 0 u]);
title(pnames(i,:),'FontSize',12,'FontWeight','bold');    
end


disp('                                                                  ');
disp(['                     ELAPSED TIME:   ', num2str(toc)]             );

elapsedtime=toc;

path(l);


