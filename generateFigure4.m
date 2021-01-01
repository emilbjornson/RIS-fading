% This Matlab script generates Figure 4 in the paper:
%
% Emil Björnson, Luca Sanguinetti, “Rayleigh Fading Modeling and Channel
% Hardening for Reconfigurable Intelligent Surfaces,” IEEE Wireless
% Communications Letters, To appear.
%
% Download article: https://arxiv.org/pdf/2009.04723.pdf
%
% This is version 1.0 (Last edited: 2021-01-01)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.


close all;
clear;

%Select if (a) or (b) should be generated
directPath = false;

%Wavelength
lambda = 0.1; %In meters

%The width and height of an RIS element
d = lambda/4;

%Number of elements in the horizontal and vertical dimensions
N_max = 40;

%Set the strength of the direct path
if directPath == true
    
    betad = db2pow(-130);
    
else
    
    betad = 0;
    
end

%Number of Monte Carlo trials
realizations = 50000;

%Prepare to save simulation results
values_RIS = zeros(realizations,N_max);
values_noOpt = zeros(realizations,N_max);


%Area of an element
A = d^2;

%Set the average intensity attenuations
mu1 = db2pow(-55);
mu2 = mu1;

%Set transmit power in dBm
PdBm = 30;

%Set the noise power in dBm
sigma2dBm = -174 + 10*log10(10e6) + 10;

%Compute the transmit power over the noise power in linear scale
Psigma2 = db2pow(PdBm - sigma2dBm);


%% Go through the different number of elements
for N_HV = 1:N_max
    
    %Output simulation progress
    disp(['Size ' num2str(N_HV) ' out of ' num2str(N_max)])
    
    
    %Generate a grid for the elements
    gridPoints = (0:N_HV-1)*d;
    
    [X,Y] = meshgrid(gridPoints,gridPoints);
    
    locations = X(:)+1i*Y(:);
    
    
    %Total number of elements
    N = length(locations);
    
    %Compute the spatial correlation matrix
    R = zeros(N,N);
    
    for m = 1:N
        for l = 1:N
            
            R(m,l) = sinc(2*abs(locations(m)-locations(l))/lambda);
            
        end
    end
    
    %Generate channel realizations
    Rsqrtm = sqrtm(R);
    hd = sqrt(betad) * (randn(1,realizations) + 1i*randn(1,realizations))/sqrt(2);
    h1 = sqrt(A*mu1) * Rsqrtm * (randn(N,realizations) + 1i*randn(N,realizations))/sqrt(2);
    h2 = sqrt(A*mu2) * Rsqrtm * (randn(N,realizations) + 1i*randn(N,realizations))/sqrt(2);
    
    %Compute the SNR with an optimized RIS
    values_RIS(:,N_HV) = Psigma2*(sum(abs(h1.*h2),1)+abs(hd)).^2;
    
    %Compute the SNR with a random RIS configuration (see Footnote 3)
    values_noOpt(:,N_HV) = Psigma2*abs(sum(h1.*h2,1)+hd).^2;
    
end


%% Compute error bars
medianValues = zeros(N_max,length(d));
lowerLimits = zeros(N_max,length(d));
upperLimits = zeros(N_max,length(d));

medianValues_noOpt = zeros(N_max,length(d));
lowerLimits_noOpt = zeros(N_max,length(d));
upperLimits_noOpt = zeros(N_max,length(d));

for N_HV = 1:N_max
    
    %Compute an interval that contains 90% of the realizations
    valuesSorted = sort(values_RIS(:,N_HV),'ascend');
    medianValues(N_HV) = median(valuesSorted);
    lowerLimits(N_HV) = medianValues(N_HV) - valuesSorted(round(0.05*length(valuesSorted)));
    upperLimits(N_HV) = valuesSorted(round(0.95*length(valuesSorted))) - medianValues(N_HV);
    
    
    %Compute an interval that contains 90% of the realizations
    valuesSorted = sort(values_noOpt(:,N_HV),'ascend');
    medianValues_noOpt(N_HV) = median(valuesSorted);
    lowerLimits_noOpt(N_HV) = medianValues_noOpt(N_HV) - valuesSorted(round(0.05*length(valuesSorted)));
    upperLimits_noOpt(N_HV) = valuesSorted(round(0.95*length(valuesSorted))) - medianValues_noOpt(N_HV);
    
end


%% Plot simulation results

N_HVrange = (1:N_max)';


figure;
hold on; box on;
errorbar(N_HVrange,medianValues(:,1),lowerLimits(:,1),upperLimits(:,1),'b-','LineWidth',2);
errorbar(N_HVrange,medianValues_noOpt(:,1),lowerLimits_noOpt(:,1),upperLimits_noOpt(:,1),'r--','LineWidth',2);
plot(1:N_max,Psigma2*A^2*mu1*mu2*(pi^2/16)*ones(1,N_max).*N_HVrange.^4,'k:','LineWidth',3);
errorbar(N_HVrange,medianValues(:,1),lowerLimits(:,1),upperLimits(:,1),'b-','LineWidth',2);
set(gca,'fontsize',16);
set(gca,'Yscale','log');
xlabel('Number of elements (per dimension)','Interpreter','latex');
ylabel('SNR','Interpreter','latex');
legend({'Optimized phases','Random phases','Asymptotic expression'},'Interpreter','latex','Location','SouthEast');
