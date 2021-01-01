% This Matlab script generates Figure 2 in the paper:
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

%Wavelength
lambda = 1; %Value doesn't matter

%Set the width and height of an RIS element
d_values = [lambda/8 lambda/4 lambda/2];

%Number of elements in the horizontal and vertical dimensions
N_HV = 40;

%Prepare to save simulation results
eigenvalues = zeros(N_HV^2,length(d_values));


%% Go through the different number of elements
for ind = 1:length(d_values)
    
    %Extract the element size
    d = d_values(ind);
    
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
    
    %Compute the sorted eigenvalues
    eigenvalues(:,ind) = sort(eig(R),'descend');
    
end


%% Plot simulation results
figure;
hold on; box on;
plot(1:N,eigenvalues(:,1),'b-','LineWidth',2);
plot(1:N,eigenvalues(:,2),'r-.','LineWidth',2);
plot(1:N,eigenvalues(:,3),'k--','LineWidth',2);
plot(1:N,ones(1,N),'k:','LineWidth',2);
set(gca,'Yscale','log');
ylim([1e-4 1e2]);

rankAsympt1 = round(pi*(N_HV*d_values(1))^2);
rankAsympt2 = round(pi*(N_HV*d_values(2))^2);
rankAsympt3 = round(pi*(N_HV*d_values(3))^2);
plot(rankAsympt1,eigenvalues(rankAsympt1,1),'bo','LineWidth',2);
plot(rankAsympt2,eigenvalues(rankAsympt2,2),'ro','LineWidth',2);
plot(rankAsympt3,eigenvalues(rankAsympt3,3),'ko','LineWidth',2);

legend({'$d_{\textrm{H}}=d_{\textrm{V}}=\lambda/8$','$d_{\textrm{H}}=d_{\textrm{V}}=\lambda/4$','$d_{\textrm{H}}=d_{\textrm{V}}=\lambda/2$','i.i.d. fading'},'Interpreter','latex','Location','NorthEast');
set(gca,'fontsize',16);
xlabel('Eigenvalue number','Interpreter','latex');
ylabel('Eigenvalue','Interpreter','latex');
