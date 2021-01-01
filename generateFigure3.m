% This Matlab script generates Figure 3 in the paper:
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
d_values = [lambda/8 lambda/4];

%Number of elements in the horizontal and vertical dimensions
N_HV_max = 40;

%Prepare to save simulation results
distanceCovariance = zeros(N_HV_max,length(d_values));
distanceCovariance_eigenvalues = zeros(N_HV_max,length(d_values));


%% Go through the different number of elements
for ind = 1:length(d_values)
    
    %Extract the element size
    d = d_values(ind);
    
    
    for N_HV = 1:N_HV_max
        
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
        
        
        %Compute the Kronecker approximation
        R_ULA = R(1:N_HV,1:N_HV);
        R_approx = kron(R_ULA,R_ULA);
        
        
        %Compute the exact and approximate eigenvalues
        eigen = sort(eig(R),'descend');
        eigen_approx = sort(eig(R_approx),'descend');
        
        
        %Compute distances according to the definition in Footnote 2
        distanceCovariance(N_HV,ind) = 1 - trace(R*R_approx)/norm(R,'fro')/norm(R_approx,'fro');
        distanceCovariance_eigenvalues(N_HV,ind) = 1 - sum(eigen.*eigen_approx)/norm(R,'fro')/norm(R_approx,'fro');
        
    end
    
end




%% Plot simulation results
figure;
hold on; box on;
plot(1:N_HV_max,distanceCovariance(:,2),'k--','LineWidth',2);
plot(1:N_HV_max,distanceCovariance(:,1),'r-.','LineWidth',2);
plot(1:N_HV_max,distanceCovariance_eigenvalues(:,2),'k:','LineWidth',2);
plot(1:N_HV_max,distanceCovariance_eigenvalues(:,1),'r','LineWidth',2);

ylim([0 0.25]);
legend({'Full matrix ($d_{\textrm{H}}=d_{\textrm{V}}=\lambda/8$)','Full matrix ($d_{\textrm{H}}=d_{\textrm{V}}=\lambda/4$)','Eigenvalues ($d_{\textrm{H}}=d_{\textrm{V}}=\lambda/8$)','Eigenvalues ($d_{\textrm{H}}=d_{\textrm{V}}=\lambda/4$)'},'Interpreter','latex','Location','NorthWest');
set(gca,'fontsize',16);
xlabel('Number of elements (per dimension)','Interpreter','latex');
ylabel('Correlation matrix distance','Interpreter','latex');
