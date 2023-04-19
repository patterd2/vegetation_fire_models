%% Matrix simulation of Poisson Forest-Grass Markov Chain
% system is spatially extended with periodic boundary conditions (circle)

clear all;

% RNG seed
%rng(1546793)

Times=[];
DO_PLOT=0;
No_Abs=1;

% Spatial Parameters
L = 5; % working on [0,L]
dist = @(a,x) angle (exp(1i*(x-a)*2*pi/L))*L/(2*pi);
J_F_fun = @(x, a, sigma_F) exp( -( dist(a,x).^2)/(2*sigma_F^2) )/sqrt(2*pi*sigma_F^2);
W_fun = @(x, a, sigma_W) exp( -( dist(a,x).^2 )/(2*sigma_W^2) )/sqrt(2*pi*sigma_W^2);

sigma_F = 0.05; % seed dispersal radius forest trees
sigma_W = 0.05; % fire radius

% Model Parameters
alpha=0.6; 

t2=0.4;
f0=0.1;
f1=0.9;
s2=0.05;

sites=500;

P=1;                % Number of Patches
N=sites*ones(1,P);    % Number of Sites / Patch
NTot=sum(N);        % Total number of sites

T=20; % length of the simulation in real-time units
dt=0.01;
t0= 0;

MC=1;

J=alpha;
W=1;
phi=@(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2))); % sigmoid definition
Ntimes=length(t0:dt:(T-dt));

tic;

nbins=100;
Histo=zeros(nbins,MC);
p_grass=1;
p_no_grass=0;
Sol_Save=zeros(MC,sites,length(0:dt:T));
%% This block of code runs one simulation of the particle system using the Gillespie algorithm
for iteration=1:MC
    %fprintf('Monte Carlo Simulation #=%d/%d\n',iteration,MC);
    Solution=zeros(NTot,1);
    Locations=sort(L*rand(sites,1));
    [A,B]=meshgrid(Locations);
    [A_left,B_left]=meshgrid(Locations-L);
    [A_right,B_right]=meshgrid(Locations+L);
%     J_Mat = (J_F_fun(A,B,sigma_F) + J_F_fun(A,B_left,sigma_F) + J_F_fun(A,B_right,sigma_F))/sites;
%     W_Mat = (W_fun(A,B,sigma_W) + W_fun(A,B_left,sigma_W) + W_fun(A,B_right,sigma_W))/sites;
% %     
%     J_Mat = (J_F_fun(A,B,sigma_F) )/sites;
%     W_Mat = (W_fun(A,B,sigma_W) )/sites;
    J_Mat = ones(sites)/sites;
    W_Mat = ones(sites)/sites;
    
    % check that the normalization of the kernels is correct
    %[X,Y] = meshgrid(0:1/sites:L);
    %[X_L,Y_L] = meshgrid(-L-(1/sites):1/sites:-(1/sites));
    %[X_R,Y_R] = meshgrid(L+(1/sites):1/sites:2*L+(1/sites));
    %J_Mat_regular = (J_F_fun(X,Y,sigma_F) + J_F_fun(X,Y_L,sigma_F) + J_F_fun(X,Y_R,sigma_F) )/sites;
    %plot(sum(J_Mat_regular));
    
    Solution(:,1)=rand(NTot,1)<p_grass;
    No_grass=find((Locations>0.2*L).*(Locations<0.5*L));
    N_no_grass=length(No_grass);
    Solution(No_grass,1)=rand(N_no_grass,1)<p_no_grass;
    
    Times=[0];
    
    k=0;
    t=0;
    while (t<T)
        k=k+1;
        BirthRates=alpha*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        %BirthRates=alpha*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        %DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        
        totalIntensity=sum(BirthRates+DeathRates);
        
        NextEvent=-log(1-rand())/totalIntensity;
        
        t=t+NextEvent;
        if (t<T)
            Times(end+1)=t;

            CDF=cumsum(BirthRates+DeathRates)/totalIntensity;
            U=rand();
            i=1;
            while U>CDF(i)
                i=i+1;
            end
            Solution(:,k+1)=Solution(:,k);
            Solution(i,k+1)=1-Solution(i,k);
        else
            Times(end+1)=T;
            Solution(:,k+1)=Solution(:,k);
        end
    end
%     Solution(:,end)=Solution(:,k);
    for i=1:sites
        Sol_Save(iteration,i,:)= interp1(Times,Solution(i,:),0:dt:Times(end));
    end
end
toc

%% Plots
figure(1);
imagesc(squeeze(mean(Sol_Save,1)));
custom_map = [1 1 1
    0 0.5 0];
colormap(custom_map);

Z=100;

V=squeeze(mean(Sol_Save,1));
U=squeeze(mean(reshape(V,sites/Z,Z,[]),1));
figure;
[x,y]=meshgrid(0:dt:Times(end),5*(1:Z)/Z);
pcolor(x,y,flipud(U))
shading interp;
custom_map = [
    linspace(1,0,100)' linspace(1,0.5,100)' linspace(1,0,100)'];
colormap(custom_map);