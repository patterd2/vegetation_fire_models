%% 1D Simulation of spatial savanna model with precipitation
% Explicit Euler in time, 2D trapezoidal rule in space
% the boundaries are absorbing
% NB fire is modelled as the integral of phi of G
tic;
%% Set options for plots/movies
STOCHASTIC_GRADIENT = 0; % randomly perturb the standard linear gradient
plot_gradient_info = 0; % plot the gradient dependent variables/functions over the spatial domain
DO_MOVIE = 0; % i.e. write movie to avi file for playback, else just plots
DO_MOVIE_SCATTER = 0;
SPACE_TIME_PLOT = 0;
%% Numerical method parameters
L = 5; % working on [0,L]
N = 199; % N+1 grid points
delta = L/N;  % spatial discretization parameter
h = 0.1; % time discretisation parameter
n = 201; % number of time steps
tau = (n-1)*h; % simulations time domain is [0,tau]
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x;
mu_p = @(mu, mu_s, p) mu + mu_s.*p;
nu_p = @(nu, nu_s, p) nu + nu_s.*p;
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
beta_p = @(beta, beta_s, p) beta + beta_s.*p;

J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2))/sqrt(2*pi*sigma_F^2);
J_T_fun = @(x, a, sigma_T) exp( -( (a-x).^2 )/(2*sigma_T^2))/sqrt(2*pi*sigma_T^2);
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2))/sqrt(2*pi*sigma_W^2) ;
omega = @(w_0_fun, w_1, g, t_1_fun, s_1) w_0_fun + (w_1-w_0_fun)./(1 + exp(-(g-t_1_fun)/s_1));
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Model parameters (values from PNAS paper)

p_left=0;
p_right=0;

p_0 = p_left;
p_1 = (p_right-p_left)/L;
mu = 0.1;
mu_s = 0;
nu = 0.05;
nu_s = 0;

%%

alpha = 0.6; 

%%

alpha_s = 0.3;
beta = 0.1;
beta_s = 0;
gamma = 0;
w_0_ref = 0.9;
w_1_ref = 0.4;
t_1_ref = 0.4;
s_1 = 0.01; % s_1 standard value is 0.01
f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05;% s_2 standard value is 0.05
sigma_F = 0.05; % seed dispersal radius forest trees
sigma_T = 0.05; % seed dispersal radius savanna trees
sigma_W = 0.05; % fire spread radius
%% Set up the initial distributions of the cover types on the grid
% each row is one time step of the simulation
% solution is a block matrix [LB; SOL; RB] where LB and RB are fixed
% boundary conditions corresponding to a periodic boundary
G0 = ones(1,N+1)-0.05*rand(1,N+1);
S0 = zeros(1,N+1);
T0 = zeros(1,N+1);
F0 = zeros(1,N+1)+0.05*rand(1,N+1);

F0(floor(0.2*(N+1)):floor(0.5*(N+1)))=1;
G0(floor(0.2*(N+1)):floor(0.5*(N+1)))=0;

CR = G0 + S0 + T0 + F0;
G0 = G0./CR;
S0 = S0./CR;
T0 = T0./CR;
F0 = F0./CR;
LB_G = G0(1,2:end);
RB_G = G0(1,1:end-1);
LB_S = S0(1,2:end);
RB_S = S0(1,1:end-1);
LB_T = T0(1,2:end);
RB_T = T0(1,1:end-1);
LB_F = F0(1,2:end);
RB_F = F0(1,1:end-1);
G = [LB_G G0 RB_G];
S = [LB_S S0 RB_S];
T = [LB_T T0 RB_T];
F = [LB_F F0 RB_F];
% compute the convolution for E
X = 0:delta:L;
X_L = X-L;
X_L = X_L(1,1:end-1);
X_R = X+L;
X_R = X_R(1,2:end);
E = ones(1,N+1);
temp_normalise = ones(1,N+1);
% Save tempW matrices to avoid computing them again
tempW = ones(N+1,3*N+1);
Trap = ones(1,3*N+1);
Trap(1,3*N+1)=0.5;
Trap(1,1)=0.5;
for i = 1:N+1
    tempW(i,:) =  delta*W_fun([X_L X X_R],(i-1)*delta,sigma_W);
    integrand = tempW(i,:).*(G(1,:)+gamma*(T(1,:)+S(1,:)));
    E(1,i) = sum(integrand.*Trap);
end

%% Compute the birth and mortality matrices as a function of rainfall
% this computes the standard (deterministic) gradient values
P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
mu_grad = mu_p(mu, mu_s, P_grad);
nu_grad = nu_p(nu, nu_s, P_grad);
alpha_grad = alpha_p(alpha, alpha_s, P_grad);
beta_grad = beta_p(beta, beta_s, P_grad);
%% Stochastic gradient generation
hetr = 0.1; % this is the proportion of stochastic variance around the deterministic gradient, e.g. 20%
if STOCHASTIC_GRADIENT == 1
    for k = 1:N+1
        P_grad(1,k) = triangular_dist(rand(),P_grad(1,k)-hetr,P_grad(1,k),P_grad(1,k)+hetr);
        % NB not doing a proportional perturbation for P_grad since it is
        % not necessarily nonnegative
        mu_grad(1,k) = triangular_dist(rand(),(1-hetr)*mu_grad(1,k),mu_grad(1,k),(1+hetr)*mu_grad(1,k));
        nu_grad(1,k) = triangular_dist(rand(),(1-hetr)*nu_grad(1,k),nu_grad(1,k),(1+hetr)*nu_grad(1,k));
        alpha_grad(1,k) = triangular_dist(rand(),(1-hetr)*alpha_grad(1,k),alpha_grad(1,k),(1+hetr)*alpha_grad(1,k));
        beta_grad(1,k) = triangular_dist(rand(),(1-hetr)*beta_grad(1,k),beta_grad(1,k),(1+hetr)*beta_grad(1,k));
    end
end
%% preallocaate some temp variables for efficiency
temp1 = ones(1,N+1);
temp2 = ones(1,N+1);
tempF = ones(N+1,3*N+1);
tempT = ones(N+1,3*N+1);
%% pre-calculate the 4D convolution matrices
for k = 1:N+1
    tempF(k,:) = delta*J_F_fun([X_L X X_R],(k-1)*delta,sigma_F);
    tempT(k,:) = delta*J_T_fun([X_L X X_R],(k-1)*delta,sigma_T);
end
%% The numerical scheme
for i = 2:n
    % compute convolutions for this time step
    progressbar(i,n);
    for k = 1:N+1
        integrand1 = tempF(k,:).*F(i-1,:);
        temp1(1,k) = sum(integrand1.*Trap);
        integrand2 = tempT(k,:).*T(i-1,:);
        temp2(1,k) = sum(integrand2.*Trap);
    end
    G(i,(N+1):2*N+1) = G(i-1,(N+1):2*N+1) + h*( mu_grad.*S(i-1,(N+1):2*N+1) + nu_grad.*T(i-1,(N+1):2*N+1)...
        - alpha_grad.*temp1.*G(i-1,(N+1):2*N+1) + ...
        phi(f_0_ref, f_1_ref, E(i-1,:), t_2_ref, s_2).*F(i-1,(N+1):2*N+1)...
        - beta_grad.*temp2.*G(i-1,(N+1):2*N+1) );
    S(i,(N+1):2*N+1) = S(i-1,(N+1):2*N+1) + h*( -mu_grad.*S(i-1,(N+1):2*N+1) + beta_grad.*temp2.*G(i-1,(N+1):2*N+1)...
        - alpha_grad.*temp1.*S(i-1,(N+1):2*N+1)...
        - omega(w_0_ref, w_1_ref, E(i-1,:), t_1_ref, s_1).*S(i-1,(N+1):2*N+1) );
    T(i,(N+1):2*N+1) = T(i-1,(N+1):2*N+1) + h*( - nu_grad.*T(i-1,(N+1):2*N+1) ...
        + omega(w_0_ref, w_1_ref, E(i-1,:), t_1_ref, s_1).*S(i-1,(N+1):2*N+1)...
        -alpha_grad.*temp1.*T(i-1,(N+1):2*N+1) );
    for k = 1:N+1
        integrand = tempW(k,:).*( G(i,:)+gamma*(T(i,:)+S(i,:)) );
        E(i,k) = sum(integrand.*Trap);
    end
    F(i,(N+1):2*N+1) = ones(1,N+1) - S(i,(N+1):2*N+1) - T(i,(N+1):2*N+1) - G(i,(N+1):2*N+1);
    % now need to update the extended parts of the solution
    G(i,1:N) = G(i,N+2:2*N+1);
    G(i,2*N+2:3*N+1) = G(i,(N+1):2*N);
    S(i,1:N) = S(i,N+2:2*N+1);
    S(i,2*N+2:3*N+1) = S(i,(N+1):2*N);
    T(i,1:N) = T(i,N+2:2*N+1);
    T(i,2*N+2:3*N+1) = T(i,(N+1):2*N);
    F(i,1:N) = F(i,N+2:2*N+1);
    F(i,2*N+2:3*N+1) = F(i,(N+1):2*N);
    % For numerical stability take max with zero and min with one
    % add error message here later
    G(i,:) = min(max(G(i,:),0),1);
    S(i,:) = min(max(S(i,:),0),1);
    T(i,:) = min(max(T(i,:),0),1);
    F(i,:) = min(max(F(i,:),0),1);
    E(i,:) = min(max(E(i,:),0),1);
end
% The following output is useful when trying to discern whether or not
% a solution is stationary in time
fprintf('\n');
fprintf('The maximum changes on the grid for each variable at the last time step were:\n');
fprintf(['G: ',num2str(max(abs(G(n,:)-G(n-1,:)))),'\n']);
fprintf(['S: ',num2str(max(abs(S(n,:)-S(n-1,:)))),'\n']);
fprintf(['T: ',num2str(max(abs(T(n,:)-T(n-1,:)))),'\n']);
fprintf(['F: ',num2str(max(abs(F(n,:)-F(n-1,:)))),'\n']);
toc
%% Visualise the solution...
% either in a movie...
if DO_MOVIE
    subplot_vis = struct('cdata',[],'colormap',[]);
    v = VideoWriter(['grass_front_formation_L=',num2str(L),'_n = ',num2str(n),'_h= ',num2str(h),'_delta=',num2str(delta),'_fire=',num2str(sigma_W),'_seedsT=',num2str(sigma_T),'_seedsF=',num2str(sigma_F),'.avi']);
    v.FrameRate = 4;
    open(v)
    f = figure;
    for j = 1:100:n
        f.Name = ['Simulation time: t = ', num2str((j-1)*h)];
        ax1 = subplot(2,2,1);
        plot(X,G(j,N+1:2*N+1),'LineWidth',2);
        title('Grass');
        
        ax2 = subplot(2,2,2);
        plot(X,S(j,N+1:2*N+1),'LineWidth',2);
        title('Saplings');
        
        ax3 = subplot(2,2,3);
        plot(X,T(j,N+1:2*N+1),'LineWidth',2);
        title('Savanna Trees');
        
        ax4 = subplot(2,2,4);
        plot(X,F(j,N+1:2*N+1),'LineWidth',2);
        title('Forest Trees');
        
        xlim([ax1 ax2 ax3 ax4],[0 L])
        ylim([ax1 ax2 ax3 ax4],[0 1])
        writeVideo(v,getframe(gcf))
        subplot_vis(j) = getframe(gcf);
    end
    % fig = figure;
    % movie(fig,subplot_vis,3,1)
    close(v)
    % or just plotting the solution
end
if DO_MOVIE_SCATTER
    subplot_vis = struct('cdata',[],'colormap',[]);
    v = VideoWriter(['Dirichlet_boundary_L=',num2str(L),'_n = ',num2str(n),'_h= ',num2str(h),'_delta=',num2str(delta),'_fire=',num2str(sigma_W),'_seedsT=',num2str(sigma_T),'_seedsF=',num2str(sigma_F),'.avi']);
    v.FrameRate = 4;
    open(v)
    f = figure;
    for j = 1:100:n
        f.Name = ['Simulation time: t = ', num2str((j-1)*h)];
        ax1 = subplot(2,2,1);
        scatter([X_L X X_R],G(j,:),'filled');
        title('Grass');
        
        ax2 = subplot(2,2,2);
        scatter([X_L X X_R],S(j,:),'filled','LineWidth',1);
        title('Saplings');
        
        ax3 = subplot(2,2,3);
        scatter([X_L X X_R],T(j,:),'filled','LineWidth',1);
        title('Savanna Trees');
        
        ax4 = subplot(2,2,4);
        scatter([X_L X X_R],F(j,:),'filled','LineWidth',1);
        title('Forest Trees');
        
        xlim([ax1 ax2 ax3 ax4],[0 L])
        ylim([ax1 ax2 ax3 ax4],[0 1])
        writeVideo(v,getframe(gcf))
        subplot_vis(j) = getframe(gcf);
    end
    % fig = figure;
    % movie(fig,subplot_vis,3,1)
    close(v)
    % or just plotting the solution
end
%% Space-time plot of dynamics
if SPACE_TIME_PLOT == 1
    fST = figure;
    fST.Name = 'Evolution over time';
    
    subplot(2,2,1)
    h1 = pcolor(G(:,N+1:2*N+1));
    shading interp
    title('Grass');
    colorbar
    set(h1, 'EdgeColor', 'none');
    ylabel('Time');
    caxis([0,1])
    xticks([1 floor((N+1)/4) floor((N+1)/2) floor(3*(N+1)/4) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
    
    subplot(2,2,2)
    h2 = pcolor(S(:,N+1:2*N+1));
    shading interp
    title('Saplings');
    colorbar
    set(h2, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/4) floor((N+1)/2) floor(3*(N+1)/4) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
    
    subplot(2,2,3)
    h3 = pcolor(T(:,N+1:2*N+1));
    shading interp
    title('Savanna Trees');
    colorbar
    set(h3, 'EdgeColor', 'none');
    ylabel('Time')
    caxis([0,1])
    xticks([1 floor((N+1)/4) floor((N+1)/2) floor(3*(N+1)/4) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
    
    subplot(2,2,4)
    h4 = pcolor(F(:,N+1:2*N+1));
    shading interp
    title('Forest Trees');
    colorbar
    set(h4, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/4) floor((N+1)/2) floor(3*(N+1)/4) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
end

figure(4);
imagesc(G(:,N+1:2*N+1)');
shading interp;
custom_map = [
    linspace(1,0,100)' linspace(1,0.5,100)' linspace(1,0,100)'];
colormap(custom_map);
colorbar;