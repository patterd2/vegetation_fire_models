%% Matrix simulation of 2D FGBA system with periodic boundary conditions
% Convention for state labels: 0 = Grass, 1 = Forest, 2 = Burning, 3 = Ash
clear all;

% Parameters
L = 5;            % Working on [0, L]
num_sites = 100;  % Total number of sites in [0, L]
T = 100;          % Simulation time length in real-time units
dt = 0.1;        % Length of time step in real-time units
t0 = 0;           % Start time

phi_A = 1e-1; % Rate of forest seeding into ash
phi_G = 1e-1; % Rate of forest seeding into grass

gamma_A = 1e1; % Rate of grass regrowth from ash

mu = 1e-3; % Rate of non-fire forest mortality

beta_F = 1e4; % Rate of fire spread through forest
beta_G = 1e5; % Rate of fire spread through grass
xi = 1e-4; % Rate of spontaneous fires

q = 1e5; % Fire quenching rate

sigma_F = 1; % width of Gaussian of forest seeding
sigma_B = 0.01;    % width of Gaussian for burning spread over grassland
sigma_G = 0.01;  % width of Gaussian for flammability of large grass

t2 = 0.4;     % Sigmoid center
g0 = 0.000;    % Sigmoid lower bound
g1 = 1;     % Sigmoid upper bound
s2 = 0.05;    % Sigmoid width

% Periodic distance function
dist = @(r1, r2) angle (exp(1i * (r1 - r2) * 2 * pi / L)) * L / (2 * pi); 

% Gaussian functions
W_F = @(r1, r2) exp( -( dist(r1, r2).^2 ) / (2 * sigma_F^2)) / sqrt(2 * pi * sigma_F^2);
W_B = @(r1, r2) exp( -( dist(r1, r2).^2 ) / (2 * sigma_B^2)) / sqrt(2 * pi * sigma_B^2);
W_G = @(r1, r2) exp( -( dist(r1, r2).^2 ) / (2 * sigma_G^2)) / sqrt(2 * pi * sigma_G^2);

% Sigmoid forest mortality function
phi = @(x) (g0 + (g1 - g0)./(1 + exp(-(x - t2) / s2)));

% Allocate states of sites
X = zeros(num_sites,1);

% Fix locations of sites
Locations = sort(L * rand(num_sites, 1));

% Initialize states of sites - one forest patch in [0.2L, 0.5L]
for i = 1:num_sites
    if (Locations(i) > 0.2*L) && (Locations(i) < 0.5*L) 
        X(i, 1) = 1;
    end
end

% Create matrix of Gaussians of distances between site locations
rate_matrix = zeros(num_sites, num_sites);

% Calculate rates of neighbor-dependent components
for i = 1:num_sites
    for j = 1:num_sites

        % X(i, j) = (G, F)
        if (X(i, 1) == 0) && (X(j, 1) == 1)
            rate_matrix(i, j) = phi_G * W_F(Locations(i), Locations(j));
        end

        % X(i, j) = (A, F)
        if (X(i, 1) == 3) && (X(j, 1) == 1)
            rate_matrix(i, j) = phi_A * W_F(Locations(i), Locations(j));
        end

        % X(i, j) = (F, B)
        if (X(i, 1) == 1) && (X(j, 1) == 2)
            rate_matrix(i, j) = beta_F * W_B(Locations(i), Locations(j));
        end
    end
end

% Calculate constant transition rates
const_rates = zeros(num_sites, 1);

for i = 1:num_sites
    if (X(i, 1) == 0) % X(i) = G
        const_rates(i) = xi;
    end

    if (X(i, 1) == 1) % X(i) = F
        const_rates(i) = mu;
    end

    if (X(i, 1) == 2) % X(i) = B
        const_rates(i) = q;
    end

    if (X(i, 1) == 3) % X(i) = A
        const_rates(i) = gamma_A;
    end
end

% Use Gillespie algorithm to model transitions over time
t = 0; % time
k = 0; % time index
Times = [0]; % array to record times of all events

tic;
while (t < T)
    k = k + 1;

    % Calculate rate of GB transition (with constants subtracted off)
    W_GG_sum = zeros(num_sites, 1);
    W_GB_sum = zeros(num_sites, 1);

    for i = 1:num_sites
        for j = 1:num_sites

            % X(i, j) = (G, G)
            if (X(i, k) == 0) && (X(j, k) == 0)
                W_GG_sum(i) = W_GG_sum(i) + W_G(Locations(i), Locations(j));
            end

            % X(i, j) = (G, B)
            if (X(i, k) == 0) && (X(j, k) == 2)
                W_GB_sum(i) = W_GB_sum(i) + W_B(Locations(i), Locations(j));
            end

        end
    end

    GB_no_const = phi(W_GG_sum) + beta_G * W_GB_sum;

    % Calculate event rate at each site
    eventRatePerSite = sum(rate_matrix, 2) + const_rates + GB_no_const;

    % Calculate total intensity over all sites
    totalIntensity = sum(eventRatePerSite);

    % Find time of next event
    NextEvent = -log(1 - rand()) / totalIntensity;
    t = t + NextEvent;

    if (t<T) % Event occurs before end of simulation
        % Append event time to Time array    
        Times(end + 1) = t;

        % Choose index of site where event occurs
        CDF = cumsum(eventRatePerSite) / totalIntensity;
        U = rand();
        site = 1;

        while U > CDF(site)
            site = site + 1;
        end
        
        % Append states of all sites to Solution array
        X(:, k + 1) = X(:, k);

        % Choose transition if site is G
        if (X(site, k) == 0)
            prob_GB = (GB_no_const(site) + const_rates(site)) / eventRatePerSite(site); % Probability of GB transition

            if (rand() < prob_GB)
                X(site, k + 1) = 2; % G to B
            else
                X(site, k + 1) = 1; % G to F
            end           
        end

        % Choose transition if site is F
        if (X(site, k) == 1)
            prob_FG = const_rates(site) / eventRatePerSite(site);

            if (rand() < prob_FG)
                X(site, k + 1) = 0; % F to G
            else
                X(site, k + 1) = 2; % F to B
            end            
        end
    
        % Choose transition if site is B
        if (X(site, k) == 2)
            X(site, k + 1) = 3; % fire must get quenched
        end
    
        % Choose transition if site is A
        if (X(site, k) == 3)
            prob_AG = const_rates(site) / eventRatePerSite(site);

            if (rand() < prob_AG)
                X(site, k + 1) = 0; % A to G
            else
                X(site, k + 1) = 1; % A to F
            end            
        end

    else % Finish simulation if event occurs after end time of simulation
        Times(end + 1) = T;
        X(:, k + 1) = X(:, k);
    end

    % Update rate_matrix 
    for i = 1:num_sites

        % Clear row at index site in rate matrix
        rate_matrix(site, :) = zeros(1, num_sites);

        % Recalculate rates in row
        if (X(site, k + 1) == 0) && (X(i, k + 1) == 1)
            rate_matrix(site, i) = phi_G * W_F(Locations(site), Locations(i));
        end

        if (X(site, k + 1) == 3) && (X(site, k + 1) == 1)
            rate_matrix(site, i) = phi_A * W_F(Locations(site), Locations(i));
        end

        if (X(i, k + 1) == 1) && (X(site, k + 1) == 2)
            rate_matrix(site, i) = beta_F * W_B(Locations(site), Locations(i));
        end

        % Reflect updates to column at index site
        rate_matrix(site, i) = rate_matrix(i, site);
    end

    % Update const_rates
    if (X(site, k + 1) == 0) % New state is G
        const_rates(site) = xi;
    end

    if (X(site, k + 1) == 1) % New state is F
        const_rates(site) = mu;
    end

    if (X(site, k + 1) == 2) % New state is B
        const_rates(site) = q;
    end

    if (X(site, k + 1) == 3) % New state is A
        const_rates(site) = gamma_A;
    end

end
toc

% Create array of full solution
Sol_Save = zeros(num_sites, length(0:dt:T));

for i = 1:num_sites
    Sol_Save(i, :)= interp1(Times, X(i,:), 0:dt:Times(end));
end

% Plot solution interpolated to a regualarly spaced grid
figure(1);
imagesc(Sol_Save);
custom_map = [0.4660 0.6740 0.1880
    0.0039 0.1953 0.1250
    0.8500 0.3250 0.0980
    0.5 0.5 0.5]; %light green, dark green, orange, gray colors
colormap(custom_map);
xlabel('time');
ylabel('space');
set(gca,'linewidth',2);
set(gca,'FontSize',20);
