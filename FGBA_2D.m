%% Matrix simulation of 2D FGBA system with periodic boundary conditions
% Convention for state labels: 0 = Grass, 1 = Forest, 2 = Burning, 3 = Ash
clear all;

% Parameters
L = 1;            % Working on [0, L]
num_sites = 25^2; % Total number of sites in [0, L] x [0, L], must be square number
T = 0.5;            % Simulation time length in real-time units
dt = 0.01;        % Length of time step in real-time units
dx = 0.005;        % Spatial resolution in x direction 
dy = 0.005;        % Spatial resolution in y direction
t0 = 0;           % Start time

varphi_A = 10;  % Rate of forest seeding into ash
varphi_G = 10;  % Rate of forest seeding into grass
gamma = 50;     % Rate of grass regrowth from ash

beta_F = 5;    % Rate of fire spread through forest
beta_G = 100;    % Rate of fire spread through grass

mu = 0.01;        % Rate of tree death due to non-fire causes
q = 150;         % Fire quenching rate

sigma_P = 0.1;   % width of Gaussian of forest seeding
sigma_B = 0.1;    % width of Gaussian for burning spread over grassland
sigma_G = 0.1;    % width of Gaussian for flammability of large grass
sigma_F = 0.1;    % width of Gaussian for flammability of large grass

theta_f = 0.4;    % Forest burning sigmoid center
f0 = 0;        % Forest burning sigmoid lower bound
f1 = 0.5;         % Forest burning sigmoid upper bound
s_f = 0.05;       % Forest burning sigmoid width

theta_g = 0.4;    % Grass burning sigmoid center
g0 = 0.5;        % Grass burning sigmoid lower bound
g1 = 1;         % Grass burning sigmoid upper bound
s_g = 0.05;       % Grass burning sigmoid width

% Periodic distance function
magnitude = @(x, y) sqrt(x^2 + y^2);
dist = @(x1, y1, x2, y2) angle (exp(1i * magnitude(x2 - x1, y2 - y1) * 2 * pi / L)) * L / (2 * pi); 

% Gaussian functions
W_P = @(x1, y1, x2, y2) exp( -( dist(x1, y1, x2, y2).^2 ) / (2 * sigma_P^2)) / sqrt(2 * pi * sigma_P^2);
W_B = @(x1, y1, x2, y2) exp( -( dist(x1, y1, x2, y2).^2 ) / (2 * sigma_B^2)) / sqrt(2 * pi * sigma_B^2);
W_G = @(x1, y1, x2, y2) exp( -( dist(x1, y1, x2, y2).^2 ) / (2 * sigma_G^2)) / sqrt(2 * pi * sigma_G^2);
W_F = @(x1, y1, x2, y2) exp( -( dist(x1, y1, x2, y2).^2 ) / (2 * sigma_F^2)) / sqrt(2 * pi * sigma_F^2);

% Sigmoid forest mortality function
phi_F = @(x) (f0 + (f1 - f0)./(1 + exp(-(x - theta_f) / s_f)));
phi_G = @(x) (g0 + (g1 - g0)./(1 + exp(-(x - theta_g) / s_g)));

% Allocate states of sites
X = zeros(num_sites, 1);

% Fix locations of sites (x, y)
locations = L * rand(num_sites, 2);  % NOT SORTED

% Fix square root of number of sites
sqrt_num_sites = sqrt(num_sites);

% Initialize states of sites - one square forest patch at [0.2L, 0.5L] x [0.2L, 0.5L] 
for i = 1:num_sites
    if (locations(i, 1) < 0.5 * L && locations(i, 1) > 0.2 * L)
        if (locations(i, 2) < 0.5 * L && locations(i, 2) > 0.2 * L)
            X(i, 1) = 1;
        end
    end
end

% Rate of interactions between each pair of sites due to forest seeding
varphi_rates = zeros(num_sites, num_sites);

% Calculate rates of neighbor-dependent spreading components
for i = 1:num_sites
    for j = 1:num_sites

        % X(i, j) = (G, F)
        if (X(i, 1) == 0) && (X(j, 1) == 1)
            varphi_rates(i, j) = varphi_G * W_P(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2))/num_sites;
        end

        % X(i, j) = (A, F)
        if (X(i, 1) == 3) && (X(j, 1) == 1)
            varphi_rates(i, j) = varphi_A * W_P(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2))/num_sites;
        end
    end
end

% Calculate constant transition rates
const_rates = zeros(num_sites, 1);

for i = 1:num_sites
    if (X(i, 1) == 0) % X(i) = G
        const_rates(i) = 0;
    end

    if (X(i, 1) == 1) % X(i) = F
        const_rates(i) = mu;
    end

    if (X(i, 1) == 2) % X(i) = B
        const_rates(i) = q;
    end

    if (X(i, 1) == 3) % X(i) = A
        const_rates(i) = gamma;
    end
end

% Use Gillespie algorithm to model transitions over time
t = 0; % time
k = 0; % time index
Times = [0]; % array to record times of all events

tic;
while (t < T)
    k = k + 1;

    % Calculate rate of G -> B transition
    W_GG_sum = zeros(num_sites, 1);
    W_GB_sum = zeros(num_sites, 1);
    GB_no_const = zeros(num_sites, 1);

    for i = 1:num_sites
        for j = 1:num_sites

            % X(i, j) = (G, G)
            if (X(i, k) == 0) && (X(j, k) == 0)
                W_GG_sum(i) = W_GG_sum(i) + W_G(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end

            % X(i, j) = (G, B)
            if (X(i, k) == 0) && (X(j, k) == 2)
                W_GB_sum(i) = W_GB_sum(i) + W_B(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end

        end

        % Compute G -> B transition rate at G sites
        if (X(i, k ) == 0)
            GB_no_const(i) = phi_G(W_GG_sum(i)/num_sites) + beta_G * W_GB_sum(i)/num_sites;
        end

    end

    % Calculate rate of F -> B transition
    W_FF_sum = zeros(num_sites, 1);
    W_FB_sum = zeros(num_sites, 1);
    FB_no_const = zeros(num_sites, 1);

    for i = 1:num_sites
        for j = 1:num_sites

            % X(i, j) = (F, G)
            if (X(i, k) == 1) && (X(j, k) == 0)
                W_FF_sum(i) = W_FF_sum(i) + W_F(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end

            % X(i, j) = (F, B)
            if (X(i, k) == 1) && (X(j, k) == 2)
                W_FB_sum(i) = W_FB_sum(i) + W_B(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end

        end

        % Compute F -> B transition rate at F sites
        if (X(i, k) == 1)
            FB_no_const(i) = phi_F(W_FF_sum(i)/num_sites) + beta_F * W_FB_sum(i)/num_sites;
        end
    end

    % Calculate event rate at each site
    eventRatePerSite = sum(varphi_rates, 2) + const_rates + GB_no_const + FB_no_const;

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
            prob_GB = GB_no_const(site) / eventRatePerSite(site); % Probability of G -> B transition

            if (rand() < prob_GB)
                X(site, k + 1) = 2; % G to B
            else
                X(site, k + 1) = 1; % G to F
            end           
        end

        % Choose transition if site is F
        if (X(site, k) == 1)
            prob_FG = const_rates(site) / eventRatePerSite(site); % Probability of F -> G transition
            fprintf([num2str(prob_FG), '\n'])

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
            prob_AG = const_rates(site) / eventRatePerSite(site); % Probability of A -> G transition

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

    % Update varphi_rates matrix 
    for i = 1:num_sites

        % Clear row at index site in rate matrix
        varphi_rates(site, :) = zeros(1, num_sites);

        % Recalculate varphi rates in row
        if (X(site, k + 1) == 0) && (X(i, k + 1) == 1)
            varphi_rates(site, i) = varphi_G * W_P(locations(site, 1), locations(site, 2), locations(i, 1), locations(i, 2))/num_sites;
        end

        if (X(site, k + 1) == 3) && (X(site, k + 1) == 1)
            varphi_rates(site, i) = varphi_A * W_P(locations(site, 1), locations(site, 2), locations(i, 1), locations(i, 2))/num_sites;
        end

        % Clear column at index site in rate matrix
        varphi_rates(:, site) = zeros(num_sites, 1);

        % Recalculate varphi rates in column
        if (X(i, k + 1) == 0) && (X(site, k + 1) == 1)
            varphi_rates(i, site) = varphi_G * W_P(locations(site, 1), locations(site, 2), locations(i, 1), locations(i, 2))/num_sites;
        end

        if (X(i, k + 1) == 3) && (X(site, k + 1) == 1)
            varphi_rates(i, site) = varphi_A * W_P(locations(site, 1), locations(site, 2), locations(i, 1), locations(i, 2))/num_sites;
        end
    end

    % Update const_rates
    if (X(site, k + 1) == 0) % New state is G
        const_rates(site) = 0;
    end

    if (X(site, k + 1) == 1) % New state is F
        const_rates(site) = mu;
    end

    if (X(site, k + 1) == 2) % New state is B
        const_rates(site) = q;
    end

    if (X(site, k + 1) == 3) % New state is A
        const_rates(site) = gamma;
    end

end
toc
%% 
% Create array of full solution
Temp_Sol = zeros(num_sites, length(0:dt:T));
Sol_Save = zeros(length(0:dx:L), length(0:dy:L), length(0:dt:T));

% Interpolate temporally
for i = 1:num_sites
    Temp_Sol(i, :)= floor(interp1(Times, X(i,:), 0:dt:Times(end)));
end

% Interpolate spatially
for i = 1:length(0:dt:T)
    F = scatteredInterpolant(locations(:,1), locations(:,2), Temp_Sol(:, i));
    F.Method = 'nearest';

    for j = 1:length(0:dx:L)
        for k = 1:length(0:dy:L)
            Sol_Save(j, k, i) = floor(F(j*dx, k*dy));
        end
    end
end

% Make directory to save images
folderName = 'FGBA_2D_test';
mkdir(folderName);

% Save images to folder
for i = 1:length(0:dt:T)

    % Plot solution interpolated to a regularly spaced grid
    figure(1);
    im = imagesc(Sol_Save(:,:,i));

    % Set up color map
    custom_map = [0.4660 0.6740 0.1880
        0.0039 0.1953 0.1250
        0.8500 0.3250 0.0980
        0.5 0.5 0.5]; %light green, dark green, orange, gray colors
    colormap(custom_map);
    clim manual
    clim([0 3]);

    % Plot labels
    xlabel('x');
    ylabel('y');
    set(gca,'linewidth',2);
    set(gca,'FontSize',20);
    
    % Save image to directory
    imdata = getframe(figure(1));
    imwrite(imdata.cdata, [folderName '/t=' num2str(dt*i, '%2.2f') '.png']);
end





