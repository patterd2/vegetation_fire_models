%% Measure time for ash patch surrounded by forest/grass to grow tree
clear
close all
set(0,'defaultaxesfontsize', 20);
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');

% Parameters
L = 100;            % Working on [0, L]
num_sites = 10000; % Total number of sites in [0, L] x [0, L]
T = 1000;            % Simulation time length in real-time units
dt = 0.001;        % Length of time step in real-time units
dx = 1;        % Spatial resolution in x direction
dy = 1;        % Spatial resolution in y direction
t0 = 0;           % Start time

% Parameter being varied (REMOVE FROM LOOP)

fig_num = 1;
prop_grass = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
varphi_A_vals = [2e-2, 4e-2, 6e-2, 8e-2, 1e-1, 1.2e-1, 1.4e-1, 1.6e-1, 1.8e-1];
num_trials = 5;
burn_times = zeros(size(prop_grass, 2), size(varphi_A_vals, 2), num_trials);


%% Begin loop
for p = 1:size(prop_grass, 2)
    for vA = 1:size(varphi_A_vals, 2)
        for n = 1:num_trials
            % Convention for state labels: 0 = Grass, 1 = Forest, 2 = Burning, 3 = Ash
            tic

            fprintf(['prop_ash=' num2str(prop_grass(p)) '  varphi_A=' num2str(varphi_A_vals(vA)) '  trial ' num2str(n)])
            
            % Parameters
            varphi_A = varphi_A_vals(vA);  % Rate of forest seeding into ash
            varphi_G = 0;  % Rate of forest seeding into grass
            
            beta_F = 10^4;    % Rate of fire spread through forest
            beta_G = 10^5;    % Rate of fire spread through grass
            gamma = 0;     % Rate of grass regrowth from ash
            %q = 10^3;         % Fire quenching rate
            q = 0;
            
            %mu = 0.01;       % Rate of tree death due to non-fire causes
            mu = 0;
            
            sigma_P = 10;   % width of Gaussian of forest seeding
            sigma_B = 10;    % width of Gaussian for burning spread over grassland
            sigma_G = 5;    % width of Gaussian for flammability of large grass
            sigma_F = 5;    % width of Gaussian for flammability of large grass
            
            theta_F = 0.4;    % Forest burning sigmoid center
            %f0 = 0.01;        % Forest burning sigmoid lower bound
            f0 = 0;
            %f1 = 0.1;         % Forest burning sigmoid upper bound
            f1 = 0;
            s_F = 0.05;       %orest burning sigmoid width
            
            theta_G = 0.4;    % Grass burning sigmoid center
            %g0 = 0.01;        % Grass burning sigmoid lower bound
            g0 = 0;
            %g1 = 1;         % Grass burning sigmoid upper bound
            g1 = 0;
            s_G = 0.05;       % Grass burning sigmoid width

            
            % Periodic distance function
            dist_1D = @(r1, r2) angle (exp(1i * (r1 - r2) * 2 * pi / L)) * L / (2 * pi); 
            dist_2D = @(x1, y1, x2, y2) sqrt(dist_1D(x1, x2)^2 + dist_1D(y1, y2)^2);
            
            % Gaussian functions
            W_B = @(x1, y1, x2, y2) (L^2 / (2 * pi * sigma_B^2 * num_sites)) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_B^2)); % For cascade effect
            W_G = @(x1, y1, x2, y2) (L^2 / (2 * pi * sigma_G^2 * num_sites)) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_G^2)); % For beta spreading
            W_F = @(x1, y1, x2, y2) (L^2 / (2 * pi * sigma_F^2 * num_sites)) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_F^2)); % For varphi spreading
            
            % Sigmoid forest mortality function
            phi_F = @(x) (f0 + (f1 - f0)./(1 + exp(-(x - theta_F) / s_F)));
            phi_G = @(x) (g0 + (g1 - g0)./(1 + exp(-(x - theta_G) / s_G)));
            
            % Allocate states of sites
            X = zeros(num_sites, 1);
            
            % Fix locations of sites (x, y)
            locations = L * rand(num_sites, 2);  % NOT SORTED
            
            %% Initialize states of sites (one ash site + mixture of forest/grass)
            ash_site = ceil(rand() * num_sites);

            for i = 1:num_sites
                if (i == ash_site)
                    X(i, 1) = 3;
                else
                    if (rand() > prop_grass(p))
                        X(i, 1) = 1;
                    else
                        X(i, 1) = 0;
                    end
                end
            end
            
            %% Create state indicator vectors
            is_G = zeros(num_sites, 1);
            is_B = zeros(num_sites, 1);
            is_F = zeros(num_sites, 1);
            is_A = zeros(num_sites, 1);
            
            for i = 1:num_sites
                if (X(i, 1) == 0)
                    is_G(i) = 1;
                elseif (X(i, 1) == 2)
                    is_B(i) = 1;
                elseif (X(i, 1) == 1)
                    is_F(i) = 1;
                else
                    is_A(i) = 1;
                end
            end
            
            % Gaussian integrals centered at each site i
            G_integral = zeros(num_sites, num_sites);
            B_integral = zeros(num_sites, num_sites);
            F_integral = zeros(num_sites, num_sites);
            
            f = waitbar(0, 'computing integrals');
            for i = 1:num_sites
                for j = 1:num_sites
                    if (X(j, 1) == 0) % Integrate over G sites
                        G_integral(i, j) = W_G(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
                    end
            
                    if (X(j, 1) == 2) % Integrate over B sites
                        B_integral(i, j) = W_B(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
                    end
                    %fprintf([num2str(sum(B_integral(1,:))) '\n'])
            
                    if (X(j, 1) == 1) % Integrate over F sites
                        F_integral(i, j) = W_F(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
                    end
                end
                waitbar(i/num_sites, f)
            end
            close(f)
            
            % Use Gillespie algorithm to model transitions over time
            t = 0; % time
            k = 0; % time index
            Times = [0]; % array to record times of all events
            
            tic;
            f = waitbar(0, 'running FGBA simulation');
        
            
        
            while (t < T)
                k = k + 1;
            
                % Calculate event rate at each site
                eventRatePerSite = is_G .* (varphi_G * sum(F_integral, 2) + phi_G(sum(G_integral, 2)) + beta_G * sum(B_integral, 2)) ...
                    + is_F .* ((mu + phi_F(sum(G_integral, 2)) + beta_F * sum(B_integral, 2))) + ...
                    + is_B .* q + ...
                    + is_A .* (gamma + varphi_A * sum(F_integral, 2));
            
                % Calculate total intensity over all sites
                totalIntensity = sum(eventRatePerSite);
            
                % Find time of next event
                NextEvent = -log(1 - rand()) / totalIntensity;
                t = t + NextEvent;
        
                time_recorded = false;
            
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
                        GB_rate = phi_G(sum(G_integral(site, :))) + beta_G * sum(B_integral(site, :));
                        prob_GB = GB_rate / eventRatePerSite(site); % Probability of G -> B transition
            
                        if (rand() < prob_GB)
                            X(site, k + 1) = 2; % G to B
            
                            % Adjust indicators
                            is_G(site) = 0;
                            is_B(site) = 1;
            
                            % Adjust spatial integrals
                            G_integral(:, site) = zeros(num_sites, 1);
                            for i = 1:num_sites
                                B_integral(i, site) = W_B(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
                        else
                            X(site, k + 1) = 1; % G to F
            
                            % Adjust indicators
                            is_G(site) = 0;
                            is_F(site) = 1;
            
                            % Adjust spatial integrals
                            G_integral(:, site) = zeros(num_sites, 1);
                            for i = 1:num_sites
                                F_integral(i, site) = W_F(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
                        end
                    end
            
                    % Choose transition if site is F
                    if (X(site, k) == 1)
                        prob_FG = mu / eventRatePerSite(site); % Probability of F -> G transition
            
                        if (rand() < prob_FG)
                            X(site, k + 1) = 0; % F to G
            
                            % Adjust indicators
                            is_F(site) = 0;
                            is_G(site) = 1;
            
                            % Adjust spatial integrals
                            F_integral(:, site) = zeros(num_sites, 1);
                            for i = 1:num_sites
                                G_integral(i, site) = W_G(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
            
                        else
                            X(site, k + 1) = 2; % F to B
            
                            % Adjust indicators
                            is_F(site) = 0;
                            is_B(site) = 1;
            
                            % Adjust spatial integrals
                            F_integral(:, site) = zeros(num_sites, 1);
                            for i = 1:num_sites
                                B_integral(i, site) = W_B(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
                        end
                    end
            
                    % Choose transition if site is B
                    if (X(site, k) == 2)
                        X(site, k + 1) = 3; % fire must get quenched
            
                        % Adjust indicators
                        is_B(site) = 0;
                        is_A(site) = 1;
            
                        % Adjust spatial integrals
                        B_integral(:, site) = zeros(num_sites, 1);
                    end
            
                    % Choose transition if site is A
                    if (X(site, k) == 3)
                        prob_AG = gamma / eventRatePerSite(site); % Probability of A -> G transition
            
                        if (rand() < prob_AG)
                            X(site, k + 1) = 0; % A to G
            
                            % Adjust indicators
                            is_A(site) = 0;
                            is_G(site) = 1;
            
                            % Adjust spatial integrals
                            for i = 1:num_sites
                                G_integral(i, site) = W_G(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
                        else
                            X(site, k + 1) = 1; % A to F

                            if (time_recorded == false)
                                burn_times(p, vA, n) = t;
                                time_recorded = true;
                            end
            
                            % Adjust indicators
                            is_A(site) = 0;
                            is_F(site) = 1;
            
                            % Adjust spatial integrals
                            for i = 1:num_sites
                                F_integral(i, site) = W_F(locations(i, 1), locations(i, 2), locations(site, 1), locations(site, 2));
                            end
                        end
                    end
            
                else % Finish simulation if event occurs after end time of simulation
                    Times(end + 1) = T;
                    X(:, k + 1) = X(:, k);
                end
            
                waitbar(Times(end)/T, f)
            
            end
            close(f)
            
            toc % time the calculations of the Gillespie algorithm
        end
    end
end