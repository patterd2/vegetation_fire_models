%% Matrix simulation of 2D FGBA system with periodic boundary conditions
clear
close all
set(0,'defaultaxesfontsize', 20);
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');
%%
% Convention for state labels: 0 = Grass, 1 = Forest, 2 = Burning, 3 = Ash

% Parameter being varied (REMOVE FROM LOOP)
param_nums = [0.01, 0.1, 1, 2];
param_name = 'sigma_B';

folderName = ['vary6_' param_name];
mkdir(folderName);
fig_num = 1;

for n = 1:size(param_nums, 2)
    tic
    
    param_value = param_nums(n);
    fprintf(['\n' param_name '=' num2str(param_value) '\n'])

    % Parameters
    L = 5;            % Working on [0, L]
    num_sites = 5000; % Total number of sites in [0, L] x [0, L]
    T = 10;            % Simulation time length in real-time units
    t0 = 0;           % Start time
    
    varphi_A = 0.1;  % Rate of forest seeding into ash
    varphi_G = 0.1;  % Rate of forest seeding into grass
    
    beta_F = 10;    % Rate of fire spread through forest
    beta_G = 200;    % Rate of fire spread through grass
    gamma = 10;     % Rate of grass regrowth from ash
    q = 10;         % Fire quenching rate
    
    mu = 0.01;       % Rate of tree death due to non-fire causes
    
    sigma_P = 0.1;   % width of Gaussian of forest seeding
    sigma_B = param_value;    % width of Gaussian for burning spread over grassland
    sigma_G = 0.1;    % width of Gaussian for flammability of large grass
    sigma_F = 0.1;    % width of Gaussian for flammability of large grass
    
    theta_F = 0.4;    % Forest burning sigmoid center
    f0 = 0.01;        % Forest burning sigmoid lower bound
    f1 = 0.5;         % Forest burning sigmoid upper bound
    s_f = 0.05;       %orest burning sigmoid width
    
    theta_G = 0.4;    % Grass burning sigmoid center
    g0 = 0.01;        % Grass burning sigmoid lower bound
    g1 = 1;         % Grass burning sigmoid upper bound
    s_g = 0.05;       % Grass burning sigmoid width

    % Save all parameters first time
    if (n == 1)
        save([folderName '/parameters.mat'], 'L', 'num_sites', 'T', 't0', ...
            'varphi_G', 'varphi_A', 'beta_G', 'beta_F', 'gamma', 'q', 'mu', ...
            'sigma_P', 'sigma_B', 'sigma_G' ,'sigma_F', 'theta_F', 'f0', 'f1', ...
            's_f', 'theta_G', 'g0', 'g1', 's_g')
    end
    
    % Periodic distance function
    dist_1D = @(r1, r2) angle (exp(1i * (r1 - r2) * 2 * pi / L)) * L / (2 * pi); 
    dist_2D = @(x1, y1, x2, y2) sqrt(dist_1D(x1, x2)^2 + dist_1D(y1, y2)^2);
    
    % Gaussian functions
    W_B = @(x1, y1, x2, y2) (1 / num_sites) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_B^2)) / sqrt(2 * pi * sigma_B^2); % For cascade effect
    W_G = @(x1, y1, x2, y2) (1 / num_sites) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_G^2)) / sqrt(2 * pi * sigma_G^2); % For beta spreading
    W_F = @(x1, y1, x2, y2) (1 / num_sites) * exp( -( dist_2D(x1, y1, x2, y2).^2 ) / (2 * sigma_F^2)) / sqrt(2 * pi * sigma_F^2); % For varphi spreading
    
    % Sigmoid forest mortality function
    phi_F = @(x) (f0 + (f1 - f0)./(1 + exp(-(x - theta_F) / s_f)));
    phi_G = @(x) (g0 + (g1 - g0)./(1 + exp(-(x - theta_G) / s_g)));
    
    % Allocate states of sites
    X = zeros(num_sites, 1);
    
    % Fix locations of sites (x, y)
    locations = L * rand(num_sites, 2);  % NOT SORTED
    
    % Fix square root of number of sites
    sqrt_num_sites = sqrt(num_sites);
    
    % Initialize states of sites - one square forest patch at [0.2L, 0.5L] x [0.2L, 0.5L]
    %for i = 1:num_sites
        %if (locations(i, 1) < 0.5 * L && locations(i, 1) > 0.2 * L)
            %if (locations(i, 2) < 0.5 * L && locations(i, 2) > 0.2 * L)
                %X(i, 1) = 1;
            %end
        %end
    %end
    
    % Create state indicator vectors
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
    
    for i = 1:num_sites
        for j = 1:num_sites
            if (X(j, 1) == 0) % Integrate over G sites
                G_integral(i, j) = W_G(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end
    
            if (X(j, 1) == 2) % Integrate over B sites
                B_integral(i, j) = W_B(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end
    
            if (X(j, 1) == 1) % Integrate over F sites
                F_integral(i, j) = W_F(locations(i, 1), locations(i, 2), locations(j, 1), locations(j, 2));
            end
        end
    end
    
    % Use Gillespie algorithm to model transitions over time
    t = 0; % time
    k = 0; % time index
    Times = [0]; % array to record times of all events
    
    tic;
    while (t < T)
        k = k + 1;
    
        % Calculate event rate at each site
        eventRatePerSite = is_G .* (varphi_G * sum(F_integral, 2) + phi_G(sum(G_integral, 2)) + beta_G * sum(B_integral, 2)) ...
            + is_F .* (mu + phi_F(sum(G_integral, 2)) + beta_F * sum(B_integral, 2)) + ...
            + is_B .* q + ...
            + is_A .* (gamma + varphi_A * sum(F_integral, 2));
    
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
    
    end
    toc % time the calculations of the Gillespie algorithm
    %% Different approach to plotting space-time dynamics
    figure(1);
    % Set up color map
    custom_map = [0.4660 0.6740 0.1880
        0.0039 0.1953 0.1250
        0.8500 0.3250 0.0980
        0.5 0.5 0.5]; %light green, dark green, orange, gray colors
    colormap(custom_map);
    grid on;
    %for i = 1:length(Times)
        %if mod(i,100) == 0
            %str = ['Time = ',num2str(Times(i))];
            %figure(1);
            %title(str);
            %scatter(locations(:,1),locations(:,2),150,X(:,i),...
                %'filled','MarkerFaceAlpha',0.2);
            %hold on;
        %end
    %end

    %% Plot the dynamics of the cover proportions
    % 0 = Grass, 1 = Forest, 2 = Burning, 3 = Ash
    figure(fig_num);
    stairs(Times,sum(X==0,1)./num_sites,'color',[0.4660 0.6740 0.1880],'LineWidth',3);
    hold on;
    stairs(Times,sum(X==1,1)./num_sites,'color',[0.0039 0.1953 0.1250],'LineWidth',3);
    stairs(Times,sum(X==2,1)./num_sites,'color',[0.8500 0.3250 0.0980],'LineWidth',3);
    stairs(Times,sum(X==3,1)./num_sites,'color',[0.5 0.5 0.5],'LineWidth',3);
    legend('Grass','Forest','Burning','Ash');
    xlabel('time');

    % Save figure
    saveas(gcf, [folderName '/' param_name '=' num2str(param_value) '.jpeg']);

    % Save X variables
    save([folderName '/X_' param_name '=' num2str(param_value) '.mat'], 'X','-v7.3')

    fig_num = fig_num + 1;
end

%%
% tic % time the plotting
% % Create array of full solution
% Temp_Sol = zeros(num_sites, length(0:dt:T));
% Sol_Save = zeros(length(0:dx:L), length(0:dy:L), length(0:dt:T));
%
% % Interpolate temporally
% for i = 1:num_sites
%     Temp_Sol(i, :)= floor(interp1(Times, X(i,:), 0:dt:Times(end)));
% end
% toc
% %%
% tic
% % Interpolate spatially
% for i = 1:length(0:dt:T)
%     F = scatteredInterpolant(locations(:,1), locations(:,2), Temp_Sol(:, i));
%     F.Method = 'nearest';
%     for j = 1:length(0:dx:L)
%         for k = 1:length(0:dy:L)
%             Sol_Save(j, k, i) = floor(F(j*dx, k*dy));
%         end
%     end
% end
% toc
% % Make directory to save images
% %folderName = 'FGBA_2D_test';
% %mkdir(folderName);
%
% %%
% tic
% % Save images to folder
% for i = 1:length(0:dt:T)
%
%     % Plot solution interpolated to a regularly spaced grid
%     figure(1);
%     im = imagesc(Sol_Save(:,:,i));
%
%     % Set up color map
%     custom_map = [0.4660 0.6740 0.1880
%         0.0039 0.1953 0.1250
%         0.8500 0.3250 0.0980
%         0.5 0.5 0.5]; %light green, dark green, orange, gray colors
%     colormap(custom_map);
%     clim manual
%     clim([0 3]);
%
%     % Plot labels
%     xlabel('x');
%     ylabel('y');
%     set(gca,'linewidth',2);
%     set(gca,'FontSize',20);
%
%     % Save image to directory
%     %imdata = getframe(figure(1));
%     %imwrite(imdata.cdata, [folderName '/t=' num2str(dt*i, '%2.2f') '.png']);
% end
%
% toc



