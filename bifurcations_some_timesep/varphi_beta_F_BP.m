% Plot estimations of the varphi-beta_F BP line

%varphi = 5;
mu_tilde = 0.11;
%beta_F = 10;
beta_G = 50;
gamma = 5;
q = 20;
g1 = 0.05;

f0 = @(varphi, beta_F) (varphi - mu_tilde) / (varphi + beta_F) ...
    - (gamma * (beta_G - q))/(beta_G * (q + gamma));
f1 = @(varphi, beta_F) (varphi - mu_tilde) / (varphi + beta_F) ...
    - (gamma * (beta_G - q))/(beta_G * (q + gamma))...
    - (q * (q + gamma) * g1^2)/(beta_G * (beta_G - q));

figure(1)
fimplicit(f0, [0, 6, 0, 20])

hold on
fimplicit(f1, [0, 6, 0, 20])
plot(data(:, 1), data(:,2))