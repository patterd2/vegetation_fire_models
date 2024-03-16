beta_G_vals = [2e4, 4e4, 6e4, 8e4, 1e5, 1.2e5, 1.4e5, 1.6e5, 1.8e5];

% Plot surface
[X,Y] = meshgrid(0.1:0.01:0.9, beta_G_vals(1):10:beta_G_vals(end));
Z = 1 ./ (X .* Y);
s = surf(X,Y,Z,'EdgeColor', 'none');
xlabel('Proportion ash', 'Interpreter','latex')
ylabel('$$\beta_G$$ (yr$^{-1}$)' , 'Interpreter', 'latex')
zlabel('Burning start time (yr)', 'Interpreter','latex')
