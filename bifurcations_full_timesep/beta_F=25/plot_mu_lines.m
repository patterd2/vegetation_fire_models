%% Plot mu lines
line_1 = false;
X = linspace(0, 1);

if (line_1)
    x_intercept = 0.609512;
    slope = -2.72533;
    Y = slope * (X - x_intercept);
else
    y_intercept = 0.989987;
    slope = -1.00003
    Y = slope * X + y_intercept;
end

hold on
plot(X,Y)