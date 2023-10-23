dx = 0.025;
dy = 0.025;

% Calculate indices of positions
location_indices(:, 1) = ceil(locations(:, 1)/dx);
location_indices(:, 2) = ceil(locations(:, 2)/dy);

% Create array of full solution
Sol_Save = zeros(ceil(L/dx), ceil(L/dy), length(0:dt:T));

% Plot locations of sites
time_index = 1;


figure(1)
Test = zeros(2,2);
Test(1,1) = 0;
Test(1,2) = 1;
Test(2,2) = 2;
Test(2,1) = 3;
custom_map = [0.4660 0.6740 0.1880
    0.0039 0.1953 0.1250
    0.8500 0.3250 0.0980
    0.5 0.5 0.5];
colormap(custom_map);

% Plot solution interpolated to a regularly spaced grid
image(Test, 'CDataMapping','scaled')

figure(2)
imagesc(Sol_Save(:, :, 1));

%im = imagesc(Sol_Save(:, :, 2));
custom_map = [0.4660 0.6740 0.1880
    0.0039 0.1953 0.1250
    0.8500 0.3250 0.0980
    0.5 0.5 0.5]; %light green, dark green, orange, gray colors
colormap(custom_map);
xlabel('time');
ylabel('space');
%set(gca,'linewidth',2);
%set(gca,'FontSize',20);

imdata = getframe(figure(1));

imwrite(imdata.cdata, 'im_test.png')