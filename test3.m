% Plot locations of sites
time_index = 8;


figure(1)
for i = 1:num_sites
    if (X(i, time_index) == 0)
        plot(locations(i, 1), locations(i, 2),'o', 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
        hold on
    end

    if (X(i, time_index) == 1)
        plot(locations(i, 1), locations(i, 2),'o', 'MarkerFaceColor', [0.0039 0.1953 0.1250], 'MarkerEdgeColor',[0.0039 0.1953 0.1250]);
        hold on
    end

    if (X(i, time_index) == 2)
        plot(locations(i, 1), locations(i, 2),'o', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor',[0.8500 0.3250 0.0980]);
        hold on
    end

    if (X(i, time_index) == 3)
        plot(locations(i, 1), locations(i, 2),'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5]);
        hold on
    end
end

% Plot solution interpolated to a regularly spaced grid
%image(Test, 'CDataMapping','scaled')
%image(Sol_Save(:, :, 1));

%im = imagesc(Sol_Save(:, :, 2));
custom_map = [0.4660 0.6740 0.1880
    0.8500 0.3250 0.0980
    0.5 0.5 0.5
    0.0039 0.1953 0.1250]; %light green, orange, gray, dark green colors
colormap(custom_map);
xlabel('time');
ylabel('space');
%set(gca,'linewidth',2);
%set(gca,'FontSize',20);

imdata = getframe(figure(1));

imwrite(imdata.cdata, 'im_test.png')