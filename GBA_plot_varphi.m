area(data(:,1),data(:,2:5))
newcolors = [0.8500 0.3250 0.0980
    0.0039 0.1953 0.1250
    0.4660 0.6740 0.1880
    0.5 0.5 0.5];
colororder(newcolors)
xlabel('$\varphi$','interpreter','latex')
xlim([2,20]);
ylim([0,1]);
ylabel('Land fraction')