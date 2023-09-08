tiledlayout(2,1)

% Top plot
nexttile
plot(HD_DISK_O2(:,1), log(HD_DISK_O2(:,2)),'-o')
title('HD Trial: Average O2 vs Time')
xlabel('Hours after plugging')
ylabel('log(fraction O2)')

% Bottom plot
nexttile
plot(LD_DISK_O2(:,1), log(LD_DISK_O2(:,2)),'-o')
title('LD Trial: Average O2 vs Time')
xlabel('Hours after plugging')
ylabel('log(fraction O2)')