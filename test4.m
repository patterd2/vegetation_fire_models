%interp2(reshape(locations(:,1),[1,num_sites]), ...
        %reshape(locations(:,2),[1,num_sites]), Temp_Sol(:, i), 0:dx:L, 0:dy:L)

test = reshape(locations(:,1),[1,num_sites]);

X = zeros(2);
Y = zeros(2);

for i = 1:length(0:dt:T)
    Sol_Save(:, :, i) = scatteredinterpolant(X(i), Y(i), Temp_Sol(:, i), 0:dx:L, 0:dy:L);
end