S_db = 10;
S = 10^(S_db/10);
I_db = 15;
I = 10^(I_db/10);

M = 16;
i = (1:M);
[re, im] = meshgrid([-3 -1 1 3], [-3 -1 1 3]);
x = re(:) + 1j*im(:);
x = x / sqrt(mean(abs(x.^2)));

[X,Y] = meshgrid(linspace(-11,11,1500), linspace(-11,11,1500));
y = X + 1j*Y;

decision_region = zeros([size(y), length(x)]);

for k = 1:M
    d_y = abs(y-sqrt(S)*x(k));
    decision_region(:,:,k) = d_y.^2 - log(besseli(0, 2*sqrt(I).*d_y));
end

[~, decision_map] = min(decision_region, [], 3);

figure;
imagesc([-11 11], [-11 11], decision_map);
axis xy;
xlabel('Re(Y)'); ylabel('Im(Y)');
title('IdB=15dB, ML Optimal decoder');
colorbar;
hold on;
plot(real(sqrt(S)*x), imag(sqrt(S)*x), 'ko', 'MarkerEdgeColor', 'none', 'MarkerFaceColor','r', 'MarkerSize', 8);
for k = 1:M
    text(real(sqrt(S)*x(k)) + 0.3, imag(sqrt(S)*x(k)), num2str(k), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
end
