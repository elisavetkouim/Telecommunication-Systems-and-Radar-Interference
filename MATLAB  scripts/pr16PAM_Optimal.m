S_db = 10;
S = 10^(S_db/10);
I_db = 15;
I = 10^(I_db/10);

M = 16;
i = (1:M);
x = (2*i-M-1);
x = x / sqrt(mean(x.^2));

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
text(real(sqrt(S)*x(1)) - 1, imag(sqrt(S)*x(1)), num2str(1), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
text(real(sqrt(S)*x(16)) + 0.3, imag(sqrt(S)*x(16)), num2str(16), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
