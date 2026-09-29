S_db = 10;
S = 10^(S_db/10);
I_db = 15;
I = 10^(I_db/10);

M = 8;
i = (1:M);
theta = 2*pi*(i-1)/M;
x = exp(1j*theta);

[X,Y] = meshgrid(linspace(-10,10,1500), linspace(-10,10,1500));

decision_region = zeros(1500, 1500);
decision_map = zeros(1500,1500);

for ix = 1:1500
    for iy = 1:1500
        y = X(ix,iy) + 1j*Y(ix,iy);
        best = 1;
        for l = 1:M
            d_yl = abs(y-sqrt(S)*x(l));
            
            better = true;
            for k = 1:M
                if k == l
                    continue;
                end
                d_yk = abs(y-sqrt(S)*x(k));

                if d_yl > sqrt(I) && d_yk > sqrt(I)
                    better = better && d_yl<d_yk;
                elseif d_yl > sqrt(I) && d_yk <= sqrt(I)
                    better = better && (d_yl+d_yk)/2<sqrt(I);
                elseif d_yl <= sqrt(I) && d_yk > sqrt(I)
                    better = better && (d_yl+d_yk)/2>=sqrt(I);
                elseif d_yl <= sqrt(I) && d_yk <= sqrt(I)
                    better = better && d_yl>=d_yk;
                end
            end
        if better
            best = l;
            break;
        end
        end
        decision_map(ix,iy) = best;
    end
end

figure;
imagesc([-10 10], [-10 10], decision_map);
axis xy;
xlabel('Re(Y)'); ylabel('Im(Y)');
title('IdB=15dB, IC Sub-Optimal decoder');
colorbar;
hold on;
plot(real(sqrt(S)*x), imag(sqrt(S)*x), 'ko', 'MarkerEdgeColor', 'none', 'MarkerFaceColor','r', 'MarkerSize', 8);
%%for k = 1:M
%%    text(real(sqrt(S)*x(k)) + 0.3, imag(sqrt(S)*x(k)), num2str(k), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
%%end
