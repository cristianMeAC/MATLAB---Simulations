function F = solveLag(x)

wc = 500;
dphi = -6.796533154601193;
da = 0.481276870689275;

F(1) = atand((wc*x(1) - x(2)*wc*x(1)) / (1 + x(2)*wc^2*x(1)^2)) - dphi; % S.129
F(2) = (sqrt(1 + wc^2*x(1)^2) / sqrt(1 + x(2)^2*wc^2*x(1)^2)) - da;     % S.129