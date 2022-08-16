% Setting ISP_s at 860 + 10% for shock
%ISP_s = 946; %m/s

Cons = 1.6;
Gamma_Inf = 1.1298;
Omega = (0.35*Gamma_Inf) - 0.32;
M_Inf = 3.719;
V_Inf = 784.3;
rho_Inf = 8.596e-3;

rho_t = 165.2 / 0.1196;


Arat = [1.1:0.1:5];
A_star = 1.64e-5;



A = Cons * Omega * M_Inf * V_Inf;
B = (Gamma_Inf - 1)* (M_Inf^2);
C = rho_Inf/rho_i;
% C = 0.75;


ISP_i = A * (1 + ((1/B)*C))/9.81;


% ISP_m = ISP_s - ISP_i; % s
% Vi = ISP_m*9.81; %m/s




