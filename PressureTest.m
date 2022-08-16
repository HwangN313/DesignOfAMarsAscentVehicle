%% Pressure Test Script
close all 
clear all

Gamma = 1.1298;
Gamma_inf = 1.121;
C = 1.6;
sigma = (0.35*Gamma_inf) - 0.32;
M_inf = 3.719;
V_inf = 784.3;
rho_inf = 0.008596;

P_inj = 4.36e3;
V = 0.1196;
MON30Mass = 165.2;
molarMass = 2082.677;
n = MON30Mass/molarMass;
R = 146.1;
T = 1752.7;
rho_t = MON30Mass/V;




P_t = [1e5:1e5:100e5];
rho_inj = rho_t * ((P_inj./P_t).^(1/Gamma));
T_t = P_t .* V ./ (n*R);

M_inj = sqrt((((P_t./P_inj).^((Gamma-1)/Gamma)) - 1)/((Gamma-1)/2));
a = sqrt(Gamma*R*T);
Vel = M_inj*a;



ISP_i = C * sigma * M_inf * V_inf.*( 1 + (1/((Gamma_inf - 1) * (M_inf^2)))*(rho_inf./rho_inj))/9.81;
ISP_m = Vel./9.81;
ISP_s = 9.81*(ISP_i + ISP_m);


plot(P_t,ISP_s)
grid on
xlabel('Tank Pressure (Pa)')
ylabel('ISP (m/s)')
