%% PE Vex
close all
clear all

ISP = 3068.568; %m/s
m_dot = 1.703;  %kg/s

% ISP = V_ex + (Pe - Pa)*Ae/m_dot


Ae = (0.53/2)^2*pi;     %57cm diameter - 2cm from each side for nozzle thickness

Pa = 630; %Pa at surface

Pe = [0:0.1e6:5e6];

V_ex = ISP - ((Pe - Pa).*(Ae/m_dot));

plot(Pe, V_ex)
xlabel('Exhaust Pressure (Pa)')
ylabel('Exhaust Velocity (m/s)')

%Thrust = m_dot*V_ex + (Pe - Pa)*Ae


