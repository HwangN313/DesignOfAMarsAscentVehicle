%% ISP check
close all
clear all

Steps = linspace(1,100);
Pe = 1e4;
m_dot = 1.703;
Pa = linspace(630,0);
vex = 1750;
Ae = (0.53/2)^2*pi; 

Thrust = m_dot*vex + (Pe - Pa).*Ae;
ISP = vex + (Pe - Pa).*Ae/m_dot;



plot(Pa,Thrust)
hold on
plot(Pa,ISP)