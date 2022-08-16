%% CEA Plots
clear all
close all

% 12
load('CEA020622V1.mat')
%10
load('CEA020622V2.mat')
%7.5
load('CEA020622V3.mat')
load('CEA020622V47barhighAEAT.mat')
load('CEA020622V57barhighAEAT.mat')

%% 30 Bar
% 30 bar, AeAt Variation

load('CEA310522V3.mat')


%% 15 Bar

% 15bar, AeAt Variation

load('CEA310522V4.mat')


%% 1 Bar
% 1 bar, AeAt Variation

load('CEA310522V2.mat')

%% Comb

figure
plot(CEA310522V2.re,CEA310522V2.Pbar,'-x','LineWidth',2) %1 
hold on
plot(CEA020622V3.re,CEA020622V3.P,'-x','LineWidth',2) %7.5
plot(CEA020622V2.re,CEA020622V2.P,'-x','LineWidth',2) %10
plot(CEA020622V1.re,CEA020622V1.P,'-x','LineWidth',2)%12
plot(CEA310522V4.re,CEA310522V4.P,'-x','LineWidth',2) %15
plot(CEA310522V3.re,CEA310522V3.P,'-x','LineWidth',2)%30
plot(CEA020622V47barhighAEAT.re,CEA020622V47barhighAEAT.P,'-x','LineWidth',2) 
grid on
xlabel('Exit radius (m)')
ylabel('Exit Pressure (bar)')
title('Exit Pressure Comparison')
yline(3.36e-3,'--')
xline(0.285)
%legend('1 bar','7.5 bar','10 bar','12 bar','15 bar','30 bar','Surface Pressure','Max Radius')
legend('1 bar','7.5 bar','10 bar','12 bar','15 bar','30 bar','7 bar high area ratio','Surface Pressure','Max Radius')


figure
plot(CEA310522V2.AeAt,CEA310522V2.rt,'--x','LineWidth',2,'Color','#0072BD')  %1
hold on
plot(CEA310522V2.AeAt,CEA310522V2.re,'-x','LineWidth',2,'Color','#0072BD')
plot(CEA020622V3.AeAt,CEA020622V3.rt,'--x','LineWidth',2,'Color','#D95319') %7.5
plot(CEA020622V3.AeAt,CEA020622V3.re,'-x','LineWidth',2,'Color','#D95319')
plot(CEA020622V2.AeAt,CEA020622V2.rt,'--x','LineWidth',2,'Color','#EDB120')%10
plot(CEA020622V2.AeAt,CEA020622V2.re,'-x','LineWidth',2,'Color','#EDB120')
plot(CEA020622V2.AeAt,CEA020622V1.rt,'--x','LineWidth',2,'Color','#7E2F8E')%12
plot(CEA020622V2.AeAt,CEA020622V1.re,'-x','LineWidth',2,'Color','#7E2F8E')
plot(CEA310522V4.AeAt,CEA310522V4.rt,'--x','LineWidth',2,'Color','#77AC30') %15
plot(CEA310522V4.AeAt,CEA310522V4.re,'-x','LineWidth',2,'Color','#77AC30')
plot(CEA310522V3.AeAt,CEA310522V3.rt,'--x','LineWidth',2,'Color','#4DBEEE') %30
plot(CEA310522V3.AeAt,CEA310522V3.re,'-x','LineWidth',2,'Color','#4DBEEE')
%plot(CEA020622V47barhighAEAT.AeAt,CEA020622V47barhighAEAT.rt,'-x','LineWidth',2,'Color','#A2142F')
%plot(CEA020622V47barhighAEAT.AeAt,CEA020622V47barhighAEAT.re,'-x','LineWidth',2,'Color','#A2142F')
yline(0.285);
grid on
title('Nozzle Radii')
xlabel('Area ratio')
ylabel('Radii (m)')
legend('1 bar - throat', '1 bar - exit','7.5 bar - throat', '7.5 bar - exit', '10 bar - throat','10 bar - exit','12 bar - throat', '12 bar - exit', '15 bar - throat', '15 bar - exit', '30 bar - throat' ,'30 bar - exit','Radius Limit')
%legend('1 bar - throat', '1 bar - exit','7.5 bar - throat', '7.5 bar - exit', '10 bar - throat','10 bar - exit','12 bar - throat', '12 bar - exit', '15 bar - throat', '15 bar - exit', '30 bar - throat' ,'30 bar - exit','7bar high throat','7bar high exit','Radius Limit')



figure

plot(CEA020622V57barhighAEAT.re,CEA020622V57barhighAEAT.P,'-x','LineWidth',2) 
hold on
grid on
xlabel('Exit radius (m)')
ylabel('Exit Pressure (bar)')
title('Exit Pressure Comparison')
yline(3.36e-3,'--')
xline(0.285)
legend('7.5 bar','Surface Pressure','Max Radius')
%% Ideal expansion



