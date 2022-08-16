%% Sizing Loop
% Loop outputting fuel and inert masses

clear all
close all

ISP         = 3068.568; %m/s
Target_DV   = 4570.045; %m/s

m0          = 400;      %kg
mpay        = 16;       %kg
m_prop_struc = 7;       %kg


F1          = 0.78;
F2          = 1-F1;

mfuel_init  = 297.5;            %kg
minert_init = m0 - mpay - (2*m_prop_struc) - mfuel_init; %kg


i = 1;
for i = 1:2000
    
    if i ==1
        mfuel = mfuel_init;
        minert = minert_init;
    else
        mfuel = mfuel - 0.5;
        minert = minert + 0.5;
    end
    
    mfuel_s1 = mfuel*F1;
    mfuel_s2 = mfuel*F2;

    minert_s1 = minert*F1;
    minert_s2 = minert*F2;

    m1 = m0 - mfuel_s1;
    DV1 = ISP*log(m0/m1);

    m1_2 = m1 - minert_s1 - m_prop_struc;
    m2 = m1_2 - mfuel_s2;
    
    DV2 = ISP*log(m1_2/m2);

    DV = DV1 + DV2
    if DV<Target_DV
        break
    end
    
    Prop_Frac = mfuel/m0;
    
    Paraffin_Mass = mfuel/5.5;
    Mon30_Mass = mfuel*4.5/5.5;
    
    Sol_Mat(i,1) = mfuel;
    Sol_Mat(i,2) = mfuel_s1;
    Sol_Mat(i,3) = mfuel_s2;
    Sol_Mat(i,4) = minert + (2*m_prop_struc);
    Sol_Mat(i,5) = DV1;
    Sol_Mat(i,6) = DV2;
    Sol_Mat(i,7) = DV;
    Sol_Mat(i,8) = Prop_Frac;
    Sol_Mat(i,9) = Paraffin_Mass;
    Sol_Mat(i,10) = Mon30_Mass;
    Sol_Mat(i,11) = minert_s1+m_prop_struc;
    Sol_Mat(i,12) = m2 - mpay;
    
end


Final_Val_Tab = array2table(Sol_Mat(end,:),'VariableNames',{'Fuel Mass','FuelS1','FuelS2','Total Inert Mass','DV1','DV2','Delta-V','Propellant Fraction','Paraffin Mass','MON30 Mass','1st Stage Inert Mass','Final Inert Mass'});









