function [AreaRatio] = ARat(M1,Gm)
%Calc for solving Area Ratio given Gamma and Mach number
% AreaRatio = (1/M1)*((1+((Gm-1)/2)*(M1^2))/((Gm+1)/2))^((Gm+1)/(2*(Gm-1)));


T1 = (Gm+1)/2;
T2 = (Gm+1)/(2*(Gm-1));
T3 = (1 + ((Gm-1)/2).*M1.^2);


AreaRatio = (T1^-T2).*((T3.^T2)./M1);
end