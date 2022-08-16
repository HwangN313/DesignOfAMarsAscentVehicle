function [x3,x1,f1] = Secant(x1,x2,f1,f2)
%Secant Method - Guesses to solve a function
x3 = x2 - ((x2 - x1)/(f2 - f1))*f2;
x1 = x2;
f1 = f2;
end

