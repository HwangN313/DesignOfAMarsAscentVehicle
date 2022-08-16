function [sol] = isentropicFlow(Gamma, M)
% 

for i = 1:length(Gamma) % 
    G = Gamma(i);
    Gp1 = G + 1;
    Gm1 = G - 1;
    a = G/(Gm1);% 
    b = 1/(Gm1);% 
    t0_t(i,:) = 1 + (Gamma(i) - 1)/2*M.^2;
    p0_p(i,:) = (t0_t(i,:)).^a;
    rho0_rho(i,:) = (t0_t(i,:)).^b;
    a_aStar(i,:) = (Gp1./2).^(-Gp1./(2.*Gm1))*(1 + Gm1/2.*M.^2).^(Gp1./((2*Gm1)))./M;
end

sol.M  = M;
sol.t0t = t0_t;
sol.p0p = p0_p;
sol.r0r = rho0_rho;
sol.aaS = a_aStar;

sol.t0t_i = 1./t0_t;
sol.p0p_i = 1./p0_p;
sol.r0r_i = 1./rho0_rho;
sol.aaS_i = 1./a_aStar;