function F = tauForceEq(tau,r,x)
    F(1) = (tau(1) - (-r(3)*x(2) + r(2)*x(3)))*(10^6);
    F(2) = tau(2) - (r(3)*x(1) - r(1)*x(3))*(10^6);
    F(3) = tau(3) - (-r(2)*x(1) + r(1)*x(2))*(10^6);
end
