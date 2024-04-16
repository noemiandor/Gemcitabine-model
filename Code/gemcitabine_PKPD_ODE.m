function dydt = gemcitabine_PKPD_ODE(t, y, theta, nu, eta, xi)

dydt = zeros(length(y), 1);

dydt(1) = -theta * y(1)/(y(1)  + nu); %@TODO: *N

dydt(2) = eta * theta * y(1)/(y(1)  + nu) - xi * y(2);
% disp(dydt)
end