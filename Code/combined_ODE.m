function dydt_combined = combined_ODE(t, y, theta, nu_PKPD, eta, xi, a, v, w1, iota, nu)

% Define dydt for the combined system
dydt_combined = zeros(length(y), 1);

% Define the first two equations from gemcitabine_PKPD_ODE
dydt_combined(1) = -theta * y(1) / (y(1) + nu_PKPD); %@TODO We need to introduce the total population in this part
dydt_combined(2) = eta * theta * y(1) / (y(1) + nu_PKPD) - xi * y(2);

% Use dydt_combined(2) as DOSE for skippedMito_ODE
DOSE = dydt_combined(2);

% Get the rest of the equations from skippedMito_ODE
for i = 3:length(y)
    mitSkipped = i - 3;
    dydt_combined(i) = alpha_p(mitSkipped, DOSE) * y(i) - k_p(mitSkipped, DOSE) * y(i) - a_p(mitSkipped, DOSE) * y(i);
    if mitSkipped >= 1
        dydt_combined(i) = dydt_combined(i) + k_p(mitSkipped - 1, DOSE) * y(i - 1);
    end
    %% sum up dead (apoptotic) cells
    dydt_combined(1) = dydt_combined(1) + a_p(mitSkipped, DOSE) * y(i);
end

% Define the functions for a, k_p and a_p
    function kp = k_p(i, x)
        kp = (nu - (1 - i)^2) * (x / (x + v));
    end

    function alphap = alpha_p(i, x)
        alphap = x / (x + a) * (i <= iota);
    end

    function ap = a_p(i, x)
        ap = (((2 - i) / 2))^2 * (x / (x + w1));
    end

end
