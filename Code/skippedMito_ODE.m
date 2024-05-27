

function dydt = skippedMito_ODE(t, y, u, v, w1, iota, nu,DOSE)
   

    % Define the system of ODEs
    dydt = zeros(length(y), 1);
    for i = 2:length(y)
        mitSkipped = i - 2;
        dydt(i) = alpha_p(mitSkipped, DOSE) * y(i) - k_p(mitSkipped, DOSE) * y(i) - a_p(mitSkipped, DOSE) * y(i);
        if mitSkipped >= 1
            dydt(i) = dydt(i) + k_p(mitSkipped - 1, DOSE) * y(i - 1);
        end
        % Sum up dead (apoptotic) cells
        dydt(1) = dydt(1) + a_p(mitSkipped, DOSE) * y(i);
    end

    % Define the functions for a, k_p and a_p
    function kp = k_p(i, x)
        kp = (nu - (1 - i)^2) * (x / (x + v));
    end

    function alphap = alpha_p(i, x)
        alphap = x / (x + u) * (i <= iota);
    end

    function ap = a_p(i, x)
        ap = (((2 - i) / 2))^2 * (x / (x + w1));
    end
end