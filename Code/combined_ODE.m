function dydt_skippedMito = combined_ODE(t, y_skippedMito, skippedMito_params,PKPD_params)


   %Set parameters
    u=skippedMito_params{1};
    v=skippedMito_params{2};
    w1=skippedMito_params{3};
    iota=skippedMito_params{4};
    nu=skippedMito_params{5};
    tspan_skippedMito=skippedMito_params{6};
    dmx_skippedMito=skippedMito_params{7};

    theta=PKPD_params{1};
    nu_PKPD=PKPD_params{2};
    eta=PKPD_params{3};
    xi=PKPD_params{4};
    tspan_PKPD=PKPD_params{5};
    dmx_PKPD=PKPD_params{6};
    initial_gemcitabine=PKPD_params{7};


    %@TODO DOSE level are computed from the same IC and passed on to Cell
    %model, the next stemp is to update the gemcitabine model differently
    %(with new data) at each step


    % Define initial conditions for gemcitabine_PKPD_ODE
    y_gemcitabine_initial = initial_gemcitabine;  % Example initial conditions for gemcitabine_PKPD_ODE

    % Solve gemcitabine_PKPD_ODE at the current time step
    [~, y_gemcitabine] = ode45(@(t_gem, y_gem) gemcitabine_PKPD_ODE(t_gem, y_gem, theta, nu_PKPD, eta, xi),  tspan_PKPD, y_gemcitabine_initial);

    % Get the last value of y_gemcitabine as DOSE
    DOSE = y_gemcitabine(end, 2);

    % Compute derivatives for skippedMito_ODE
    dydt_skippedMito = skippedMito_ODE(t, y_skippedMito, u, v, w1, iota, nu, DOSE);

    % Nested function for gemcitabine_PKPD_ODE
    function dydt_gemcitabine = gemcitabine_PKPD_ODE(t_gem, y_gem, theta, nu, eta, xi)
        dydt_gemcitabine = zeros(size(y_gem));
        dydt_gemcitabine(1) = -theta * y_gem(1) / (y_gem(1) + nu);
        dydt_gemcitabine(2) = eta * theta * y_gem(1) / (y_gem(1) + nu) - xi * y_gem(2);
    end

    % Nested function for skippedMito_ODE
    function dydt = skippedMito_ODE(t, y, a, v, w1, iota, nu, DOSE)
        dydt = zeros(length(y), 1);
        for i = 2:length(y)
            mitSkipped = i - 2;
            dydt(i) = alpha_p(mitSkipped, DOSE) * y(i) - k_p(mitSkipped, DOSE) * y(i) - a_p(mitSkipped, DOSE) * y(i);
            if mitSkipped >= 1
                dydt(i) = dydt(i) + k_p(mitSkipped - 1, DOSE) * y(i - 1);
            end
            %% sum up dead (apoptotic) cells
            dydt(1) = dydt(1) + a_p(mitSkipped, DOSE) * y(i);
        end

        % Define the functions for a, k_p, and a_p
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
end
