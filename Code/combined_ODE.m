function dydt = combined_ODE(t, y, pars)

    % Extract parameters
    theta = pars(1); 
    nu = pars(2); 
    eta = pars(3); 
    xi = pars(4);
    a = pars(5); 
    v = pars(6); 
    w1 = pars(7); 
    iota = pars(8); 

    % Unpack the state variables
    N = 2; % Number of state variables in gemcitabine_PKPD_ODE
    G = y(1); % State variable for gemcitabine_PKPD_ODE
    x = y(2); % State variable for gemcitabine_PKPD_ODE
    P = y(N+1:end); % State variables for skippedMito_ODE

    % Initialize the derivative vector
    dydt = zeros(length(y), 1);

    % gemcitabine_PKPD_ODE
    dydt(1) = -theta * G / (G + nu) * sum(P(2:end)); % dG/dt
    dydt(2) = eta * theta * G / (G + nu) - xi * x; % dx/dt

    % skippedMito_ODE
    % Initialize the derivative for skippedMito_ODE
    dydt_P = zeros(length(P), 1);
    for i = 2:length(P)
        mitSkipped = i - 2;
        dydt_P(i) = alpha_p(mitSkipped, x) * P(i) - k_p(mitSkipped, x) * P(i) - a_p(mitSkipped, x) * P(i);
        if mitSkipped >= 1
            dydt_P(i) = dydt_P(i) + k_p(mitSkipped - 1, x) * P(i - 1);
        end
        % Sum up dead (apoptotic) cells
        dydt_P(1) = dydt_P(1) + a_p(mitSkipped, x) * P(i);
    end

    % Assign the derivatives for skippedMito_ODE to the combined dydt
    dydt(N+1:end) = dydt_P;

      % Define the functions for k_p, alpha_p, and a_p
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


% function dydt = combined_ODE(t, y, pars)
% 
%  % Extract parameters
%     theta = pars(1); nu= pars(2); eta = pars(3); xi = pars(4);
%     a = pars(5); v = pars(6); w1 = pars(7); iota = pars(8); 
% 
%     % Define ODE model parameters
%     % params = [theta, nu_gemcitabine, eta, xi, a, v, w1, iota, nu_skippedMito];
% 
%     % Unpack the state variables
%     N = 2; % Number of state variables in gemcitabine_PKPD_ODE
%     y1 = y(1:N); % State variables for gemcitabine_PKPD_ODE
%     y2 = y(N+1:end); % State variables for skippedMito_ODE
% 
%     % Initialize the derivative vector
%     dydt = zeros(length(y), 1);
% 
%     % gemcitabine_PKPD_ODE
% 
%     dydt(1) = -theta * y1(1) / (y1(1) + nu); %@TODO: *N
%     dydt(2) = eta * theta * y1(1) / (y1(1) + nu) - xi * y1(2); % This is `DOSE` in skippedMito_ODE
% 
%     % Define DOSE as y1(2) for skippedMito_ODE
%     DOSE = y1(2);
% 
%     % skippedMito_ODE
%     % Initialize the derivative for skippedMito_ODE
%     dydt_skippedMito = zeros(length(y2), 1);
%     for i = 2:length(y2)
%         mitSkipped = i - 2;
%         dydt_skippedMito(i) = alpha_p(mitSkipped, DOSE) * y2(i) - k_p(mitSkipped, DOSE) * y2(i) - a_p(mitSkipped, DOSE) * y2(i);
%         if mitSkipped >= 1
%             dydt_skippedMito(i) = dydt_skippedMito(i) + k_p(mitSkipped - 1, DOSE) * y2(i - 1);
%         end
%         % Sum up dead (apoptotic) cells
%         dydt_skippedMito(1) = dydt_skippedMito(1) + a_p(mitSkipped, DOSE) * y2(i);
%     end
% 
%     % Assign the derivatives for skippedMito_ODE to the combined dydt
%     dydt(N+1:end) = dydt_skippedMito;
% 
%     % Define the functions for k_p, alpha_p, and a_p
%     function kp = k_p(i, x)
%         kp = (nu - (1 - i)^2) * (x / (x + v));
%     end
% 
%     function alphap = alpha_p(i, x)
%         alphap = x / (x + a) * (i <= iota);
%     end
% 
%     function ap = a_p(i, x)
%         ap = (((2 - i) / 2))^2 * (x / (x + w1));
%     end
% end
