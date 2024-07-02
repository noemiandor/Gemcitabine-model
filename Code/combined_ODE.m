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
gamma_a = pars(9);
lambda = pars(10);


% Unpack the state variables
N = 2; % Number of state variables in gemcitabine_PKPD_ODE
P = y(N+1:end); % State variables for skippedMito_ODE
G = y(1); % State variable for gemcitabine_PKPD_ODE
x = y(2); %/sum(P(2:end)); % State variable for gemcitabine_PKPD_ODE
%disp([t,y(1),y(2), x])

% Initialize the derivative vector
dydt = zeros(length(y), 1);

% skippedMito_ODE
% Initialize the derivative for skippedMito_ODE
dydt_P = zeros(length(P), 1);
a_p_tot = 0;
for i = 2:length(P)
    mitSkipped = i - 2;
    dydt_P(i) = alpha_p(mitSkipped,x) * P(i) - k_p(mitSkipped, x) * P(i) - a_p(mitSkipped, x) * P(i);
    if mitSkipped >= 1
        dydt_P(i) = dydt_P(i) + k_p(mitSkipped - 1, x) * P(i - 1);
    end
    % Sum up dead (apoptotic) cells
    a_p_ = a_p(mitSkipped, x)* P(i);
    a_p_tot = a_p_tot + a_p_;
    dydt_P(1) = dydt_P(1) + a_p_ ;
end
a_p_tot = a_p_tot/ sum(P(2:end));
%disp([t,dydt_P'])

% Assign the derivatives for skippedMito_ODE to the combined dydt
dydt(N+1:end) = dydt_P;

% gemcitabine_PKPD_ODE
dydt(1) = max(-theta, -G); % * G * sum(P(2:end)); % dG/dt
dydt(2) = -dydt(1) - xi * a_p_tot * y(2); % dx/dt
% disp([t, dydt(2),  - xi * a_p_tot * y(2)]);

% Define the functions for k_p, alpha_p, and a_p
    function kp = k_p(i, x)
        kp = nu * (x / (x + v));
    end

    function alphap = alpha_p(i, x)
        alphap = (i <= iota) * eta * (1 - x / (x + a))^2 ;
    end

    function ap = a_p(i, x)
        ap = gamma_a * (x / (x + w1))^2 * exp(-lambda*i);
    end
end
