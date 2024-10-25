function dydt = combined_ODECellCycleWithoutS(t, y, pars)

% Extract parameters
theta = pars(1);
nu = pars(2);
xi = pars(3);
a = pars(4);
v = pars(5);
w1 = pars(6);
iota = pars(7);
gamma_a = pars(8);
lambda = pars(9);
k_1 = pars(10);
k_2 = pars(11);
K = pars(12);
mu = pars(13);


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

%The Logistic term should only apply to the growth and cell cycle dynamics
%of the cells, not their drug-induced death or endoreplication.

%G1/S
dydt_P(2) = (1 - (P(2) + P(3)) / K) * (2 * k_1 * P(3) - k_2 * P(2)) - a_p(0, x) * P(2);

%G2/M
dydt_P(3) = (1 - (P(2) + P(3)) / K) * (k_2 * P(2) - k_1 * P(3)) - k_p(0, x) * P(3);

%P1
%Removed alpha_p(1, x) * P(4) because alpha always = 0 for P1 and P2
dydt_P(4) = - k_p(1, x) * P(4) + k_p(0, x) * P(3) - a_p(1, x) * P(4);
%Including influx from endoreplication

%P2
%Removed alpha_p(1, x) * P(5) because alpha always = 0 for P1 and P2
dydt_P(5) = - k_p(2, x) * P(5) + k_p(1, x) * P(4) - a_p(2, x) * P(5);

mitSkipped = [0, 0, 1, 2];
for i = 2:length(P)
    % Sum up dead (apoptotic) cells
    a_p_ = a_p(mitSkipped(i - 1), x)* P(i);
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
        kp = (i < 2) * nu * (exp(1) ^ (mu * -i)) * (x / (x + v));
    end

    function ap = a_p(i, x)
        ap = gamma_a * (x / (x + w1))^2 * exp(-lambda*i);
    end
end
