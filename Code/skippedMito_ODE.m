function dydt = skippedMito_ODE(t, y, a, v, w1, iota, nu)
global DOSE
% Define the system of ODEs
dydt = zeros(length(y), 1);
for i = 2:length(y)
    mitSkipped=i-2;
    dydt(i) = alpha_p(mitSkipped,DOSE)*y(i) - k_p(mitSkipped,DOSE)*y(i) - a_p(mitSkipped,DOSE)*y(i);
    if mitSkipped>=1
        dydt(i) = dydt(i) + k_p(mitSkipped-1,DOSE)*y(i-1);
    end
    %% sum up dead (apoptotic) cells
    dydt(1) = dydt(1) + a_p(mitSkipped,DOSE)*y(i);
end

% Define the functions for a, k_p and a_p
    function kp = k_p(i,x)
%       kp= (1-abs(1-i))*(x/(x+v));
        kp= (nu-(1-i)^2)*(x/(x+v));
    end

    function alphap = alpha_p(i,x)
        alphap = x/(x+a) *(i<=iota);
%         alphap = (1-x/(x+a)) *(i<=iota); %% goodness of fit is same but
%         converges slower -- @TODO: matters once we fit to >1 drug
%         concentrations
    end

    function ap = a_p(i,x)
        ap= (((2-i)/2))^2 *(x/(x+w1));
        %     ap = w1*x*(2-i)^2;
        %     ap = (w1^(1/(i+1)^2)* x/(x+w2));
    end
end