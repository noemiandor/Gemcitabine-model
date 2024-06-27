function cost = combined_cost(params, dmx)
   
    % Initial conditions
    initial_conditions_N2 = [dmx.N2(:,1); dmx.N2dFdCTP(1)];
    initial_conditions_N4 = [dmx.N4(:,1); dmx.N4dFdCTP(1)];

    % Simulate the combined ODE model for N2 population
    [~, y_N2] = ode45(@(t, y) combined_ODE(t, y, params), dmx.cellTime, initial_conditions_N2);
    model_N2 = y_N2(:, 1);
    model_drugN2 = y_N2(:, 4);

    % Simulate the combined ODE model for N4 population
    [~, y_N4] = ode45(@(t, y) combined_ODE(t, y, params), dmx.cellTime, initial_conditions_N4);
    model_N4 = y_N4(:, 1);
    model_drugN4 = y_N4(:, 4);

    % @TODO: we should probably use the wasserstein distance again
    % Calculate the cost as the sum of squared differences between model and data
    cost_N2 = sum(sum((model_N2 - dmx.N2').^2));
    cost_N4 = sum(sum((model_N4 - dmx.N4').^2));
    cost_drugN2 = sum((model_drugN2- dmx.N2dFdCTP').^2);
    cost_drugN4 = sum((model_drugN4 - dmx.N4dFdCTP').^2);

    % Combine the costs for N2 and N4
    cost = cost_N2 + cost_N4 + cost_drugN2 + cost_drugN4;

    % @TODO: plotting needs to be restored
end

