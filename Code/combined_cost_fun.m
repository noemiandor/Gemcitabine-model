%To perform the global optimization over all doses, we  modify 
% the script to incorporate a combined cost function that sums the costs for
% each dose in the set of all dose levels considered


%This function computes the sum of costs for all doses
function total_cost = combined_cost_fun(pars, dmx, DOSE)
    total_cost = 0;
    for dose = DOSE
        total_cost = total_cost + cost(pars, dmx, dose);
    end
end
