function cost = combined_cost(params, dmx)

% Initial conditions
initial_conditions = [ dmx.dFdCTP(1); 0; dmx.cells(:,1)];

% Simulate the combined ODE model for population
[t_N, y_N] = ode15s(@(t, y) combined_ODE(t, y, params), dmx.cellTime, initial_conditions);
model_drugN= y_N(:, 1:2);
model_N  = y_N(:, 3:end);
if size(model_N,1)<size(dmx.cells,2)
    cost = inf;
    return;
end

% @TODO: we should probably use the wasserstein distance again
% Calculate the cost as the sum of squared differences between model and data
cost_N = sum(sum((model_N(:,1:end) - dmx.cells(1:end,:)').^2)) / max(max((dmx.cells)))^2;
cost_drugN = sum((model_drugN(:,1)- dmx.dFdCTP').^2) / max(dmx.dFdCTP)^2;


%% plot
mycolors = flip(colormap(parula(size(dmx.cells,1))));
subplot(1,2,1);
hold off;
for cI = 1:size(dmx.cells,1)
    plot(t_N, model_N(:,cI),'Color',mycolors(cI,:),'LineWidth',3)
    %         set(gca, 'XScale', 'log')
    %         set(gca, 'YScale', 'log')
    %         ylim([1,1000])
    hold on;
    plot(t_N, dmx.cells(cI,:),'*--','Color',mycolors(cI,:),'LineWidth',3)
    %         ylim([0,65])
    %         title(type{1})

    xlabel('Time (hours)')
    ylabel('Number of cells')
    %     legend('D','$$\hat{D}$$','Interpreter','Latex')
    prefix='D';
    legend('A','$$\hat{A}$$',['$$',prefix,'_0$$'],['$$\hat{',prefix,'_0}$$'],['$$',prefix,'_1$$'],['$$\hat{',prefix,'_1}$$'],['$$',prefix,'_2$$'],['$$\hat{',prefix,'_2}$$'],'Interpreter','Latex')
end
%% plot dFdCTP
subplot(1,2,2); hold on;
plot(t_N, dmx.dFdCTP,'*--','Color','black','LineWidth',3)
plot(t_N, model_drugN(:,1),'Color','black','LineWidth',3)
yyaxis right
plot(t_N, model_drugN(:,2),'Color','red','LineWidth',3)
legend({'dFdCTP__data','dFdCTP__model','dna__dfdctp'})

% Combine the costs for cells and drug
cost = cost_drugN + cost_N; 


end

