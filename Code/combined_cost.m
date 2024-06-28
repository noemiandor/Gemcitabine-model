function cost = combined_cost(params, dmx, GemcitabineConc_nM)

% Initial conditions
initial_conditions_N2 = [GemcitabineConc_nM; dmx.N2dFdCTP(1); dmx.N2(:,1)];
initial_conditions_N4 = [GemcitabineConc_nM; dmx.N4dFdCTP(1); dmx.N4(:,1)];

% Simulate the combined ODE model for N2 population
[t_N2, y_N2] = ode45(@(t, y) combined_ODE(t, y, params), dmx.cellTime, initial_conditions_N2);
model_N2 = y_N2(:, 1);
model_drugN2 = y_N2(:, 4);

% Simulate the combined ODE model for N4 population
[t_N4, y_N4] = ode45(@(t, y) combined_ODE(t, y, params), dmx.cellTime, initial_conditions_N4);
model_N4 = y_N4(:, 1);
model_drugN4 = y_N4(:, 4);

% @TODO: we should probably use the wasserstein distance again
% Calculate the cost as the sum of squared differences between model and data
cost_N2 = sum(sum((model_N2 - dmx.N2').^2));
cost_N4 = sum(sum((model_N4 - dmx.N4').^2));
cost_drugN2 = sum((model_drugN2- dmx.N2dFdCTP').^2);
cost_drugN4 = sum((model_drugN4 - dmx.N4dFdCTP').^2);

%% plotting 
fits=struct('N2',{t_N2,y_N2}, 'N4',{t_N4,y_N4});
idx=1;
typecolors=struct('N2','black','N4','red');
for type= {'N2', 'N4'}
    t=getfield(fits,{1},type{1});
    y=getfield(fits,{2},type{1});
    dmx_=getfield(dmx,type{1});

    %% plot
    mycolors = flip(colormap(parula(size(dmx_,1))));
    subplot(1,3,idx); idx=idx+1;
    hold off;
    for cI = 1:size(dmx_,1)
        plot(t, y(:,cI),'Color',mycolors(cI,:),'LineWidth',3)
        %         set(gca, 'XScale', 'log')
        %         set(gca, 'YScale', 'log')
        %         ylim([1,1000])
        hold on;
        plot(t, dmx_(cI,:),'*--','Color',mycolors(cI,:),'LineWidth',3)
        %         ylim([0,65])
        %         title(type{1})
    end
    xlabel('Time (hours)')
    ylabel('Number of cells')
    %     legend('D','$$\hat{D}$$','Interpreter','Latex')
    prefix='D';
    if strcmp(type{1},'N4')
        prefix='T';
    end
    legend('A','$$\hat{A}$$',['$$',prefix,'_0$$'],['$$\hat{',prefix,'_0}$$'],['$$',prefix,'_1$$'],['$$\hat{',prefix,'_1}$$'],['$$',prefix,'_2$$'],['$$\hat{',prefix,'_2}$$'],'Interpreter','Latex')
    %% plot dFdCTP
    subplot(1,3,3)
    plot(t, y(:,end),'Color',getfield(typecolors,type{1}),'LineWidth',3)
    hold on;
    plot(t, getfield(dmx,[type{1},'dFdCTP']),'*--','Color',getfield(typecolors,type{1}),'LineWidth',3)

end
legend([fieldnames(typecolors);fieldnames(typecolors)])


% Combine the costs for N2 and N4
cost = cost_N2 + cost_N4 + cost_drugN2 + cost_drugN4;


end

