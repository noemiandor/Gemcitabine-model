function J  = cost(pars)
global dmx;

%% Parameters
u=pars(1);
w1=pars(2);
v=pars(3);
nu = pars(4);

%% @TODO: initial point should be a parameter

%% Define the initial conditions
tspan=(4:4:size(dmx.N2,2)*4);

%% Solve
[t_2N,y_2N] = ode45(@(t,y) skippedMito_ODE(t,y,u,v,w1, 0, nu), tspan, dmx.N2(:,1));
[t_4N,y_4N] = ode45(@(t,y) skippedMito_ODE(t,y,u,v,w1, 1, nu), tspan, dmx.N4(:,1));
fits=struct('N2',{t_2N,y_2N}, 'N4',{t_4N,y_4N});

%% Bring simulations and measurements to same dimensions
wsdAll=[];
idx=1;
for type=fieldnames(dmx)'
    t=getfield(fits,{1},type{1});
    y=getfield(fits,{2},type{1});
    dmx_=getfield(dmx,type{1});

    y_=y';
    t_=t;
    dmx__=dmx_;

    wsd=[];
    for ii =1:size(dmx__,1)
        wsd(ii) = ws_distance(y_(ii,:), dmx__(ii,:));
    end
    wsdAll=[wsdAll,mean(wsd)];

    %% plot
    mycolors = flip(colormap(parula(size(dmx__,1))));
    subplot(1,2,idx); idx=idx+1;
    hold off;
    for cI = 1:size(dmx__,1)
        plot(t_, y_(cI,:),'Color',mycolors(cI,:),'LineWidth',3)
%         set(gca, 'XScale', 'log')
%         set(gca, 'YScale', 'log')
%         ylim([1,1000])
        hold on;
        plot(t_, dmx__(cI,:),'*--','Color',mycolors(cI,:),'LineWidth',3)
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

end

J = mean(wsdAll);
end

