function [cost,drug_dose] = cost_PKPD(pars,dmx_PKPD,GemcitabineConc_nM)
%global dmx;
%global GemcitabineConc_nM

%% Parameters
theta=pars(1);
nu=pars(2);
eta=pars(3);
xi = pars(4);

disp('Cost fxn entered')


%% @TODO: initial point should be a parameter

%% Solve
[t_high,y_high]= ode45(@(t,y) gemcitabine_PKPD_ODE(t,y,theta, nu, eta, xi), dmx_PKPD.time, [GemcitabineConc_nM.high, dmx_PKPD.high(:,1)]);
[t_low,y_low] = ode45(@(t,y) gemcitabine_PKPD_ODE(t,y,theta, nu, eta, xi), dmx_PKPD.time, [GemcitabineConc_nM.low, dmx_PKPD.low(:,1)]);
fits=struct('low',{t_low,y_low}, 'high',{t_high,y_high});


%% Bring simulations and measurements to same dimensions
wsdAll=[];
idx=1;
% for type={'high'} %
for type=fieldnames(fits)'
    t=getfield(fits,{1},type{1});
    y=getfield(fits,{2},type{1});
    dmx_=getfield(dmx_PKPD,type{1});

    y_=y';
    t_=t;
    dmx__=dmx_;

    ii=find(~isnan(dmx__));
    wsd=ws_distance(y_(2,ii), dmx__(1,ii));
    wsdAll=[wsdAll,wsd/mean(dmx__(1,ii))];

    %% plot
    mycolors = flip(colormap(parula(size(dmx__,1))));
    subplot(1,2,idx); idx=idx+1;
    hold off;
    plot(t_, y_(2,:),'Color',mycolors(1,:),'LineWidth',3)
    %         set(gca, 'XScale', 'log')
    %         set(gca, 'YScale', 'log')
    %         ylim([1,1000])
    hold on;
    plot(t_, dmx__(1,:),'*--','Color',mycolors(1,:),'LineWidth',3)
    %         ylim([0,65])
    %         title(type{1})

    xlabel('Time (hours)')
    ylabel('dFdCTP')
    %     legend('D','$$\hat{D}$$','Interpreter','Latex')
    prefix='D';
    if strcmp(type{1},'high')
        prefix='T';
    end
    legend('A','$$\hat{A}$$',['$$',prefix,'_0$$'],['$$\hat{',prefix,'_0}$$'],['$$',prefix,'_1$$'],['$$\hat{',prefix,'_1}$$'],['$$',prefix,'_2$$'],['$$\hat{',prefix,'_2}$$'],'Interpreter','Latex')
    

end
drug_dose=y();
cost= mean(wsdAll);
%dose=getfield(fits,{2},type{1});
%dose=t
end

