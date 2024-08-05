function J = combined_cost(pars, dmx_skippedMito, dmx_PKPD, GemcitabineConc_nM)
 %disp('Entering combined_cost');
    %% Parameters
    % Parameters for skippedMito_ODE
    u = pars(1);
    v = pars(2);
    w1 = pars(3);
    nu = pars(4);
    
    % Parameters for gemcitabine_PKPD_ODE
    theta = pars(5);
    nu_PKPD=pars(6);
    eta = pars(7);
    xi = pars(8);
    
    DOSE=200;

    %% Define the initial conditions and time spans
    tspan_skippedMito = 2:2:size(dmx_skippedMito.N2, 2) * 2;
    tspan_PKPD = dmx_PKPD.time;

    %% Solve skippedMito_ODE
    [t_2N,y_2N] = ode45(@(t,y) skippedMito_ODE(t,y,u,v,w1, 0, nu,DOSE), tspan_skippedMito, dmx_skippedMito.N2(:,1));
    [t_4N,y_4N] = ode45(@(t,y) skippedMito_ODE(t,y,u,v,w1, 1, nu,DOSE), tspan_skippedMito, dmx_skippedMito.N4(:,1));
%     [t_2N, y_2N] = ode45(@(t,y) skippedMito_ODE(t, y, theta, nu, eta, xi, u, v, w1, 0), tspan_skippedMito, [GemcitabineConc_nM.low, dmx_skippedMito.N2(:,1)']);
%     [t_4N, y_4N] = ode45(@(t,y) skippedMito_ODE(t, y, theta, nu, eta, xi, u, v, w1, 1), tspan_skippedMito, [GemcitabineConc_nM.low, dmx_skippedMito.N4(:,1)']);
    fits_skippedMito = struct('N2', [t_2N, y_2N], 'N4', [t_4N, y_4N]);


    %% Solve gemcitabine_PKPD_ODE
%      [t_high, y_high] = ode45(@(t,y) combined_ODE(t, y, theta, nu, eta, xi, u, v, w1, 0), tspan_PKPD, [GemcitabineConc_nM.high, dmx_PKPD.high]);
%     [t_low, y_low] = ode45(@(t,y) combined_ODE(t, y, theta, nu, eta, xi, u, v, w1, 0), tspan_PKPD, [GemcitabineConc_nM.low, dmx_PKPD.low]);
%     fits_PKPD = struct('low', [t_low, y_low], 'high', [t_high, y_high]);
% 
    [t_high, y_high] = ode45(@(t,y) gemcitabine_PKPD_ODE(t, y, theta,nu_PKPD, eta, xi), tspan_PKPD, [GemcitabineConc_nM.high, dmx_PKPD.high]);
    [t_low, y_low] = ode45(@(t,y)   gemcitabine_PKPD_ODE(t, y, theta, nu_PKPD,eta, xi), tspan_PKPD, [GemcitabineConc_nM.low, dmx_PKPD.low]);
    fits_PKPD = struct('low', [t_low, y_low], 'high', [t_high, y_high]);

  
    %% Combine fits
    fits_combined = struct('skippedMito', fits_skippedMito, 'PKPD', fits_PKPD);

    %% Calculate weighted sum of distances for skippedMito
    wsdAll = [];
    
        
    idx = 1;
    for type=fieldnames(dmx_skippedMito)'
    %t=getfield(fits_skippedMito,{1},type{1});
    t=fits_skippedMito.(type{1})(:,1);
    %y=getfield(fits_skippedMito,{2},type{1})
    y=fits_skippedMito.(type{1})(:,2:end);
    %fprintf('y: %d\n',t );
    %dmx_=getfield(dmx_skippedMito,type{1});
    dmx_=dmx_skippedMito.(type{1});


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
        disp('supposed to be plotting')
        
        plot(t_', y_(cI,:),'Color',mycolors(cI,:),'LineWidth',3)
%         set(gca, 'XScale', 'log')
%         set(gca, 'YScale', 'log')
%         ylim([1,1000])
        hold on;
        plot(t_', dmx__(cI,:),'*--','Color',mycolors(cI,:),'LineWidth',3)
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
   

% 
   J = mean(wsdAll);

    
    
end
