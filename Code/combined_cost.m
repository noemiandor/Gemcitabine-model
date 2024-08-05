function J = combined_cost(params, dmx_skippedMito, dmx_PKPD, GemcitabineConc_nM)
    
%% Parameters
    theta = params(1);
    nu_PKPD=params(2);
    eta = params(3);
    xi = params(4);

    u = params(5);
    v = params(6);
    w1 = params(7);
    iota=params(8);
    nu = params(9);
    
    % Parameters for gemcitabine_PKPD_ODE
    
    

   
    %% Define the initial conditions and time spans
    tspan_skippedMito = 2:2:size(dmx_skippedMito.N2, 2) * 2;
    tspan_PKPD = dmx_PKPD.time;
    % Define initial conditions for gemcitabine_PKPD_ODE
    initial_gemcitabine = [GemcitabineConc_nM.low, dmx_PKPD.low];

    % Define initial conditions for skippedMito_ODE
    initial_skippedMito_N2 = dmx_skippedMito.N2(:,1); % Initial conditions for N2
    initial_skippedMito_N4 = dmx_skippedMito.N4(:,1); % Initial conditions for N4

    % Combine initial conditions for combined_ODE
    initial_conditions_N2 = [initial_gemcitabine, initial_skippedMito_N2'];
    initial_conditions_N4 = [initial_gemcitabine, initial_skippedMito_N4'];

   skippedMito_params={u, v, w1, iota, nu,tspan_skippedMito,dmx_skippedMito};
   PKPD_params={theta, nu_PKPD, eta, xi,tspan_PKPD,dmx_PKPD, initial_gemcitabine}; 

    % Call ode45 for combined_ODE with N2 initial conditions
    [t_2N, y_2N] = ode45(@(t, y) combined_ODE(t, y, skippedMito_params,PKPD_params),...
        tspan_skippedMito, initial_skippedMito_N2);

    % Call ode45 for combined_ODE with N4 initial conditions
    [t_4N, y_4N] = ode45(@(t, y) combined_ODE(t, y, skippedMito_params,PKPD_params),...
        tspan_skippedMito, initial_skippedMito_N4);

    fits_skippedMito = struct('N2', [t_2N, y_2N], 'N4', [t_4N, y_4N]);




    %% Combine fits
  %    fits_combined = struct('skippedMito', fits_skippedMito, 'PKPD', fits_PKPD);

    %% Calculate weighted sum of distances for skippedMito
    wsdAll = [];
    idx = 1;
    for type = fieldnames(dmx_skippedMito)'
        t = fits_skippedMito.(type{1})(:, 1);
        y = fits_skippedMito.(type{1})(:, 2:end)';
        dmx_ = dmx_skippedMito.(type{1});

        y_=y;
        t_=t;
        dmx__=dmx_;

        wsd = zeros(size(dmx_, 1), 1);
        for ii = 1:size(dmx_, 1)
            wsd(ii) = mean(abs(y(ii, :) - dmx_(ii, :)));
        end
        wsdAll = [wsdAll; wsd];

        % Plotting code
        % Note: Consider using a separate function or section for plotting
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
    % @TODO include an optimization step for Gemcitabine. 

    % %% Calculate weighted sum of distances for PKPD
    % wsdAll_PKPD=[];
    % for type = fieldnames(fits_PKPD)'
    %     t = fits_PKPD.(type{1})(:, 1);
    %     y = fits_PKPD.(type{1})(:, 2:end)';
    %     dmx_ = dmx_PKPD.(type{1});
    % 
    %     wsd = mean(abs(y - dmx_));
    %     %wsdAll = [wsdAll; wsd];
    %     wsdAll_PKPD=[wsdAll_PKPD; wsd];
    % 
    %     % Plotting code
    %     % Note: Consider using a separate function or section for plotting
    % end
    % 
    % J = mean(wsdAll);
    % J_PKPD=mean(wsdAll_PKPD);
    J = mean(wsdAll);
end
