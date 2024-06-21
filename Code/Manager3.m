

cd('/Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Data/M00_GemcitabinePKPD_101823')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incucyte data: treated conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('/Users/4477116/My Drive/IMO/Anaconda/Git/Gemcitabine-model/Code')
addpath /Users/4477116/My Drive/IMO/Anaconda/Git/Gemcitabine-model/Code/
addpath /Users/4477116/My Drive/IMO/Anaconda/Git/Gemcitabine-model/Code/wassersteinFun/
addpath /Users/4477116/My Drive/IMO/Anaconda/Git/Gemcitabine-model/Data/matlab
addpath /Users/4477116/My Drive/IMO/Anaconda/Git/Gemcitabine-model/Data/M00_GemcitabinePKPD_101823
replicates=struct('N2',{'B_4'},...
                 'N4',{'E_4'});

DOSE=20;

% Initialize shared data structures
dmx_data = cell(length(dir('nM1000_*txt')), 1);
pars_data = cell(length(dir('nM1000_*txt')), 1);
S_data = struct();

% Iterate across replicates
for k = 1:length(dir('nM1000_*txt'))
    % Load data
    dm = readtable(['nM1000_', num2str(k), '.txt']);

    % Prepare data structure
    dmx = struct();
    dmx.high = dm.dFdCTP___ng_mL_';   
    dmx.low = dm.dFdCTP___ng_mL__low';
    dmx.time = dm.time';

    % Fit model to data
    pars = fit_model(dmx, GemcitabineConc_nM);
    
    % Store computed data
    dmx_data{k} = dmx;
    pars_data{k} = pars;
    
    % Plot best fit
    plot_best_fit(pars, dmx, GemcitabineConc_nM);
    
    % Update S_data
    fname = ['nM1000_', num2str(k)];
    S_data.(fname) = pars;
end

% Process Incucyte data
for k = 1:length(replicates)
    % Load data
    dm = readtable([replicates(k).N2, '.txt']);
    dmx = struct();
    dmx.N2 = table2array(dm);
    
    % Fit model to data
    pars = fit_model(dmx, DOSE);
    
    % Store computed data
    dmx_data{k + length(dir('nM1000_*txt'))} = dmx;
    pars_data{k + length(dir('nM1000_*txt'))} = pars;
    
    % Update S_data
    S_data.(replicates(k).N2) = pars;
    
    % Plot best fit
    plot_best_fit(pars, dmx, DOSE);
end

% Functions definitions
function pars = fit_model(dmx, conc_data)
    % Define parameters and options
    theta = 50;
    nu = 1500;
    eta = 0.005;
    xi = 0.05;
    pars = [theta, nu, eta, xi];
    lb = pars / 5000;
    ub = pars * 5000;
    opts = optimoptions(@fmincon);

    % Create optimization problem
    problem = createOptimProblem('fmincon', 'objective', @(pars) cost_PKPD(pars, dmx, conc_data), ...
                                  'x0', pars, 'lb', lb, 'ub', ub, 'options', opts);

    % Run optimization
    ms = MultiStart('UseParallel', true);
    [pars, ~, ~, ~, ~] = run(ms, problem);
end

function plot_best_fit(pars, dmx, conc_data)
    figure;
    cost_PKPD(pars, dmx, conc_data);
    % Customize plot if needed
end
