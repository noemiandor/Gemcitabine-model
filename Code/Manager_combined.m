% path2celldata='/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_SkippedMitosisClassification_042523';
path2celldata='~/Repositories/Gemcitabine-model/Data/matlab/';
path2PKPDdata='/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/M00_GemcitabinePKPD_101823';
path2root='/Users/4470246/Repositories/Gemcitabine-model/Code';
cd(path2root)
addpath(path2root)

global dmx
global GemcitabineConc_nM
GemcitabineConc_nM=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit model to data %%%
S = struct();
% @TODO: we need to read in the plate map here 
% @TODO: we need to somehow normalize the drug concentrations in the PKPD data to the cell concentration in the PKPD experiment (which is different from the cell concentration in the Incucyte experiment)
% @TODO: we need to use wells that have the same (or similar) normalized drug concentration as the starting concentration in the PKPD experiment
replicates = struct('N2', {'B_4', 'B_5', 'B_6'}, ...
                    'N4', {'E_4', 'E_5', 'E_6'});

%% Iterate across replicates
for k = 1:length(replicates.N2)

    % Read PKPD data
    cd(path2PKPDdata)
    dmx = struct();
    pkpd_file = dir(['nM1000_2N-', num2str(k), '.txt']);
    pkpd_data = readtable(pkpd_file.name);
    dmx.N2dFdCTP= pkpd_data.dFdCTP___ng_mL_';
    pkpd_file = dir(['nM1000_4N-', num2str(k), '.txt']);
    pkpd_data = readtable(pkpd_file.name);
    dmx.N4dFdCTP = pkpd_data.dFdCTP___ng_mL_';
    dmx.PKPDtime = pkpd_data.time';

    % Read cell data for N2 and N4
    cd(path2celldata)
    dmx.N2 = readtable([replicates(k).N2, '.txt']);
    dmx.N2 = table2array(dmx.N2(:, 2:end)); % assuming data starts from the second column
    dmx.N4 = readtable([replicates(k).N4, '.txt']);
    dmx.N4 = table2array(dmx.N4(:, 2:end)); % assuming data starts from the second column
    dmx.cellTime = 2:2:size(dmx.N4, 2) * 2;
    cd(path2root)

    % Interpolate PKPD data to match cellTime points
    N2_interpolated = interp1(dmx.PKPDtime, dmx.N2dFdCTP, dmx.cellTime, 'linear', 'extrap');
    N4_interpolated = interp1(dmx.PKPDtime, dmx.N4dFdCTP, dmx.cellTime, 'linear', 'extrap');
    N2_interpolated(N2_interpolated<0)=0;
    N4_interpolated(N4_interpolated<0)=0;

    % Visualize the interpolation
    figure;
    subplot(1, 2, 1);
    plot(dmx.PKPDtime, dmx.N2dFdCTP, 'o-', 'DisplayName', 'Original High');
    hold on;
    plot(dmx.cellTime, N2_interpolated, '*-', 'DisplayName', 'Interpolated High');
    title('High Concentration Interpolation');
    xlabel('Time');
    ylabel('Concentration (ng/mL)');
    legend;

    subplot(1, 2, 2);
    plot(dmx.PKPDtime, dmx.N4dFdCTP, 'o-', 'DisplayName', 'Original Low');
    hold on;
    plot(dmx.cellTime, N4_interpolated, '*-', 'DisplayName', 'Interpolated Low');
    title('Low Concentration Interpolation');
    xlabel('Time');
    ylabel('Concentration (ng/mL)');
    legend;
    
    % Update dmx with interpolated values
    dmx.N2dFdCTP = N2_interpolated;
    dmx.N4dFdCTP = N4_interpolated;

    %% Do the fitting
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % Combined parameters for both models
    theta = 50; nu_gemcitabine = 1500; eta = 0.005; xi = 0.05;
    a = 150; v = 65; w1 = 310; iota = 0;
    pars = {theta, nu_gemcitabine, eta, xi, a, v, w1, iota};
    bounds = cell2mat(cellfun(@(x) [x/5000; x*5000], pars, 'UniformOutput', false));
    lb = bounds(1,:)';
    ub = bounds(2,:)';

    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon', 'objective', ...
        @(pars) combined_cost(pars, dmx, GemcitabineConc_nM), 'x0', cell2mat(pars), 'lb', lb, 'ub', ub, 'options', opts);

    rs = RandomStartPointSet('NumStartPoints', 2);
    points = list(rs, problem);
    ms = MultiStart('UseParallel', true);
    [pars_, fval, exitflag, output, solutions] = run(ms, problem, CustomStartPointSet(points));

    fname = [replicates(k).N2,'_and_', replicates(k).N4];
    S = setfield(S, fname, pars_);

    %% Plot best fit:
    % close all hidden
    figure('name', ['~/Downloads/Combined_model_', fname], 'Position', [100 100 1000 400])
    combined_cost(getfield(S, fname), dmx)
end
