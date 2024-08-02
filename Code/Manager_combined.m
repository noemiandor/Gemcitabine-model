% path2celldata='/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_SkippedMitosisClassification_042523';
path2celldata='~/Repositories/Gemcitabine-model/Data/matlab/';
path2PKPDdata='/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/M00_GemcitabinePKPD_101823';
path2root='/Users/4470246/Repositories/Gemcitabine-model/Code';
cd(path2root)
addpath(path2root)

global dmx

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit model to data %%%
S = struct();
% @TODO: we need to read in the plate map here
% @TODO: we need to use wells that have the same (or similar) normalized drug concentration as the starting concentration in the PKPD experiment
replicates = struct('N2', {'B_4', 'B_5', 'B_6'}, ...
    'N4', {'E_4', 'E_5', 'E_6'});
dfdctpConc_incucyte_nM=15; %@TODO: we need to figure out to which intra-cellular dFdCTP concentration each well corresponds at timepoint zero 

%% Iterate across replicates
for k = 1:length(replicates.N2)

    % Read PKPD data
    cd(path2PKPDdata)
    dmx = struct();
    pkpd_file = dir(['nM1000_2N-', num2str(k), '.txt']);
    pkpd_data = readtable(pkpd_file.name);
    dmx.N2dFdCTP= pkpd_data.dFdCTP___ng_mL_';
    dmx.PKPDtime = pkpd_data.time';

    % Read cell data for N2
    cd(path2celldata)
    dmx.N2 = readtable([replicates(k).N2, '.txt']);
    dmx.N2 = table2array(dmx.N2(:, 2:end)); % assuming data starts from the second column
    dmx.cellTime = 2:2:size(dmx.N2, 2) * 2;
    % plot(dmx.cellTime, dmx.N2(2,:)); xlabel('hour'); ylabel('# diploid cells')
    cd(path2root)

    %% correct data
    ii = find(dmx.cellTime<55);
    N2_interpolated = interp1(dmx.cellTime(ii), dmx.N2(2,ii), dmx.cellTime, 'linear', 'extrap');
    dmx.N2(2,:)=N2_interpolated;

    % Interpolate PKPD data to match cellTime points
    N2_interpolated = interp1(dmx.PKPDtime, dmx.N2dFdCTP, dmx.cellTime, 'linear', 'extrap');
    N2_interpolated(N2_interpolated<0)=0;
   
    % Visualize the interpolation
    figure;
    plot(dmx.PKPDtime, dmx.N2dFdCTP, 'o-', 'DisplayName', 'Original High');
    hold on;
    plot(dmx.cellTime, N2_interpolated, '*-', 'DisplayName', 'Interpolated High');
    title('High Concentration Interpolation');
    xlabel('Time');
    ylabel('Concentration (ng/mL)');
    legend;

    % Match pkpd drug concentration to incucyte drug concentration
    [~,ia]=min(abs(dfdctpConc_incucyte_nM-N2_interpolated));
    N2_interpolated = interp1(dmx.cellTime(1:(1+end-ia)), N2_interpolated(ia:end), dmx.cellTime, 'linear', 'extrap');
    dmx.N2dFdCTP = N2_interpolated;

    %% Do the fitting for one cell type at a time
    ii = 1:length(dmx.cellTime); %find(dmx.cellTime<=55);
    dmx.dFdCTP = dmx.N2dFdCTP(ii);
    dmx.cells = dmx.N2(:,ii);
    % need to add zero for cells of highest ploidy + 1 orelse the model will not penalize their existence
    dmx.cells = [dmx.cells; zeros(1, size(dmx.cells,2))];
    dmx.cellTime = dmx.cellTime(ii);
   
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % Combined parameters for both models
    theta = 5e-1; nu = 4.0378; 
    eta =  0.3340; xi = 7.3125; a = 1.0000; 
    v = 381.4379; w1 = 20; iota = 0;
    gamma_a = 0.1000; lambda = 1.0000;
    pars = {theta, nu, eta, xi, a, v, w1, iota, gamma_a, lambda };
    bounds = cell2mat(cellfun(@(x) [x/10; x*10], pars, 'UniformOutput', false));
    % bounds(:,2) = [pars{2}/10; pars{2}*10];
    % bounds(:,6) = [pars{6}/10; pars{6}*10];
    lb = bounds(1,:)';
    ub = bounds(2,:)';

    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon', 'objective', ...
        @(pars) combined_cost(pars, dmx), 'x0', cell2mat(pars), 'lb', lb, 'ub', ub, 'options', opts);

    rs = RandomStartPointSet('NumStartPoints', 5);
    points = list(rs, problem);
    ms = MultiStart('UseParallel', true);
    [pars_, fval, exitflag, output, solutions] = run(ms, problem, CustomStartPointSet(points));

    S = setfield(S, replicates(k).N2, pars_);
    % save('solutions.mat','pars_');
    % pars=arrayfun(@(x) {x}, pars_);

    %% Plot best fit:
    % close all hidden
    figure('name', ['~/Downloads/Combined_model_', replicates(k).N2], 'Position', [100 100 1000 400])
    combined_cost(getfield(S, replicates(k).N2), dmx)
end
