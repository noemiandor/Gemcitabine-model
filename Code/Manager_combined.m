clear workspace

% path2celldata='/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_SkippedMitosisClassification_042523';
path2celldata='C:\Users\chayn\OneDrive\Documents\USF Moffitt Cancer\AndorLabCode\Gemcitabine-model\Data\matlab';
path2PKPDdata='C:\Users\chayn\OneDrive\Documents\USF Moffitt Cancer\AndorLabCode\Gemcitabine-model\Data\M00_GemcitabinePKPD_101823';
path2root='C:\Users\chayn\OneDrive\Documents\USF Moffitt Cancer\AndorLabCode\Gemcitabine-model\Code';
path2cellCycledata='C:\Users\chayn\OneDrive\Documents\USF Moffitt Cancer\AndorLabCode\Gemcitabine-model\Data\FilteredAndMultipliedData';
cd(path2root)
addpath(path2root)

global dmx

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit model to data %%%
S = struct();
% @TODO: we need to read in the plate map here
% @TODO: we need to use wells that have the same (or similar) normalized drug concentration as the starting concentration in the PKPD experiment

N2Cells = {'A_2', 'A_3', 'A_4'};
N4Cells = {'E_2', 'E_2', 'E_4'};

concentrationMap = [NaN, 0.001, 3.125, 6.25, 12.5, 25, 50, 100, 200, 400, 800, NaN];

replicates = struct('N2', N2Cells, ...
    'N4', N4Cells);
dfdctpConc_incucyte_nM=15; %@TODO: we need to figure out to which intra-cellular dFdCTP concentration each well corresponds at timepoint zero 


cd(path2PKPDdata)
dmx = struct();
pkpdFile1 = dir(['nM1000_2N-', num2str(1), '.txt']);
pkpdFile2 = dir(['nM1000_2N-', num2str(2), '.txt']);
pkpdFile3 = dir(['nM1000_2N-', num2str(3), '.txt']);

pkpdData1 = readtable(pkpdFile1.name);
pkpdData2 = readtable(pkpdFile2.name);
pkpdData3 = readtable(pkpdFile3.name);

%Now that the three pharmacodynamics replicates are read in, we can take
%their average to find the average pharmacodynamics of gemcitabine in the
%2N population: N2dFdCTP

dmx.N2dFdCTP = (pkpdData1.dFdCTP___ng_mL_' + pkpdData2.dFdCTP___ng_mL_' + pkpdData3.dFdCTP___ng_mL_') / 3;
dmx.N2dFdCTP = dmx.N2dFdCTP - dmx.N2dFdCTP(1);

negativedFdCTPTimes = find(dmx.N2dFdCTP <= 0);

dmx.N2dFdCTP(negativedFdCTPTimes) = 0;

dmx.PKPDtime = pkpdData1.time';

length(N2Cells)

%% Iterate across replicates
for k = 1:1

    % Read PKPD data
    

    % Read cell data for N2
    cd(path2cellCycledata)
    
    N2Well = replicates(k).N2;

    dmx.N2 = readtable([N2Well, '.txt']);
    dmx.N2 = table2array(dmx.N2(:, 2:end)); % assuming data starts from the second column
    dmx.cellTime = 2:2:size(dmx.N2, 2) * 2;
    
    %At a point near the end, all of the cell population data drops off
    %almost instantaneously for reasons that we don't quite understand.
    %Thus, we'll focus more on the data up through t = 60 until the data is
    %corrected.
    dmx.cellTime = dmx.cellTime(1:35);

    % plot(dmx.cellTime, dmx.N2(2,:)); xlabel('hour'); ylabel('# diploid cells')
    cd(path2root)
    
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
    
    %Returns well number through converting char to int through ASCII
    wellNumber = typecast(N2Well(3), "int16");

    %Thus, subtract 48
    wellNumber = wellNumber - 48;
    
    wellConcentration = concentrationMap(wellNumber);

    dmx.dFdCTP = (wellConcentration / dmx.N2dFdCTP(2)) * dmx.dFdCTP;

    %Disregard P3 Cells until we discuss how to address them - should we
    %just add them to P2?
    dmx.cells = dmx.N2(1:end - 1,ii);

    % need to add zero for cells of highest ploidy, when the data is not included, orelse the model will not penalize their existence
    if (size(dmx.cells, 1) < 5)
        addingZeros = "Adding Row of Zeros because data for highest ploidy is not present"
        dmx.cells = [dmx.cells; zeros(1, size(dmx.cells,2))];
    end

    dmx.cellTime = dmx.cellTime(ii);
   
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % Combined parameters for both models
    theta = 1; nu = 10.0; 
    xi = 7.3125; a = 1.00; 
    v = 381.438; w1 = 20; iota = 0;
    gamma_a = 0.1; lambda = 1;
    k_1 = 0.1925; k_2 = 0.1925;
    K = 2400; mu = 1;
    
    %Initial guess for k_2 comes from a doubling time of 36 hours
    %ln(2) / k_2 = 36 => k_2 = ln(2) / 36 = 0.01925
    %Initial guess for k_1 comes from the fact that the actual ratio b/w
    %G1/S and G2/M is rather large, so k_1 should be significantly larger
    %than k_2 in order to allow for a similar effect in the data.

    pars = {theta, nu, xi, a, v, w1, iota, gamma_a, lambda, k_1, k_2, K, mu};
    bounds = cell2mat(cellfun(@(x) [x / 10; x * 10], pars, 'UniformOutput', false));
    % bounds(:,2) = [pars{2}/10; pars{2}*10];
    % bounds(:,6) = [pars{6}/10; pars{6}*10];
    lb = bounds(1,:)';
    ub = bounds(2,:)';

    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon', 'objective', ...
        @(pars) combined_cost(pars, dmx), 'x0', cell2mat(pars), 'lb', lb, 'ub', ub, 'options', opts);

    rs = RandomStartPointSet('NumStartPoints', 100);
    points = list(rs, problem);
    ms = MultiStart('UseParallel', true);
    [pars_, fval, exitflag, output, solutions] = run(ms, problem, CustomStartPointSet(points));

    S = setfield(S, replicates(k).N2, pars_);
    % save('solutions.mat','pars_');
    % pars=arrayfun(@(x) {x}, pars_);
end


%{
cd(path2PKPDdata)
dmx = struct();
pkpdFile1_4N = dir(['nM1000_4N-', num2str(1), '.txt']);
pkpdFile2_4N = dir(['nM1000_4N-', num2str(2), '.txt']);
pkpdFile3_4N = dir(['nM1000_4N-', num2str(3), '.txt']);

pkpdData1_4N = readtable(pkpdFile1_4N.name);
pkpdData2_4N = readtable(pkpdFile2_4N.name);
pkpdData3_4N = readtable(pkpdFile3_4N.name);

dmx.N4dFdCTP = (pkpdData1_4N.dFdCTP___ng_mL_' + pkpdData2_4N.dFdCTP___ng_mL_' + pkpdData3_4N.dFdCTP___ng_mL_') / 3;
dmx.N4dFdCTP = dmx.N4dFdCTP - dmx.N4dFdCTP(1);

negativedFdCTPTimes = find(dmx.N4dFdCTP <= 0);

dmx.N4dFdCTP(negativedFdCTPTimes) = 0;

dmx.PKPDtime = pkpdData1.time';

for k = 1:1
    % Read PKPD data

    % Read cell data for N4
    cd(path2cellCycledata)
    N4Well = replicates(k).N4
    
    %Returns well number through converting char to int through ASCII
    wellNumber = typecast(N4Well(3), "int16");

    %Thus, subtract 48
    wellNumber = wellNumber - 48;
    
    wellConcentration = concentrationMap(wellNumber);

    dmx.N4 = readtable([N4Well, '.txt']);
    dmx.N4 = table2array(dmx.N4(:, 2:end)); % assuming data starts from the second column
    dmx.cellTime = 2:2:size(dmx.N4, 2) * 2;
    % plot(dmx.cellTime, dmx.N2(2,:)); xlabel('hour'); ylabel('# diploid cells')
    cd(path2root)
    
    %At a point near the end, all of the cell population data drops off
    %almost instantaneously for reasons that we don't quite understand.
    %Thus, we'll focus more on the data up through t = 60 until the data is
    %corrected.
    dmx.cellTime = dmx.cellTime(1:35);

    % Interpolate PKPD data to match cellTime points
    N4_interpolated = interp1(dmx.PKPDtime, dmx.N4dFdCTP, dmx.cellTime, 'linear', 'extrap');
    N4_interpolated(N4_interpolated<0)=0;
   
    % Visualize the interpolation
    figure;
    plot(dmx.PKPDtime, dmx.N4dFdCTP, 'o-', 'DisplayName', 'Original High');
    hold on;
    plot(dmx.cellTime, N4_interpolated, '*-', 'DisplayName', 'Interpolated High');
    title('High Concentration Interpolation');
    xlabel('Time');
    ylabel('Concentration (ng/mL)');
    legend;

    % Match pkpd drug concentration to incucyte drug concentration
    [~,ia]=min(abs(dfdctpConc_incucyte_nM - N4_interpolated));
    N4_interpolated = interp1(dmx.cellTime(1:(1+end-ia)), N4_interpolated(ia:end), dmx.cellTime, 'linear', 'extrap');
    dmx.N4dFdCTP = N4_interpolated;

    %% Do the fitting for one cell type at a time
    ii = 1:length(dmx.cellTime); %find(dmx.cellTime<=55);
    
    %This scales the concentration data to be such that the first
    %timepoint of the interpolated pKpD data, which is then input into 
    % the simulation as an initial condition, is putatively the concentration
    % of the drug that the cells are being exposed to.
    dmx.dFdCTP = dmx.N4dFdCTP(ii);
    dmx.dFdCTP = (wellConcentration / dmx.N4dFdCTP(1)) * dmx.dFdCTP;
       
    % Visualize the interpolation
    figure;
    plot(dmx.cellTime, dmx.dFdCTP, 'o-', 'DisplayName', 'Original High');
    hold on;
    xlabel('Time');
    ylabel('Concentration (ng/mL)');
    legend;
    
    %Disregard P3 Cells until we discuss how to address them - should we
    %just add them to P2?
    dmx.cells = dmx.N4(1:end - 1,ii);


    % need to add zero for cells of highest ploidy, when the data is not included, orelse the model will not penalize their existence
    if (size(dmx.cells, 1) < 5)
        addingZeros = "Adding Row of Zeros because data for highest ploidy is not present"
        dmx.cells = [dmx.cells; zeros(1, size(dmx.cells,2))];
    end

    dmx.cellTime = dmx.cellTime(ii);
   
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % Combined parameters for both models
    theta = 1; nu = 10.0; 
    xi = 7.3125; a = 1.00; 
    v = 381.438; w1 = 20; iota = 0;
    gamma_a = 0.1; lambda = 1;
    k_1 = 0.1925; k_2 = 0.1925;
    K = 2400; mu = 1;
    
    %Initial guess for k_2 comes from a doubling time of 36 hours
    %ln(2) / k_2 = 36 => k_2 = ln(2) / 36 = 0.01925
    %Initial guess for k_1 comes from the fact that the actual ratio b/w
    %G1/S and G2/M is rather large, so k_1 should be significantly larger
    %than k_2 in order to allow for a similar effect in the data.

    pars = {theta, nu, xi, a, v, w1, iota, gamma_a, lambda, k_1, k_2, K, mu};
    bounds = cell2mat(cellfun(@(x) [x / 10; x * 10], pars, 'UniformOutput', false));
    % bounds(:,2) = [pars{2}/10; pars{2}*10];
    % bounds(:,6) = [pars{6}/10; pars{6}*10];
    lb = bounds(1,:)';
    ub = bounds(2,:)';

    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon', 'objective', ...
        @(pars) combined_cost(pars, dmx), 'x0', cell2mat(pars), 'lb', lb, 'ub', ub, 'options', opts);

    rs = RandomStartPointSet('NumStartPoints', 100);
    points = list(rs, problem);
    ms = MultiStart('UseParallel', true);
    [pars_, fval, exitflag, output, solutions] = run(ms, problem, CustomStartPointSet(points));

    S = setfield(S, replicates(k).N4, pars_);
    % save('solutions.mat','pars_');
    % pars=arrayfun(@(x) {x}, pars_);

end
%}

%{
%% Plot Data
figure('name', ['C:\Users\chayn\OneDrive\Documents\USF Moffitt Cancer\AndorLabCode', replicates(1).N4], 'Position', [100 100 1000 400])
combined_cost(getfield(S, replicates(1).N4), dmx)
%}