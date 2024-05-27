
%%%%%%%%%%%%%%%%%
%%% PKPD data %%%
%%%%%%%%%%%%%%%%%
 
cd /Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/M00_GemcitabinePKPD_101823
% dm=readtable('dm_2N.txt');

global dmx
global GemcitabineConc_nM
GemcitabineConc_nM=struct(low=100,high=1000);

% dmx=struct(); 
% dmx.high=dm.dFdCTP___ng_mL_';
% dmx.low=dm.dFdCTP___ng_mL__low';
% dmx.time=dm.time';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit model to data %%%
S=struct();
%% Iterate across replicates
for f=dir('nM1000_*txt')';
    dm=readtable(f.name);

    dmx=struct();
    dmx.high=dm.dFdCTP___ng_mL_';   
    dmx.low=dm.dFdCTP___ng_mL__low';
    dmx.time=dm.time';

    %% Do the fitting
    A=[];b= [];Aeq=[];beq=[];

    theta=50; nu=1500; eta=0.005; xi=0.05;
    pars = {theta, nu, eta, xi};
    bounds= cell2mat(cellfun(@(x) [x/5000;x*5000], pars, 'UniformOutput', false));
    lb = bounds(1,:)';
    ub = bounds(2,:)';
    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon','objective',...
        @cost_PKPD,'x0',cell2mat(pars),'lb',lb,'ub',ub,'options',opts);

    rs = RandomStartPointSet('NumStartPoints',25);
    points = list(rs,problem);
    ms = MultiStart('UseParallel',true);
    % [pars_,fval,exitflag,output,solutions]  = run(ms,problem,CustomStartPointSet(points));

    fname=strrep(extractBefore(f.name,12), '-','_');
    % S=setfield(S,fname, pars_);


    %% plot best fit:
    % close all hidden
    figure('name',['~/Downloads/Gemcitabine_PKPD_model_',fname],'Position',[100 100 1000 400])
    cost_PKPD(getfield(S,fname))
end
%% @TODO: decide which parameter values to go with

%% @TODO: normalize all parameters to number of cells (1 Million)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incucyte data: untreated conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% @TODO next


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incucyte data: treated conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Incucyte data: treated conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_SkippedMitosisClassification_042523/');
cd('/Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code')
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/wassersteinFun/
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Data/matlab
replicates=struct('N2',{'B_4'},...
                 'N4',{'E_4'});

global dmx
%global DOSE
DOSE=20;

% Define the initial conditions
P0 = [1; zeros(2,1)]; % P_0 = 1, all other P_i = 0
tspan = [0 5]; % Time interval

% Define the parameters
u = 150;
v = 65;
w1 = 310;
w2=190;
nu=1.095;
iota = 0;

S=struct();

    for k = 1:length({replicates.N2})
        dmx = struct();
        
        for type = {'N2', 'N4'}
            % Load data
            if strcmp(type{1}, 'N4') == 1
                dm = readtable([replicates(k).N4, '.txt']); % 4N
            else
                dm = readtable([replicates(k).N2, '.txt']); % 2N
            end
            
            dm = dm(:, 2:size(dm, 2));
            dmx = setfield(dmx, type{1}, table2array(dm));
        end

        % Solve the problem
        pars = {150, 65, 310, 1.095}; % Example parameters
        bounds = cell2mat(cellfun(@(x) [x / 1000; x * 1500], pars, 'UniformOutput', false));
        lb = bounds(1, :)';
        ub = bounds(2, :)';
        ub(4) = min(2, ub(4));
        lb(4) = max(1, lb(4));
        
        problem = createOptimProblem('fmincon', 'objective',@(pars) cost(pars, dmx,DOSE), 'x0', cell2mat(pars), 'lb', lb, 'ub', ub);
        rs = RandomStartPointSet('NumStartPoints', 250);
        points = list(rs, problem);
        ms = MultiStart('UseParallel', true);
        [pars_, ~, ~, ~, ~] = run(ms, problem, CustomStartPointSet(points));
        
        % Update results
        S = setfield(S, replicates(k).N2, pars_([1, 2, 3, 4]));
    end

    
    
    %% view parameter differences between 2N and 4N
% replicates:
S_=cell2mat(struct2cell(S));
figure('name',['~/Downloads/Gemcitabine_model_params'],'Position',[100 100 600 500])
subplot(2,2,4)
boxplot(S_,{'u','v','w','nu'})
% boxplot(S_(:,2),repmat({'N2','N4'},1,size(S_,1)/2))
% set(gca,'xticklabel',fieldnames(dmx))
set(gca, 'YScale', 'log')
ylabel('value')
subplot(2,2,1)
bar(i,[ alpha_p(i,S.A6(1),0,DOSE);alpha_p(i,S.A6(1),1,DOSE)]','BaseValue',-0.5E-4); %, i,k_p(i,DOSE,S.E2(2),S.E2(4)))
ylim([-0.5E-4,12E-4])
xlabel('number of skipped mitoses')
ylabel('proliferation rate')
legend({'D','T'})
subplot(2,2,3)
bar(i,k_p(i,DOSE,S.A6(2), S.A6(4)),'BaseValue',-0.5E-4); %, i,k_p(i,DOSE,S.E2(2),S.E2(4)))
xlabel('already skipped mitoses')
ylabel('rate of another mitotic slippage')
legend({'D & T'})
subplot(2,2,2)
bar(i, a_p(i,DOSE,S.A6(3),S.A6(3)),'BaseValue',-0.01);%, i,a_p(i,DOSE,S.E2(3),S.E2(3)))
xlabel('number of skipped mitoses')
ylabel('death rate')
legend({'D & T'})
ylim([-0.01,0.12])
% savefigs()


