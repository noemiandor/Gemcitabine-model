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
% cd('/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_SkippedMitosisClassification_042523/');
cd('/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K01_WGDClassification_050124/');
addpath /Users/4470246/Repositories/Gemcitabine-model/Code/
addpath /Users/4470246/Repositories/Gemcitabine-model/Code/wassersteinFun/
addpath  /Users/4470246/Documents/Matlab-workspace/SelectionForces_GastricCLs/Code/utils
addpath  /Users/4470246/Documents/Matlab-workspace/SelectionForces_GastricCLs/Code/utils/export_fig/
replicates=struct('N2',{'A6','B6','C6','D6'},'N4',{'E6','F6','G6','H6'});

global dmx
global DOSE
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

% figure('name',['~/Downloads/Gemcitabine_model_functions'],'Position',[100 100 900 200])
% % Test k_p function
% subplot(1,3,1)
% % k_p= @(i,x,v) repmat(v*x,length(i),1); %% constant
% % k_p = @(i,x,v) [(v*x).^(i+1)]; %% decreasing
% k_p = @(i,x,v, nu)  (nu-(1-i).^2).*(x/(x+v));
% i=0:(length(P0) - 1);
% params=[v,v*2];
% plot(i, k_p(i,DOSE,params(1), nu), i,k_p(i,DOSE,params(2), nu));%,i,k_p(i,DOSE,params(3), nu))
% xlabel('number of already skipped mitoses')
% ylabel('rate of another mitotic slippage')
% legend({'D','T'})
% % legend(strcat('v=',num2str(params')))
% 
% % Test a_p function
% subplot(1,3,2)
% % a_p= @(i,x,w1,w2) (w1.^(1./(i+1).^2)* x./(x+w2));
% % a_p= @(i,x,w1,w2) w1*x*(3-i).^2;
% a_p= @(i,x,w1,w2) ((1/(length(P0)-1))*(length(P0)-1-i)).^2 .*(x/(x+w1));
% i=0:(length(P0) - 1);
% params=[w1/2,w1,w1*2];
% plot(i, a_p(i,DOSE,params(1),params(1)),'black');%, i,a_p(i,DOSE,params(2)),i,a_p(i,DOSE,params(3)))
% plot(i, a_p(i,DOSE,params(1),w2), i,a_p(i,DOSE,params(2),w2),i,a_p(i,DOSE,params(3),w2))
% xlabel('number of skipped mitoses')
% ylabel('death rate')
% legend({'D and T'})
% % legend(strcat('w=',num2str(params')))
% 
% % Test alpha_p function
% subplot(1,3,3)
% alpha_p= @(i,u,iota, x)  (x./(x+u)) .*(i<=iota);
% % alpha_p= @(i,u,iota, x)  (1-x./(x+u)) .*(i<=iota); %% goodness of fit is same but converges slower -- @TODO: matters once we fit to >1 drug concentrations
% i=0:(length(P0) - 1);
% plot(i, alpha_p(i,u,0,DOSE), i,alpha_p(i,u,1,DOSE));
% ylim([0,u*1.3])
% xlabel('number of skipped mitoses')
% ylabel('proliferation rate')
% legend({'D','T'})
% % savefigs()
% 
% % Solve the ODEs
% [t,y] = ode45(@(t,y) skippedMito_ODE(t,y,Inf,v,w1,  iota, nu), tspan, [0; P0]);
% 
% % Plot the results
% subplot(1,3,3)
% plot(t, y)
% xlabel('Time')
% ylabel('P')
% legend('Dead','P_0','P_1','P_2','P_3','P_4','P_5','P_6','P_7','P_8','P_9')
% % set(gca, 'YScale', 'log')




%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit model to data %%%
S=struct();
%% Iterate across replicates
for k=1:length({replicates.N2})

    dmx=struct();
    for type={'N2','N4'}
        f=dir(".");
        f(1:2)=[];
        f=struct2cell(f);
        f=f(1,:);
        % use just one Gemcitabine concentration for time being (250) -- column 6!
        f=f(cellfun(@(x) ~isempty(strfind(x,'2')), f)==1);
        if strcmp(type{1},'N4')==1
            dm=readtable([replicates(k).N4,'.txt']); %% 4N
        else
            dm=readtable([replicates(k).N2,'.txt']); %% 2N
        end
        dm=dm(:,2:size(dm,2));
        % % use only the first 5 days: @TODO remove and model drug decay
        % ii=find((1:size(dm,2))*4/24 <3);
        % dm=dm(1:(1+length(P0)),ii);
        dmx=setfield(dmx,type{1},table2array(dm));
    end
    % %% Correct image based classification: @TODO -- remove <- should be done by HALO + r
    % dmx.N4(2:3,:)=dmx.N4(1:2,:);
    % dmx.N4(1,:)=0;

    %% @TODO: read in dead cell counts and model dead cell compartment as well <-- continue here


    %% Do the fitting
    A=[];b= [];Aeq=[];beq=[];

    pars = {u,v, w1, nu};
    bounds= cell2mat(cellfun(@(x) [x/1000;x*1500], pars, 'UniformOutput', false));
    lb = bounds(1,:)';
    ub = bounds(2,:)';
    ub(4)=min(2,ub(4));
    lb(4)=max(1,lb(4));

    opts = optimoptions(@fmincon);

    problem = createOptimProblem('fmincon','objective',...
        @cost,'x0',cell2mat(pars),'lb',lb,'ub',ub,'options',opts);
    rs = RandomStartPointSet('NumStartPoints',250);
    points = list(rs,problem);
    ms = MultiStart('UseParallel',true);
    [pars_,fval,exitflag,output,solutions]  = run(ms,problem,CustomStartPointSet(points));
    S=setfield(S,replicates(k).N2,pars_([1,2, 3,4]));
    % S=setfield(S,replicates(k).N4,pars_([1,3,2]));


    %% plot best fit:
    close all hidden
    figure('name','~/Downloads/Gemcitabine_model','Position',[100 100 1000 400])
    cost(pars_)
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
