



cd('/Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Data/M00_GemcitabinePKPD_101823')
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/wassersteinFun/

GemcitabineConc_nM=struct('low',100,'high',1000);


files = dir('nM1000_*txt');
disp({files.name}); % Display filenames


for f = files'
    disp(['Processing file: ', f.name]); % Display current file being processed
    
    % Read data from the file
    dm = readtable(f.name);
    disp(dm(1:5, :)); % Display first few rows of the table
    
    % Extract relevant columns into the structure dmx_PKPD
    dmx_PKPD = struct();
    dmx_PKPD.high = dm.dFdCTP___ng_mL_';
    dmx_PKPD.low = dm.dFdCTP___ng_mL__low';
    dmx_PKPD.time = dm.time';


cd('/Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code')
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Code/wassersteinFun/
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Data/matlab
addpath /Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/Data/M00_GemcitabinePKPD_101823


replicates=struct('N2',{'B_4','B_5','B_6'},...
                 'N4',{'E_4','E_5','E_6'});


% Define the parameters
u = 150;
v = 65;
w1 = 310;
w2=190;
nu=1.095;
iota = 1;

S = struct();
%% Iterate across replicates
%drug=drugModel()

for k = 1:length({replicates.N2})

    dmx_skippedMito = struct();
    for type = {'N2', 'N4'}
        f = dir(".");
        f(1:2) = [];
        f = struct2cell(f);
        f = f(1,:);
        % Use just one Gemcitabine concentration for time being (250) -- column 6!
        f = f(cellfun(@(x) ~isempty(strfind(x, '2')), f) == 1);
        if strcmp(type{1}, 'N4') == 1
            dm = readtable([replicates(k).N4, '.txt']); %% 4N
        else
            dm = readtable([replicates(k).N2, '.txt']); %% 2N
        end
        dm = dm(:, 2:size(dm, 2));
        dmx_skippedMito = setfield(dmx_skippedMito, type{1}, table2array(dm));
    end

    % Do the global fitting for all doses in the set DOSE
    A = []; b = []; Aeq = []; beq = [];

    pars = {u, v, w1, iota,nu,theta,nu_PKPD,eta,xi};
    bounds = cell2mat(cellfun(@(x) [x / 1000; x * 1500], pars, 'UniformOutput', false));
    lb = bounds(1,:)';
    ub = bounds(2,:)';
    ub(4) = min(2, ub(4));
    lb(4) = max(1, lb(4));

    % if any(lb >= ub)
    %         error('Infeasible bounds detected: some lower bounds are >= upper bounds.');
    % end
    
%%    
    
    
     opts = optimoptions(@fmincon);
    problem = createOptimProblem('fmincon','objective',...
        @(pars) combined_cost(pars, dmx_skippedMito,dmx_PKPD, GemcitabineConc_nM),'x0',cell2mat(pars),'lb',lb,'ub',ub,'options',opts);
     rs = RandomStartPointSet('NumStartPoints',250);
     
    points = list(rs,problem);
    ms = MultiStart('UseParallel',true);
   [pars_,fval,exitflag,output,solutions]  = run(ms,problem,CustomStartPointSet(points));

    %opts = optimoptions(@fmincon);

   % Use fmincon for global optimization across all doses (Without using parallel computing)
  %   [pars_, fval, exitflag, output] = fmincon(@(pars) combined_cost(pars, dmx_skippedMito,dmx_PKPD, GemcitabineConc_nM), ...
   %    cell2mat(pars), A, b, Aeq, beq, lb, ub, [], opts);

    % Store global results
    S.(replicates(k).N2).global = pars_([1, 2, 3, 4]);

    % Plot best fit for each dose using the global optimum
    for dose = DOSE
        figure('name', ['~/Downloads/Gemcitabine_model_Dose_', num2str(dose)], 'Position', [100, 100, 1000, 400]);
        cost(pars_, dmx, dose);
    end



   end
end



























