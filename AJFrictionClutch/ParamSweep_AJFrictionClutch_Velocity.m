% This script perfoms a set of simulations with the AJ Friction Clutch
% model in which the parameter speed V0 is swept and for each value of V0
% the following linkage configurations are simulated:
%   1: 50 Total Linkages, 100% Type 1
%   2: 100 Total Linakges, 100% Type 1
%   3: 50 Total Linkages, 75% Type 1 / 25% Type 3
%   4: 50 Total Linkages, 50% Type 1 / 50% Type 3
%   5: 50 Total Linkages, 25% Type 1 / 75% Type 3
%   6: 50 Total Linkages, 100% Type 3
%   7: 50 Total Linkages, 100% Type 2
%
% Refer to Supplementary Note 1 of Shoyer et al. and to the comprehensively
% commented function "runSimulation_AJClutch.m" for complete documentation.
%

close all; clear all; clc; 

%% Setup Output Directory
ParentDir = 'G:\AJ_FA_ClutchModels_Demo';

ModelType = 'AJ_FrictionClutch';
StudyType = 'VelocitySweep';

ChildDir = [ModelType '_' StudyType];
OutDir = fullfile(ParentDir,ChildDir);
if ~exist(OutDir)
    mkdir(OutDir)
end

%% Simulation Params (see runSimulation_AJClutch.m for complete details)
% SimulationParams: A data structure storing params that control how the
% model is simulated and how simulation data is collected; has the
% following fields:
%   tend: total simulation time [s]
%   maxAllowedTimeStep: max allowed timestep [s] in Gillespie SSA. If dt
%       larger than this value, advance time and update mechanics without
%       executing a linkage transition.
%   interval_LogLinkageData: time step [s] for logging timecourse data
%       on the clutch ensemble (DataLog_ClutchEnsemble) and on the state of
%       a single linkage (DataLog_SingleLinkage).
%   timeSpan_LogTransitionData: time interval [start end] in [s] for
%       logging linkage state transitions (DataLog_Transitions).
%
% Read documentation in "runSimulation_AJClutch.m" to properly set these
% simulation parameters. For a given parameter combination, the user should
% verify that the results with the SSA match those obtained using a
% discrete time-step algorithm with a sufficiently small time-step.
SimulationParams.tend = 1000;
SimulationParams.interval_LogLinkageData = 1E-3; 
SimulationParams.maxAllowedTimeStep = .01;
SimulationParams.timeSpan_LogTransitionData = [.5 1]*SimulationParams.tend;

%% Base Model Parameters (see runSimulation_AJClutch.m for complete details)
% ModelParams: A data structure storing params for the model; has the
% following fields:
%
%   V0: speed [nm/s]. For friction clutch, the relative speed between
%       the cell (actin cytoskeleton) and the external surface (adjacent
%       cell for AJ and substrate for FA).
%
%   springConstExternal: Spring constant [pN/nm] representing the stiffness
%       of the external surface. For FA model, this is the substrate. For
%       AJ model, this is the adjacent cell.
%
%   springConstLinkage: Spring constant [pN/nm] representing the stiffness
%   of a single linkage.
%
%   RateConstants: Structure array storing the kinetic rate constants for
%   all binding/unbinding state transitions.
%           1: Cadherin:Cadherin [binding/unbinding]
%           2: a-Catenin:F-actin [binding/unbinding]
%           3: Vinculin(Vt):F-actin [binding/unbinding]
%           4: Vinculin(Vh):a-Catenin [binding/unbinding]
%       See "DefineForceDepRateConstantModels.m" for definition of
%       transition rate constant models.
%
%   numLinkagesByType: 1x3 double specifying the number of each type of
%   linkage in the clutch ensemble, i.e. [numLinkagesType1,
%   numLinkagesType2, numLinkagesType3] where:
%           Linkage Type 1 = Non-reinforced, i.e. aCat:Factin is only bond
%               at Factin interface --> Cad:Cad-aCat:Factin)
%           Linkage Type 2 = Vinculin Only, i.e. Vcl:Factin is only bond at
%               Factin interface --> Cad:Cad-aCat:Vcl:Factin
%           Linkage Type 3 = Reinforced, i.e. both aCat:Factin and
%               Vcl:Factin bonds at Factin interface -->
%               Cad:Cad-aCat:(Factin AND Vcl:Factin)
%
%   InitialStates: Structure containing the initial state for each
%       bond/domain for each of the three linkage types. Each of the
%       following fields is a 1x3 double storing the state for
%       LinkageType1, LinkageType2, and LinkageType3, respectively.
%           Bond_CadCad: State of Cad:Cad bond. 0=unbound, 1=bound.
%           Bond_aCatFactin: State of aCat:Factin bond. 0=unbound, 1=bound.
%           Bond_VclFactin: State of Vcl:Factin bond. 0=unbound, 1=bound.
%           Bond_VclaCat: State of Vcl:aCat bond. 0=unbound, 1=bound.

% Mechanics Parameters
BaseParams.springConstExternal = [10^(0.5)]; % Kext [pN/nm] {Corresponds to stiffness of cell monolayer 20-33 kPa}
BaseParams.numTotalLinkages = [50];
BaseParams.springConstLinkage = [5];
BaseParams.V0 = [10*1E3/60/60]; % V0 [nm/s] {10 um/hr ~ typical cell speed for collectively migrating MDCK cells --> Estimate for cell-cell relative motion}

% Unbinding RateConstants
% Import Force-Dependent Rate Constant Models
%   RateConstants: Structure array storing the kinetic rate constants for
%   all binding/unbinding state transitions.
%           1: Cadherin:Cadherin [binding/unbinding]
%           2: a-Catenin:F-actin [binding/unbinding]
%           3: Vinculin(Vt):F-actin [binding/unbinding]
%           4: Vinculin(Vh):a-Catenin [binding/unbinding]
%       See "DefineForceDepRateConstantModels.m" for definition of
%       transition rate constant models.
DefineForceDepRateConstantModels_AJ;
BaseParams.RateConstants{1} = Ecad_Ecad_WT_DudkoHummerSzaboParams;
BaseParams.RateConstants{2} = aCatFactin_TwoBoundStateParams;
BaseParams.RateConstants{3} = Vcl_Factin_NegPull_TwoBoundStateParams;
BaseParams.RateConstants{4} = Vcl_aCat_BellModelWithMaxLifetimeLimitParams;

% Binding Rate Constants
BaseParams.RateConstants{1}.kon = [2]; % on rate for Ecad_Ecad. [1/s]
BaseParams.RateConstants{2}.k01 = [2]; % on rate for aCat_Factin. [1/s]
BaseParams.RateConstants{2}.k02 = [0]; % No transition from unbound to Strong Bound State in 2-Bound-State Model
BaseParams.RateConstants{3}.k01 = [2]; % on rate for Vcl_Factin. [1/s]
BaseParams.RateConstants{3}.k02 = [0]; % No transition from unbound to Strong Bound State in 2-Bound-State Model
BaseParams.RateConstants{4}.kon = [1]; % on rate for Vcl_aCat. [1/s]
   
% Set Number of each type of linkage
%   numLinkagesByType: 1x3 double specifying the number of each type of
%   linkage in the clutch ensemble, i.e. [numLinkagesType1,
%   numLinkagesType2, numLinkagesType3] where:
%           Linkage Type 1 = Non-reinforced, i.e. aCat:Factin is only bond
%               at Factin interface --> Cad:Cad-aCat:Factin)
%           Linkage Type 2 = Vinculin Only, i.e. Vcl:Factin is only bond at
%               Factin interface --> Cad:Cad-aCat:Vcl:Factin
%           Linkage Type 3 = Reinforced, i.e. both aCat:Factin and
%               Vcl:Factin bonds at Factin interface -->
%               Cad:Cad-aCat:(Factin AND Vcl:Factin)
BaseParams.fractionLinkageType1 = [1]; % fraction of linkages that are Type 1
BaseParams.fractionLinkageType2 = [0]; % fraction of linkages that are Type 2
BaseParams.fractionLinkageType3 = [0]; % fraction of linkages that are Type 3
BaseParams.numLinkagesByType(1) = round(BaseParams.numTotalLinkages*BaseParams.fractionLinkageType1);
BaseParams.numLinkagesByType(2) = round(BaseParams.numTotalLinkages*BaseParams.fractionLinkageType2);
BaseParams.numLinkagesByType(3) = round(BaseParams.numTotalLinkages*BaseParams.fractionLinkageType3); 

%% Initial Conditions
% Summary: All linkages are started with all bonds unbound, meaning 0
% linkages are engaged.
%
% Initial State Definitions for Each Link Type
% State for [LinkType1 LinkType2 LinkType3]
%   InitialStates: Structure containing the initial state for each
%       bond/domain for each of the three linkage types. Each of the
%       following fields is a 1x3 double storing the state for
%       LinkageType1, LinkageType2, and LinkageType3, respectively.
%           Bond_CadCad: State of Cad:Cad bond. 0=unbound, 1=bound.
%           Bond_aCatFactin: State of aCat:Factin bond. 0=unbound, 1=bound.
%           Bond_VclFactin: State of Vcl:Factin bond. 0=unbound, 1=bound.
%           Bond_VclaCat: State of Vcl:aCat bond. 0=unbound, 1=bound.
BaseParams.InitialStates.Bond_CadCad = [0 0 0];
BaseParams.InitialStates.Bond_aCatFactin = [0 0 0];
BaseParams.InitialStates.Bond_VclFactin = [0 0 0];
BaseParams.InitialStates.Bond_VclaCat = [0 0 0];

%% Parameter Sweep
V0_Sweep = [logspace(-1,2,13)]; % Velocity Sweep
baseSimID = 0; % Each V0 value has a unique baseSimID
simID = 0; % Each simulation (V0 value + linkage configuration) has a unique simID
for V0 = V0_Sweep % Velocity Sweep

    baseSimID = baseSimID + 1;
    
    % Set Model Params to Base Params
    ModelParams = BaseParams;
    
    % Set V0 to current sweep value
    ModelParams.V0 = V0;
    
    % For each value of above defined parameter sweep, simulate multiple
    % different linkage configurations:
    %   1: 50 Total Linkages, 100% Type 1
    %   2: 100 Total Linakges, 100% Type 1
    %   3: 50 Total Linkages, 75% Type 1 / 25% Type 3
    %   4: 50 Total Linkages, 50% Type 1 / 50% Type 3
    %   5: 50 Total Linkages, 25% Type 1 / 75% Type 3
    %   6: 50 Total Linkages, 100% Type 3
    %   7: 50 Total Linkages, 100% Type 2
    LinkConfigSweep_numTotalLinkages = [50 100 50 50 50 50 50];
    LinkConfigSweep_fractionLinkageType1 = [1 1 .75 .5 .25 0 0];
    LinkConfigSweep_fractionLinkageType2 = [0 0 0 0 0 0 1];
    LinkConfigSweep_fractionLinkageType3 = [0 0 .25 .5 .75 1 0];
    LinkConfigSweep_singleBondOnRateConstant_Vcl_aCat_Sweep = [0 0 1 1 1 1 1];
    
    for linkConfigSweepID = 1:length(LinkConfigSweep_numTotalLinkages)
        
        simID = simID + 1;
        
        % Setup linkage configuration as defined above.
        ModelParams.numTotalLinkages = LinkConfigSweep_numTotalLinkages(linkConfigSweepID);
        ModelParams.fractionLinkageType1 = LinkConfigSweep_fractionLinkageType1(linkConfigSweepID);
        ModelParams.fractionLinkageType2 = LinkConfigSweep_fractionLinkageType2(linkConfigSweepID);
        ModelParams.fractionLinkageType3 = LinkConfigSweep_fractionLinkageType3(linkConfigSweepID);
        ModelParams.singleBondOnRateConstant_Vcl_aCat = LinkConfigSweep_singleBondOnRateConstant_Vcl_aCat_Sweep(linkConfigSweepID);
        ModelParams.numLinkagesByType(1) = round(ModelParams.numTotalLinkages*ModelParams.fractionLinkageType1);
        ModelParams.numLinkagesByType(2) = round(ModelParams.numTotalLinkages*ModelParams.fractionLinkageType2);
        ModelParams.numLinkagesByType(3) = round(ModelParams.numTotalLinkages*ModelParams.fractionLinkageType3); 
        
        % Run Simulation
        fileCheck = dir(fullfile(OutDir,['Data_simID_' num2str(simID) '.mat']));
        if isempty(fileCheck)
            tic;
            [DataLog_ClutchEnsemble,DataLog_SingleLinkage,DataLog_Transitions,AllVariableNames,SimulationPerformanceStats] = runSimulation_AJClutch(ModelParams,SimulationParams);
            t=toc;
            disp(['simID ' num2str(simID) ' complete in ' num2str(t) ' sec.']);
            % Save simulation data
            SimulationPerformanceStats.computerTime = t;
            save(fullfile(OutDir,['Data_simID_' num2str(simID) '.mat']), 'DataLog_ClutchEnsemble', 'DataLog_SingleLinkage', 'DataLog_Transitions', 'AllVariableNames', 'SimulationPerformanceStats', 'ModelParams', 'SimulationParams');   
        end
        
        % Make record of parameters
        numTotalLinkages = ModelParams.numTotalLinkages;
        fractionLinkageType1 = ModelParams.fractionLinkageType1;
        fractionLinkageType2 = ModelParams.fractionLinkageType2;
        fractionLinkageType3 = ModelParams.fractionLinkageType3;
        % V0 set by for loop. V0 = ModelParams.V0;
        springConstExternal = ModelParams.springConstExternal;
        springConstLinkage = ModelParams.springConstLinkage;
       
        ParamSweep(simID,:) = table(simID, baseSimID, linkConfigSweepID,...
            numTotalLinkages,...
            fractionLinkageType1,...
            fractionLinkageType2,...
            fractionLinkageType3,...
            V0,...
            springConstExternal,...
            springConstLinkage);
        
    end
    
end

save(fullfile(OutDir,'AllParamSweep.mat'), 'ParamSweep');
