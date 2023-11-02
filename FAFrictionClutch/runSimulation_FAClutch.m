function [DataLog_ClutchEnsemble,DataLog_SingleLinkage,DataLog_Transitions,AllVariableNames,SimulationPerformanceStats] = runSimulation_FAClutch(ModelParams,SimulationParams)

%% DOCUMENTATION
% runSimulation_FAClutch: A function that simulates the FA friction clutch
% model.
%
% ----- INPUTS -----
%
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
%           1: Intg:FN [binding/unbinding]
%           2: Talin:F-actin [binding/unbinding]
%           3: Vinculin(Vt):F-actin [binding/unbinding] 
%           4: Vinculin(Vh):Talin [binding/unbinding]
%       See "DefineForceDepRateConstantModels_FA.m" for definition of
%       transition rate constant models.
%
%   numLinkagesByType: 1x3 double specifying the number of each type of
%   linkage in the clutch ensemble, i.e. [numLinkagesType1,
%   numLinkagesType2, numLinkagesType3] where:
%           Linkage Type 1 = Non-reinforced, i.e. Tal:Factin is only bond
%               at Factin interface --> FN:Intg-Tal:Factin)
%           Linkage Type 2 = Vinculin Only, i.e. Vcl:Factin is only bond at
%               Factin interface --> FN:Intg-Tal:Vcl:Factin
%           Linkage Type 3 = Reinforced, i.e. both Tal:Factin and
%               Vcl:Factin bonds at Factin interface -->
%               FN:Intg-Tal:(Factin AND Vcl:Factin)
%
%   InitialStates: Structure containing the initial state for each
%       bond/domain for each of the three linkage types. Each of the
%       following fields is a 1x3 double storing the state for
%       LinkageType1, LinkageType2, and LinkageType3, respectively.
%           Bond_IntgFN: State of Intg:FN bond. 0=unbound, 1=bound.
%           Bond_TalFactin: State of Tal:Factin bond. 0=unbound, 1=bound.
%           Bond_VclFactin: State of Vcl:Factin bond. 0=unbound, 1=bound.
%           Bond_VclTal: State of Vcl:Tal bond. 0=unbound, 1=bound.
%
% SimulationParams: A data structure storing params that control how the
% model is simulated and how simulation data is collected; has the
% following fields:
%   tend: total simulation time [s]
%   maxAllowedTimeStep: max allowed timestep [s] in Gillespie SSA. If dt
%       larger than this value, advance time and update mechanics without
%       executing a linkage transition. See comments directly preceding
%       Gillespie SSA implementation for more information.
%   interval_LogLinkageData: time step [s] for logging timecourse data
%       on the clutch ensemble (DataLog_ClutchEnsemble) and on the state of
%       a single linkage (DataLog_SingleLinkage).
%   timeSpan_LogTransitionData: time interval [start end] in [s] for
%       logging linkage state transitions (DataLog_Transitions).
%
%
% ----- OUTPUTS -----
% 
% DataLog_ClutchEnsemble: Timecourse data for the clutch ensemble. Each row
%   is a timepoint. Column headers output in
%   "AllVariableNames.ClutchEnsemble". Defined below where created.
%   All times in [sec], forces in [pN], and extensions in [nm].
%
% DataLog_SingleLinkage: Timecourse data for a single linkage in the clutch
%   (arbitrarily Linkage 1). Each row is a timepoint. Column headers output
%   in "AllVariableNames.SingleLinkage". Defined below where created.
%   All times in [sec], forces in [pN], and extensions in [nm].
% 
% DataLog_Transitions: Data for all linkage state transitions. Each row is
%   a state transition. Column headers output in
%   "AllVariableNames.Transitions". Defined below where created.
%   All times in [sec], forces in [pN], and extensions in [nm].
%
% SimulationPerformanceStats: Performance statistics on the stochastic
%   simulation algorithm (number of accepted transitions, number of
%   rejected transitions, fraction of rejected transitions, and average
%   timestep).
%

%% Setup "DataLog_Transitions"
% DataLog_Transitions: Data for all linkage state transitions. Each row is
%   a state transition. Column headers output in
%   "AllVariableNames.Transitions". Defined below where created.
%   All times in [sec], forces in [pN], and extensions in [nm].
Toggle_LogTransitions = 1; % on/off
initHeight = 1E5;
VariableNames_DataLog_Transitions = {'t',... % time
    'idLinkage',... % ID specifying which linkage undergoes state transition
    'idStateVariable',... % ID specifying which state variable undergoes state transition (Bond_IntgFN = 1, Bond_TalFactin = 2, Bond_VclFactin = 3, Bond_VclTal = 4)
    'StateValueFrom',... % Value of state variable before transition (0=unbound/1=bound for single-bound-state-bonds, 0=unbound/1=weak_bound/2=strong_bound for two-bound-state bonds)
    'StateValueTo',... % Value of state variable after transition (0=unbound/1=bound for single-bound-state-bonds, 0=unbound/1=weak_bound/2=strong_bound for two-bound-state bonds)
    'RateConstant',... % Value of rate constant [1/s] for transition (depends on kinetic model and force)
    'RateConstantForce',... % force across bond/domain used to compute rate constant
    'EngagedBeforeTransition',... % Engaged status of linkage before transition (0=disengaged, 1=engaged)
    'EngagedAfterTransition',... % Engaged status of linkage after transition (0=disengaged, 1=engaged)
    'LinkageEngagedLifetime'}; % Cumalitive lifetime since last transition to become engaged (used to compute engagement lifetime using engaged-->disengaged transitions)
PreAllocationBlocks.DataLog_Transitions = zeros(initHeight,length(VariableNames_DataLog_Transitions));
DataLog_Transitions = PreAllocationBlocks.DataLog_Transitions;
iDataLog_Transitions = 0; % Index of current data log. Incremented in simulation.

%% Setup Ouput "DataLog_SingleLinkage"
% DataLog_SingleLinkage: Timecourse data for a single linkage in the clutch
%   (arbitrarily Linkage 1). Each row is a timepoint. Column headers output
%   in "AllVariableNames.SingleLinkage". Defined here below.
%   All times in [sec], forces in [pN], and extensions in [nm].
initHeight = 1E5;
VariableNames_DataLog_SingleLinkage = {'t',... % time
        'Force',... % force on linkage
        'Extension',... % linkage extension
        'Engaged',... % engaged status of linkage (0=disengaged, 1=engaged)
        'Bond_IntgFN',... % status of Intg:FN bond (0=unbound, 1=bound)
        'Bond_TalFactin',... % status of Tal:Factin bond (0=unbound, 1=bound)
        'Bond_VclFactin',... % status of Vcl:Factin bond (0=unbound, 1=bound weak, 2=bound strong)
        'Bond_VclTal',... % status of Vcl:Tal bond (0=unbound, 1=bound)
        'F_Bond_IntgFN',... % force across Intg:FN bond
        'F_Bond_TalFactin',... % force across Tal:Factin bond
        'F_Bond_VclFactin',... % force across Vcl:Factin bond
        'F_Bond_VclTal',... % force across Vcl:Tal bond
        };
PreAllocationBlocks.DataLog_SingleLinkage = zeros(initHeight,length(VariableNames_DataLog_SingleLinkage));
DataLog_SingleLinkage = PreAllocationBlocks.DataLog_SingleLinkage;

%% Setup Output "DataLog_ClutchEnsemble"
% DataLog_ClutchEnsemble: Timecourse data for the clutch ensemble. Each row
%   is a timepoint. Column headers output in
%   "AllVariableNames.ClutchEnsemble", defined here below.
%   All times in [sec], forces in [pN], and extensions in [nm]
initHeight = 1E5;
VariableNames_DataLog_ClutchEnsemble = {'t',... % time
        'numEngaged',... % # engaged linkages (complete mechanical connection between external spring outside cell and actin inside cell)
        'Actin_Speed',... % speed of actin (loading the linkages); constant for friction clutch and time-varying for motor-clutch
        'Total_Force',... % total force on external spring (equal to force on linkages in series)
        'sum_Link_Force',... % total force across linkage ensemble (equal to force on external spring in series)
        'mean_F_Bond_IntgFN_overAllBonds',... % mean force across all Intg:FN bonds
        'mean_F_Bond_TalFactin_overAllBonds',... % mean force across all Tal:Factin bonds
        'mean_F_Bond_VclFactin_overAllBonds',... % mean force across all Vcl:Factin bonds
        'mean_F_Bond_VclTal_overAllBonds',... % mean force across all Vcl:Tal bonds
        'sum_Bond_IntgFN_Bound_overAllBonds',... % # Intg:FN in bound state (1) across all linkages
        'sum_Bond_TalFactin_Bound_overAllBonds',... % # Tal:Factin in bound state (1) across all linkages
        'sum_Bond_VclFactin_WeakBound_overAllBonds',... % # Vcl:Factin in weak bound state (1) across all linkages
        'sum_Bond_VclFactin_StrongBound_overAllBonds',... % # Vcl:Factin in strong bound state (2) across all linkages
        'sum_Bond_VclTal_Bound_overAllBonds',... % # Vcl:Tal in bound state (1) across all linkages
        };
PreAllocationBlocks.DataLog_ClutchEnsemble = zeros(initHeight,length(VariableNames_DataLog_ClutchEnsemble));
DataLog_ClutchEnsemble = PreAllocationBlocks.DataLog_ClutchEnsemble;
iDataLog_Linkages = 0; % Index of current data log. Incremented in simulation.
t_lastDataLogLinkages = 0;

%% Simulation Parameters
tend = SimulationParams.tend; % total simulation time [s]
interval_LogLinkageData = SimulationParams.interval_LogLinkageData; % time step [s] for logging timecourse data on the clutch ensemble (DataLog_ClutchEnsemble) and on the state of a single linkage (DataLog_SingleLinkage).
timeSpan_LogTransitionData = SimulationParams.timeSpan_LogTransitionData; % time interval [start end] in [s] for logging linkage state transitions (DataLog_Transitions).
maxAllowedTimeStep = SimulationParams.maxAllowedTimeStep; % max allowed timestep [s] in Gillespie SSA. If dt larger than this value, advance time and update mechanics without executing a linkage transition.
rateMax = 1E15; % Max value for force-dependent rate constant evaluation
kT = 4.114; % k_B*T in [pN*nm] at T = 298 K

%% Model Parameters
Kext = ModelParams.springConstExternal; % spring constant [pN/nm] representing the stiffness of the external surface. For FA model, this is the substrate. For AJ model, this is the adjacent cell.
Kc = ModelParams.springConstLinkage; % spring constant [pN/nm] representing the stiffness of a single linkage.
V0 = ModelParams.V0; % speed [nm/s]. For friction clutch, the relative speed between the cell (actin) and the external surface.

% Import Force-Dependent Rate Constant Models
%   RateConstants: Structure array storing the kinetic rate constants for
%   all binding/unbinding state transitions.
%           1: Intg:FN [binding/unbinding]
%           2: Talin:F-actin [binding/unbinding]
%           3: Vinculin(Vt):F-actin [binding/unbinding] 
%           4: Vinculin(Vh):Talin [binding/unbinding]
%       See "DefineForceDepRateConstantModels_FA.m" for definition of
%       transition rate constant models.
RateConstants_Bond_IntgFN = ModelParams.RateConstants{1}; 
RateConstants_Bond_TalFactin = ModelParams.RateConstants{2};
RateConstants_Bond_VclFactin = ModelParams.RateConstants{3};
RateConstants_Bond_VclTal = ModelParams.RateConstants{4};

% Set Number of each type of linkage
%   numLinkagesByType: 1x3 double specifying the number of each type of
%   linkage in the clutch ensemble, i.e. [numLinkagesType1,
%   numLinkagesType2, numLinkagesType3] where:
%           Linkage Type 1 = Non-reinforced, i.e. Tal:Factin is only bond
%               at Factin interface --> FN:Intg-Tal:Factin)
%           Linkage Type 2 = Vinculin Only, i.e. Vcl:Factin is only bond at
%               Factin interface --> FN:Intg-Tal:Vcl:Factin
%           Linkage Type 3 = Reinforced, i.e. both Tal:Factin and
%               Vcl:Factin bonds at Factin interface -->
%               FN:Intg-Tal:(Factin AND Vcl:Factin)
numLinkagesType1 = ModelParams.numLinkagesByType(1);
numLinkagesType2 = ModelParams.numLinkagesByType(2);
numLinkagesType3 = ModelParams.numLinkagesByType(3);
numLinkagesTotal = sum(ModelParams.numLinkagesByType);
Link_Types = [repmat(1,[numLinkagesType1 1]); repmat(2,[numLinkagesType2 1]); repmat(3,[numLinkagesType3 1])];

%% Initial Conditions
% Summary: All linkages are started with all bonds unbound, meaning 0
% linkages are engaged.

% Initialize Clutch Mechanics (forces and extensions)
Link_ExtensionTotal = zeros(numLinkagesTotal,1); % total extension of each linkage, including component from the external spring that is in series with the linkages: x_tot,i = x_link,i + x_ext
Link_ExtensionLinkOnly = zeros(numLinkagesTotal,1); % component of total extension that applies to the linkage only, this is equal to the total extension minus the component of external extension (of the external spring in series): x_link,i
Link_Force = zeros(numLinkagesTotal,1); % force across each linkage: F_link,i = Kc*x_link,i = Kc*(x_tot,i - x_ext)
External_Extension = 0; % extension of external spring in series with the linkages: x_ext
External_Force = 0; % total force on external spring in series with the linkages: F_ext or F_tot

% Initial State Definitions for Each Bond for each Link Type
%   InitialStates: Structure containing the initial state for each
%       bond/domain for each of the three linkage types. Each of the
%       following fields is a 1x3 double storing the state for
%       LinkageType1, LinkageType2, and LinkageType3, respectively.
%           Bond_IntgFN: State of Intg:FN bond. 0=unbound, 1=bound.
%           Bond_TalFactin: State of Tal:Factin bond. 0=unbound, 1=bound.
%           Bond_VclFactin: State of Vcl:Factin bond. 0=unbound, 1=bound.
%           Bond_VclTal: State of Vcl:Tal bond. 0=unbound, 1=bound.
initState_Bond_IntgFN = ModelParams.InitialStates.Bond_IntgFN;
initState_Bond_TalFactin = ModelParams.InitialStates.Bond_TalFactin;
initState_Bond_VclFactin = ModelParams.InitialStates.Bond_VclFactin;
initState_Bond_VclTal = ModelParams.InitialStates.Bond_VclTal;
% Initialize Linkages with each bond in initial condition and 0 force on
% all bonds.
for ii = 1:numLinkagesTotal
    thisLinkType = Link_Types(ii);
    Bond_IntgFN(ii,:) = initState_Bond_IntgFN(thisLinkType); % stores state of Intg:FN bond for each linkage (0=unbound,1=bound)
    Bond_TalFactin(ii,:) = initState_Bond_TalFactin(thisLinkType); % stores state of Tal:Factin bond for each linkage (0=unbound, 1=bound)
    Bond_VclFactin(ii,:) = initState_Bond_VclFactin(thisLinkType); % stores state of Vcl:Factin bond for each linkage (0=unbound, 1=weak_bound, 2=strong_bound)
    Bond_VclTal(ii,:) = initState_Bond_VclTal(thisLinkType); % stores state of Vcl:Tal bond for each linkage (0=unbound,1=bound)
end
F_Bond_IntgFN = zeros(numLinkagesTotal,1); % stores force on Intg:FN bond for each linkage
F_Bond_TalFactin = zeros(numLinkagesTotal,1); % stores force on Tal:Factin bond for each linkage
F_Bond_VclFactin = zeros(numLinkagesTotal,1); % stores force on Vcl:Factin bond for each linakge
F_Bond_VclTal = zeros(numLinkagesTotal,1); % stores force on Vcl:Tal bond for each linkage

% Compute Engaged Status (all linkages should be disengaged to start)
Engaged = Bond_IntgFN & (Bond_TalFactin | (Bond_VclFactin & Bond_VclTal)); % stores engaged status of each linkage (0=disengaged,1=engaged) 

% Structure to track engagement lifetime for each linkage
Link_EngagedLifetimes = zeros(numLinkagesTotal,1); % stores time since last transition to become engaged (used for computing engagement lifetime statistics)

% Initialize Time
t = 0; % time [s], updated in stochastic simulation algorithm below
dt = 0; % time step [s], updated in stochastic simulation algorithm below

%% Stochastic Simulation Algorithm
numAcceptedTransitions = 0; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt, so a reaction will occur along with update to forces
numRejectedTransitions = 0; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt, so no reaction will occur, just update forces
avgTimeStep = 0; % Record simulation performance stats: avg length of time step
while t<tend
    
    %% [1] Determine and Execute Next State Transition and Update Time (Gillespie Direct SSA)
    % Summary of Approach:
    % [1A] Calculate rate constants for each available bond state
    % transition in each linkage based on the kinetic model for the bond
    % and the force across the bond.
    % [1B] Determine and execute earliest state transition and update time
    % using the Direct Gillespie SSA.
    % Note: The stochastic simulation algorithm (SSA) used here is the
    % Gillespie SSA, which was previously used to simulate FA motor-clutch
    % models (e.g. Bangasser et al., Biophys. J., 2013). As linkage forces,
    % and thus force-dependent unbinding rate constants, are time-varying,
    % the Gillespie SSA is an approximation which can be improved by
    % implementing a max allowed timestep. If time-to-next-transition
    % exceeds the max timestep, time is advanced by the max timestep and
    % forces are updated without executing a transition (making use of the
    % memoryless property of Markov processes). For the parameters used in
    % Shoyer et al., it was verified that results with this SSA matched
    % those obtained using a discrete time-step algorithm with a
    % sufficiently small time-step. For new model parameters and
    % configurations, the user is advised to perform analagous
    % verifications.
    
    % [1A] Find all available state transitions (all transitions open to
    % each linkage) and compute each transition rate.
    
    % Possible Transitions is a Nx5 matrix where each row is a possible
    % transition available at the current point in time. Columns store the
    % following:
    % Col 1: Linkage ID
    % Col 2: State Variable ID:
        % Bond_IntgFN = 1
        % Bond_TalFactin = 2
        % Bond_VclFactin = 3
        % Bond_VclTal = 4
    % Col 3: Current/From State Variable Value
    % Col 4: New/To State Variable Value
    % Col 5: Rate Constant
    % Col 6: Force used to compute Rate Constant
    
    allTransitions = zeros(numLinkagesTotal*10,6); % Initialize allTransitions to save time
    count=0;
    for thisLinkID = 1:numLinkagesTotal
        % Check Intg_FN
        thisStateVariableID = 1;
        thesePossibleTransitions = getTransitions_1BoundStateModel_BellwUpperLimit(F_Bond_IntgFN(thisLinkID),Bond_IntgFN(thisLinkID),RateConstants_Bond_IntgFN);
        allTransitions(count+1:count+size(thesePossibleTransitions,1),:) = ...
            [repmat([thisLinkID thisStateVariableID],[size(thesePossibleTransitions,1) 1]) thesePossibleTransitions];
        count = count+size(thesePossibleTransitions,1);
        % Check Talin_Factin (applies only to Linkage Types 1 and 3)
        if Link_Types(thisLinkID) == 1 || Link_Types(thisLinkID) == 3 
            thisStateVariableID = 2;
            thesePossibleTransitions = getTransitions_1BoundStateModel_BellwUpperLimit(F_Bond_TalFactin(thisLinkID),Bond_TalFactin(thisLinkID),RateConstants_Bond_TalFactin);
            allTransitions(count+1:count+size(thesePossibleTransitions,1),:) = ...
                [repmat([thisLinkID thisStateVariableID],[size(thesePossibleTransitions,1) 1]) thesePossibleTransitions];
            count = count+size(thesePossibleTransitions,1);          
        end
        % Check Vcl_Factin and Vcl_Talin (applies only to Linkage Types 2 and 3)
        if Link_Types(thisLinkID) == 2 || Link_Types(thisLinkID) == 3 
            % Check Vcl_Factin
            thisStateVariableID = 3;
            thesePossibleTransitions = getTransitions_2BoundStateModel(F_Bond_VclFactin(thisLinkID),Bond_VclFactin(thisLinkID),RateConstants_Bond_VclFactin);
            allTransitions(count+1:count+size(thesePossibleTransitions,1),:) = ...
                [repmat([thisLinkID thisStateVariableID],[size(thesePossibleTransitions,1) 1]) thesePossibleTransitions];
            count = count+size(thesePossibleTransitions,1);  
            % Check Vcl_Talin
            thisStateVariableID = 4;
            thesePossibleTransitions = getTransitions_1BoundStateModel_BellwUpperLimit(F_Bond_VclTal(thisLinkID),Bond_VclTal(thisLinkID),RateConstants_Bond_VclTal);
            allTransitions(count+1:count+size(thesePossibleTransitions,1),:) = ...
                [repmat([thisLinkID thisStateVariableID],[size(thesePossibleTransitions,1) 1]) thesePossibleTransitions];
            count = count+size(thesePossibleTransitions,1);
        end          
    end   
    allTransitions(allTransitions(:,1)==0,:) = [];
    
    % [1B] Determine and execute earliest state transition and update time
    % using the Direct Gillespie SSA [with maxAllowedTimeStep].
    % Method: 
    % (1) Determine time to next state transition using inverse
    % transform sampling. tau = 1/sum(rateConstants)*-log(r1) wheere r1 is
    % pseudorandom number on [0,1]
    % (2a) If tau <= maxAllowedTimeStep: Determine and execute next state
    % transition based on value of rate constants and using r2, a second
    % pseuodorandom number on [0,1], and advance time by dt=tau. 
    % (2b) Else: Execute no state transition and advance time by
    % dt=maxAllowedTimeStep.
    
    rateConstants = allTransitions(:,5);
    GillespieActivity = sum(rateConstants); % Activity = sum of transition rates
    r1 = rand();
    dt = -log(r1)/GillespieActivity; % Time to next event. This time is exponentiatlly distributed with rate parameter equal to sum of all rates, and it is computed here via inverse transform sampling.
    
    % If time-to-next-transition exceeds maxAllowedTimeStep
    % Then execute no state transition and advance time by maxAllowedTimeStep
    if dt>maxAllowedTimeStep
        dt=maxAllowedTimeStep;
        executeTransition = 0;
        numRejectedTransitions = numRejectedTransitions + 1; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt
    
    % Else (if time-to-next-transition is less than maxAllowedTimeStep)
    % Then execute the state transition and advance time by dt from Gillespie SSA    
    else
        executeTransition = 1;
        numAcceptedTransitions = numAcceptedTransitions + 1; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt
        % Determine next event.
        probvec = cumsum(rateConstants)/GillespieActivity;
        r2 = rand();
        idNextTransition = find(r2<=probvec,1);
        % Execute next transition
        nextTransition = allTransitions(idNextTransition,:);
        thisLinkID = nextTransition(1);
        thisStateVariableID = nextTransition(2);
        thisStateVariableValue = nextTransition(4);
        if thisStateVariableID==1
            Bond_IntgFN(thisLinkID) = thisStateVariableValue;         
        elseif thisStateVariableID==2
            Bond_TalFactin(thisLinkID) = thisStateVariableValue;
        elseif thisStateVariableID==3
            Bond_VclFactin(thisLinkID) = thisStateVariableValue;
        elseif thisStateVariableID==4
            Bond_VclTal(thisLinkID) = thisStateVariableValue;
            % If Vh unbinds Tal, Vt forced to unbind F-actin
            if Bond_VclTal(thisLinkID)==0
                Bond_VclFactin(thisLinkID) = 0;
            end
        end
                  
    end
    
    % Record linkage engagement and force before transition (for use in
    % linkage transition logging)
    EngagedBeforeTransition = Engaged(thisLinkID);
    LinkForceBeforeTransition = Link_Force(thisLinkID);
    EngagedAfterTransition = Bond_IntgFN(thisLinkID) &&(Bond_TalFactin(thisLinkID) || (Bond_VclFactin(thisLinkID)&&Bond_VclTal(thisLinkID)));
    
    % Update record of engaged lifetimes for all linkages
    Link_EngagedLifetimes = Engaged .* (Link_EngagedLifetimes + dt);

    % Update time
    t = t + dt; 
    avgTimeStep = (dt + avgTimeStep*(numRejectedTransitions+numAcceptedTransitions))/(numRejectedTransitions+numAcceptedTransitions+1); % Record simulation performance stats: avg length of time step
    
    % [1C] Optional: Log Linkage Transition Data
    % DataLog_Transitions: Data for all linkage state transitions. Each row is
    %   a state transition. Column headers output in
    %   "AllVariableNames.Transitions", defined above where created.
    %   All times in [sec], forces in [pN], and extensions in [nm].
    if executeTransition && Toggle_LogTransitions && t>timeSpan_LogTransitionData(1) && t<=timeSpan_LogTransitionData(2)
        iDataLog_Transitions = iDataLog_Transitions + 1;
        % Pre-allocate new block if needed
        if iDataLog_Transitions>size(DataLog_Transitions,1)
            DataLog_Transitions = [DataLog_Transitions; PreAllocationBlocks.DataLog_Transitions];
        end
            DataLog_Transitions(iDataLog_Transitions,:) = [t, nextTransition, EngagedBeforeTransition, EngagedAfterTransition, Link_EngagedLifetimes(thisLinkID)];
    end 
    
    % [1D] Update Engagement Status of All Linkages
    Engaged = Bond_IntgFN & (Bond_TalFactin | (Bond_VclFactin & Bond_VclTal));
    
    %% [2] Update Clutch Mechanics (Forces and Extensions)
    % Note 1: Clutch mechanics are handled as described previously in
    % Bangasser et al, Biophys. J., 2013, except in the Friction Clutch
    % the actin speed is set to an imposed value V0 corresponding to the
    % relative speed between the cell (actin) and external surface.
    
    % Note 2: Clutch mechanics are updated after linkage transition to ensure
    % that mechanical equilibrium is re-established and linkage mechanics
    % are updated before logging timecourse data.
    
    % Variable Definitions:
    % Va = actin speed
    % Link_ExtensionTotal = total extension of each linkage, including component from the external spring that is in series with the linkages: x_tot,i = x_link,i + x_ext
    % Link_ExtensionLinkOnly = component of total extension that applies to the linkage only, this is equal to the total extension minus the component of external extension (of the external spring in series): x_link,i
    % Link_Force = force across each linkage: F_link,i = Kc*x_link,i = Kc*(x_tot,i - x_ext)
    % External_Extension = extension of external spring in series with the linkages: x_ext
    % External_Force = total force on external spring in series with the linkages (i.e. Traction Force in previous FA Motor-Clutchs): F_ext or F_tot

    % [2A] Actin speed
    % Imposed/constant Va=V0 corresponding to relative speed between cell
    % (actin cytoskeleton) and the external surface (adjacent cell for AJ
    % and substrate for FA).
    Va = V0; % Note: Mathematically equivalent to Va = V0*(1-Ftot/(Nm*Fm)) in limit that (Nm*Fm)-->Infinity
    
    % [2B] Update total extension of each ENGAGED linkage
    Link_ExtensionTotal = Link_ExtensionTotal + Engaged*Va*dt;

    % [2C] Calculate extension of external spring
    % Based on a force balance across linkages and external system:
    % Ftot = Flinks = Fext
    % Flinks = sum(Force on engaged links)
    % = Kc*sum(Link_ExtensionTotal - External_Extension)
    % = Kc*sum(Link_ExtensionTotal(Engaged)) - Kc*numEngaged*External_Extension
    % Fext = Kext * External_Extension
    % Solve for External_Extension: 
    % External_Extension = Kc*sum(Link_ExtensionTotal)/(Kext + Kc*sum(Engaged))
    External_Extension = Kc*sum(Link_ExtensionTotal(Engaged))/(Kext + Kc*sum(Engaged));
    External_Force = Kext*External_Extension;

    % [2D] Update total extension of each DISENGAGED linkage
    % For disengaged linkages, simply equal to extension of external spring
    Link_ExtensionTotal(~Engaged) = External_Extension;

    % [2E] Calculate force across each linkage
    % Flink = Kc*(Link_ExtensionTotal - External_Extension)
    Link_Force = Kc*(Link_ExtensionTotal - External_Extension);
    Link_ExtensionLinkOnly = Link_Force/Kc; % component of total extension that applies to the linkage (Link_ExtensionTotal - External_Extension)
   
    % [2F] Determine force across each bond in each linkage
    for thisLinkID = 1:numLinkagesTotal
    % If linkage is not engaged, then all forces are 0    
    if ~Engaged(thisLinkID)
        F_Bond_IntgFN(thisLinkID) = 0;
        F_Bond_TalFactin(thisLinkID) = 0;
        F_Bond_VclFactin(thisLinkID) = 0;
        F_Bond_VclTal(thisLinkID) = 0;
        
    % Else (if linkage is engaged)
    else
        % Intg_FN always takes full force
        F_Bond_IntgFN(thisLinkID) = Link_Force(thisLinkID);
        
        % Force across adapters depends on which bonds are bound
        % Because this else statement deals only with engaged linkages,
        % this means that either (1) Tal only, (2) Vcl only, or (3) both
        % Tal and Vcl form force-bearing link to F-actin. For Tal to
        % form force-bearing link, must have
        % Bond_TalFactin(thisLinkID)~=0. For Vcl to form force-bearing
        % link, must have (Bond_VclFactin(thisLinkID)~=0 &&
        % Bond_VclTal(thisLinkID)~=0).
        % Case 1 for Adapters: Tal is only force-bearing link to F-actin
        if ~(Bond_VclFactin(thisLinkID) && Bond_VclTal(thisLinkID)) % Already confirmed whole linkage engaged, so this means Tal is only force-bearing link to F-actin
            F_Bond_TalFactin(thisLinkID) = Link_Force(thisLinkID);
            F_Bond_VclFactin(thisLinkID) = 0;
            F_Bond_VclTal(thisLinkID) = 0;
        % Case 2 for Adapters: Vcl is only force-bearing link to F-actin
        elseif ~Bond_TalFactin(thisLinkID) % Already confirmed linkage engaged, so this means Vcl is only force-bearing link to F-actin
            F_Bond_TalFactin(thisLinkID) = 0;
            F_Bond_VclFactin(thisLinkID) = Link_Force(thisLinkID);
            F_Bond_VclTal(thisLinkID) = Link_Force(thisLinkID);   
        % Case 3 for Adapters: Both Tal and Vcl form force-bearing links to F-actin
        % Then: Force split equally across two bonds
        else % Already confirmed link engaged, so this means both Tal and Vcl form force-bearing links to F-actin
            F_Bond_TalFactin(thisLinkID) = Link_Force(thisLinkID)/2;
            F_Bond_VclFactin(thisLinkID) = Link_Force(thisLinkID)/2;
            F_Bond_VclTal(thisLinkID) = Link_Force(thisLinkID)/2;   
        end
    end
    end
                                        
    %% [3] Record Timecourse Data   
    
    % DataLog_ClutchEnsemble: Timecourse data for the clutch ensemble. Each row
    %   is a timepoint. Column headers output in
    %   "AllVariableNames.ClutchEnsemble", defined above where created.
    %   All times in [sec], forces in [pN], and extensions in [nm].
    %
    % DataLog_SingleLinkage: Timecourse data for a single linkage in the clutch
    %   (arbitrarily Linkage 1). Each row is a timepoint. Column headers output
    %   in "AllVariableNames.SingleLinkage", defined above where created.
    %   All times in [sec], forces in [pN], and extensions in [nm].
    %

    % Log linkage data at regular time intervals
    if (t-t_lastDataLogLinkages) >= interval_LogLinkageData
        iDataLog_Linkages = iDataLog_Linkages + 1;
        t_lastDataLogLinkages = t;
        % Pre-allocate new block if needed
        if iDataLog_Linkages>size(DataLog_SingleLinkage,1)
            DataLog_SingleLinkage = [DataLog_SingleLinkage; PreAllocationBlocks.DataLog_SingleLinkage];
            DataLog_ClutchEnsemble = [DataLog_ClutchEnsemble; PreAllocationBlocks.DataLog_ClutchEnsemble];
        end
        
        % DataLog_SingleLinkage: Timecourse data for a single linkage in the clutch
        %   (arbitrarily Linkage 1). Each row is a timepoint. Column headers output
        %   in "AllVariableNames.SingleLinkage", defined above where created.
        %   All times in [sec], forces in [pN], and extensions in [nm].
        DataLog_SingleLinkage(iDataLog_Linkages,:) = [...
                    t,...
                    Link_Force(1),...
                    Link_ExtensionLinkOnly(1),...
                    Engaged(1),...
                    Bond_IntgFN(1),...
                    Bond_TalFactin(1),...
                    Bond_VclFactin(1),...
                    Bond_VclTal(1),...
                    F_Bond_IntgFN(1),...
                    F_Bond_TalFactin(1),...
                    F_Bond_VclFactin(1),...
                    F_Bond_VclTal(1)];      
        % DataLog_ClutchEnsemble: Timecourse data for the clutch ensemble. Each row
        %   is a timepoint. Column headers output in
        %   "AllVariableNames.ClutchEnsemble", defined above where created.
        %   All times in [sec], forces in [pN], and extensions in [nm].
        DataLog_ClutchEnsemble(iDataLog_Linkages,:) = [...
                    t,...
                    sum(Engaged),...
                    Va,...
                    External_Force,...
                    sum(Link_Force),...
                    mean(F_Bond_IntgFN),...
                    mean(F_Bond_TalFactin),...
                    mean(F_Bond_VclFactin),...
                    mean(F_Bond_VclTal),...
                    sum(Bond_IntgFN==1),...
                    sum(Bond_TalFactin==1),...
                    sum(Bond_VclFactin==1),...
                    sum(Bond_VclFactin==2),...
                    sum(Bond_VclTal==1),...
                    ];
    end
            
end
    
%% Clean-up Data Structures for Export
% Remove unused preallocated storage
DataLog_Transitions(iDataLog_Transitions+1:end,:) = [];
DataLog_SingleLinkage(iDataLog_Linkages+1:end,:) = [];
DataLog_ClutchEnsemble(iDataLog_Linkages+1:end,:) = [];
 
% Package up variable names for export
AllVariableNames.Transitions = VariableNames_DataLog_Transitions;
AllVariableNames.SingleLinkage = VariableNames_DataLog_SingleLinkage;
AllVariableNames.ClutchEnsemble = VariableNames_DataLog_ClutchEnsemble;

% Compute Simulation Performance Stats
SimulationPerformanceStats.numAcceptedTransitions = numAcceptedTransitions; % Record simulation performance stats: # time steps where Gillespie waiting time is less than max allowed dt, so a reaction will occur along with update to forces
SimulationPerformanceStats.numRejectedTransitions = numRejectedTransitions; % Record simulation performance stats: # time steps where Gillespie waiting time is greater than max allowed dt, so no reaction will occur, just update forces
SimulationPerformanceStats.probabilityRejectTransition = numRejectedTransitions/(numAcceptedTransitions+numRejectedTransitions);
SimulationPerformanceStats.avgTimeStep = avgTimeStep; % Record simulation performance stats: avg length of time step

%% XXXXX Support Functions XXXXX
% (Below)

%% One-pathway Bell Model
% A function to get the force-dependent rate constant according to a
% one-pathway bell model: k = k01*exp(x1*f). Bell model parameters are
% input as a 2 element row vector: [k01 x1].
function k = kbell(f,BellParams)
%     k01 = BellParams(1);
%     x1 = BellParams(2); 
    k = BellParams(1)*exp(BellParams(2)*f/kT);
    k = min(k,rateMax);
end

%% Two-pathway Bell Model
% A function to get the force-dependent rate constant according to a
% two-pathway bell model: k = k01*exp(x1*f) + k02*exp(x2*f). Bell model
% parameters are input as a 4 element row vector: [k01 x1 k02 x2].
function k = kbell2(f,BellParams)
%     k01 = BellParams(1);
%     x1 = BellParams(2); 
%     k02 = BellParams(3); 
%     x2 = BellParams(4); 
    k = BellParams(1)*exp(BellParams(2)*f/kT) + BellParams(3)*exp(BellParams(4)*f/kT);
    k = min(k,rateMax);
end

%% Functions to Determine Possible State Transitions and Associated Rates
% (1) Bonds w/ 1-Bound-State
% A function to get all possible transitions based on the current state of
% a bond, in a one bound state system using Bell Model (1 or 2 Pathway).
% Outputs a matrix where each row corresponds to a possible transition:
% [fromState toState rateConst]. Unbinding rate constant follows 1 or 2
% Pathway Bell Model
function possibleTransitions = getTransitions_1BoundStateModel_Bell(F,currentState,Params)

    if currentState==0 % Unbound State
        rate_0to1 = Params.kon;
        possibleTransitions = [0 1 rate_0to1 F];
    elseif currentState==1 % Bound State
        rate_1to0 = kbell2(F, [Params.k01 Params.x1 Params.k02 Params.x2]);
        possibleTransitions = [1 0 rate_1to0 F];
    end
    
end

% (2) Bonds w/ 1-Bound-State Bell Model with Upper Lifetime Limit.
% A function to get all possible transitions based on the current state of
% a bond, in a one bound state system using Bell Model (1 or 2 Pathway)
% with Upper Lifetime Limit. Outputs a matrix where each row corresponds to
% a possible transition: [fromState toState rateConst] Unbinding rate
% constant follows 1 or 2 Pathway Bell Model.
function possibleTransitions = getTransitions_1BoundStateModel_BellwUpperLimit(F,currentState,Params)

    if currentState==0 % Unbound State
        rate_0to1 = Params.kon;
        possibleTransitions = [0 1 rate_0to1 F];
    elseif currentState==1 % Bound State
        rate_1to0 = max(1/Params.maxLifetime, kbell2(F, [Params.k01 Params.x1 Params.k02 Params.x2]));
        possibleTransitions = [1 0 rate_1to0 F];
    end
    
end

% (3) Bonds w/ 2-Bound-State Model
% A function to get all possible transitions based on the current state of
% a bond, in the two bound state model. Outputs a matrix where each row
% corresponds to a possible transition: [fromState toState rateConst]
% Rate constants follow Ideal Bond or 1 Pathway Bell Models
function possibleTransitions = getTransitions_2BoundStateModel(F,currentState,Params)

    if currentState==0 % Unbound State
        rate_0to1 = Params.k01;
        rate_0to2 = Params.k02;
        possibleTransitions = [0 1 rate_0to1 F; 0 2 rate_0to2 F];
    elseif currentState==1 % Weak Bound State
        rate_1to0 = kbell(F, [Params.k0_10 Params.x10]);
        rate_1to2 = kbell(F, [Params.k0_12 Params.x12]); 
        possibleTransitions = [1 0 rate_1to0 F; 1 2 rate_1to2 F];
    elseif currentState==2 % Strong Bound State
        rate_2to0 = kbell(F, [Params.k0_20 Params.x20]);
        rate_2to1 = kbell(F, [Params.k0_21 Params.x21]); 
        possibleTransitions = [2 0 rate_2to0 F; 2 1 rate_2to1 F];
    end
    
end

% (4) Bonds w/ 1-Bound-State Dudko-Hummer-Szabo Model Slip Bond
% A function to get all possible transitions based on the current state of
% a bond, in any one bound state model. Outputs a matrix where each row
% corresponds to a possible transition: [fromState toState rateConst].
% Unbinding rate constant follows Dudko-Hummer-Szabo Model Slip Bond 
function possibleTransitions = getTransitions_1BoundStateModel_Dudko(F,currentState,Params)

    v = Params.nu; % scaling factor that specifies the free-energy profile. nu=0.5 for harmonic energy well.
    x = Params.x; % [nm] distance between the bound state and the transition state
    DeltaG_0 = Params.DeltaG_ActivationNoForce; % [kT] free energy of activation in absence of external force
    k0 = Params.k0; % intrinsic off rate
    
    koff = k0*(1-v*F*x/DeltaG_0)^(1./v-1)*exp(DeltaG_0/kT*(1-(1-v*F*x/DeltaG_0)^(1/v)));
    if currentState==0 % Unbound State
        rate_0to1 = Params.kon;
        possibleTransitions = [0 1 rate_0to1 F];
    elseif currentState==1 % Bound State
        rate_1to0 = koff;
        possibleTransitions = [1 0 rate_1to0 F];
    end
    
end

end