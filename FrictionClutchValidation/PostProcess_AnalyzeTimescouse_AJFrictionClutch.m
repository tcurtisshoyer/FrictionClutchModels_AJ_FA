% This script perfoms analyses on timecourse data from a set of simulations
% of the AJ friction clutch model. The specified directory
% "SimulationDirectory" contains the data for each simulation in a
% parameter sweep. See ParamSweep_AJFrictionClutch_Velocity.m for an
% example of a parameter sweep.
%
% Refer to Supplementary Note 1 of Shoyer et al. and to the comprehensively
% commented function "runSimulation_AJClutch.m" for complete documentation.
%

close all; clear all; clc; 

%% Input Directory
% The specified directory "SimulationDirectory" contains the data for each
% simulation in a parameter sweep. See
% ParamSweep_AJFrictionClutch_Velocity.m for an example of a parameter
% sweep.
SimulationDirectory = 'G:\AJ_FA_ClutchModels_Demo\FrictionClutch_Validation_VelocitySweep';

%% Analysis Settings
% Time window to analyze DataLog_ClutchEnsemble outputs, specified as [min,
% max] fraction of total time. User must set to properly sample steady
% state.
AnalysisTimeWindow_ClutchEnsemble = [.2 1]; % Time window to analyze DataLog_ClutchEnsemble outputs, specified as [min, max] fraction of total time. User must set to properly sample steady state.

%% Analyze Simulations
load(fullfile(SimulationDirectory,'AllParamSweep.mat'), 'ParamSweep');
fileCheck = dir(fullfile(SimulationDirectory,['SummaryTable.mat']));
if isempty(fileCheck)
for simID = ParamSweep.simID'
    
    load(fullfile(SimulationDirectory,['Data_simID_' num2str(simID) '.mat']));   

    % [1] Time-averaged metrics for Clutch Ensemble (from "DataLog_ClutchEnsemble")
    % Set analysis time window
    idt = find(strcmp('t',AllVariableNames.ClutchEnsemble));
    time = DataLog_ClutchEnsemble(:,idt);
    insideAnalysisTimeWindow = time>=max(time)*AnalysisTimeWindow_ClutchEnsemble(1) &...
        time<=max(time)*AnalysisTimeWindow_ClutchEnsemble(2);
    DataLog_ClutchEnsemble = DataLog_ClutchEnsemble(insideAnalysisTimeWindow,:);
    time = DataLog_ClutchEnsemble(:,idt);
    
    % [1A] # Engaged Linkages and Fraction Engaged Linkages
    idnumEngaged = find(strcmp('numEngaged',AllVariableNames.ClutchEnsemble));
    numEngaged = DataLog_ClutchEnsemble(:,idnumEngaged);
    numTotalLinkages = ParamSweep(simID,:).numTotalLinkages;
    fractionEngaged = numEngaged./numTotalLinkages;
    mean_numEngaged = mean(timeseries(numEngaged,time),'Weighting','time'); % time-weighted mean
    mean_fractionEngaged = mean(timeseries(fractionEngaged,time),'Weighting','time'); % time-weighted mean

    % [1B] Total Force
    idTotal_Force = find(strcmp('Total_Force',AllVariableNames.ClutchEnsemble)); 
    totalForce = DataLog_ClutchEnsemble(:,idTotal_Force);
    mean_totalForce = mean(timeseries(totalForce,time),'Weighting','time'); % time-weighted mean
    
    % [1C] Ensemble Vinculin Molecular Tension
    % Mean force over vinculin in all linkages (including Type 1 linkages where Vcl is not activated)
    % Force over Vcl = Force over Vcl:Factin bond = Force over Vcl:aCat bond.
    idmean_F_Bond_VclFactin_overAllBonds = find(strcmp('mean_F_Bond_VclFactin_overAllBonds',AllVariableNames.ClutchEnsemble));
    mean_F_Bond_VclFactin_overAllBonds = DataLog_ClutchEnsemble(:,idmean_F_Bond_VclFactin_overAllBonds);
    mean_F_Bond_VclFactin_overAllBonds = mean(timeseries(mean_F_Bond_VclFactin_overAllBonds,time),'Weighting','time');
    EnsembleVinculinTension = mean_F_Bond_VclFactin_overAllBonds;
    
    % [2] Mean Linkage Engagement Lifetime (from "DataLog_Transitions")
    % Mean engagement lifetime is computed from "DataLog_Transitions", taking
    % the mean engagement lifetime across all linkage disengagement events.
    % Time window for DataLog_Transitions pre-set as a simulation param. 
    % First, determine disengagement events
    idEngagedBeforeTransition = find(strcmp('EngagedBeforeTransition',AllVariableNames.Transitions));
    idEngagedAfterTransition = find(strcmp('EngagedAfterTransition',AllVariableNames.Transitions));
    DisengageEvents = find(DataLog_Transitions(:,idEngagedBeforeTransition)==1 & DataLog_Transitions(:,idEngagedAfterTransition)==0);
    % Second, compute average engagement lifetime over disengagement events
    idLinkEngagedLifetime = find(strcmp('LinkageEngagedLifetime',AllVariableNames.Transitions));
    LinkEngagedLifetimes = DataLog_Transitions(DisengageEvents,idLinkEngagedLifetime);
    mean_LinkageEngagedLifetime = mean(LinkEngagedLifetimes);
    
    % [3] Compile Summary Stats in Summary Table
    SummaryTable(simID,:) = table(mean_numEngaged,...
        mean_fractionEngaged,...
        mean_totalForce,...
        EnsembleVinculinTension,...
        mean_LinkageEngagedLifetime);
       
end

SummaryTable = [ParamSweep, SummaryTable];
save(fullfile(SimulationDirectory,['SummaryTable.mat']), 'SummaryTable');
writetable(SummaryTable,fullfile(SimulationDirectory,['SummaryTable.xlsx']));

OutDir = 'Example_SummaryOutputs';
if ~exist(OutDir)
    mkdir(OutDir)
end
save(fullfile(OutDir,['SummaryTable.mat']), 'SummaryTable');
writetable(SummaryTable,fullfile(OutDir,['SummaryTable.xlsx']));

end