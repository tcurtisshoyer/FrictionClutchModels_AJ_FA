% This script generates a minimal set of summary figures. It requires a
% SummaryTable that was generated using
% PostProcess_AnalyzeTimescouse_FrictionClutch.m.

close all; clear all; clc;

%% Figure Output Directory
PlotOutDir = 'Example_SummaryOutputs';
if ~exist(PlotOutDir)
    mkdir(PlotOutDir)
end

%% Summary Table Input Directory
SimulationDirectory = 'G:\AJ_FA_ClutchModels_Demo\FrictionClutch_Validation_VelocitySweep';
load(fullfile(SimulationDirectory,['SummaryTable.mat']), 'SummaryTable');
load(fullfile(SimulationDirectory,['AllParamSweep.mat']), 'ParamSweep');

% Load Simulation Data for 5 Linkage Configurations
%   1: Ideal Bond, koff = 10/s
    %   2: Ideal Bond, koff = 5/s
    %   3: Bell Slip Bond, koff0 = 10/s, Fb = 5 pN
    %   4: Bell Slip Bond, koff0 = 5/s, Fb = 5 pN
    %   5: Bell Slip Bond, koff0 = 10/s, Fb = 10 pN  
for linkConfigSweepID = 1:5
    thisSummaryTable = SummaryTable(SummaryTable.linkConfigSweepID==linkConfigSweepID,:);
    Fv_Simulated{linkConfigSweepID} = [thisSummaryTable.V0, thisSummaryTable.mean_totalForce, thisSummaryTable.mean_numEngaged, thisSummaryTable.mean_LinkageEngagedLifetime];
end

%% [Analytical Exp for Sim 1] Ideal Bond Case, Inf Rigid External Spring (Sens, EPL, 2013, Eqn 1)
% Base Params
kon = 2;
koff0 = 10;
Fb = 5;
Nlink = 50;
Klink = 5;
% F = Nlink*Klink*v*kon/(koff0*(kon+koff0))
v_sweep01 = [0 .0001 .001 .01 .1:.1:200]'; % 2x the velocity sweep from FA and AJ friction clutch simulations
F_sweep01 = Nlink*Klink*kon/(koff0*(kon+koff0)) .* v_sweep01;
Fv_Analytical{1} = [v_sweep01, F_sweep01];

%% [Analytical Exp for Sim 2] Ideal Bond Case, Inf Rigid External Spring (Sens, EPL, 2013, Eqn 1)
% Modified Params
kon = 2;
koff0 = 10/2; % 1/2 x koff0
Fb = 5;
Nlink = 50;
Klink = 5;
% F = Nlink*Klink*v*kon/(koff0*(kon+koff0))
v_sweep01 = [0 .0001 .001 .01 .1:.1:200]'; % 2x the velocity sweep from FA and AJ friction clutch simulations
F_sweep01 = Nlink*Klink*kon/(koff0*(kon+koff0)) .* v_sweep01;
Fv_Analytical{2} = [v_sweep01, F_sweep01];

%% [Analytical Exp for Sim 3] Slip Bond Case, Inf Rigid External Spring (Sens, EPL, 2013, Eqns 2-6)
% Base Params
kon = 2;
koff0 = 10;
Fb = 5;
Nlink = 50;
Klink = 5;

v_beta = koff0*Fb/Klink;
r_on = kon/koff0;

% v_star (matches Stochastic Theory) and F_max (overpredicts Stochastic
% Theory), Mean Field Analysis [Sens EPL 2013, Eqn 2-4]
F_max = Nlink * Fb * lambertw(r_on/exp(1));
v_star = v_beta*r_on*(1 + 1/lambertw(r_on/exp(1)) );

% Force-Velocity Relationship, Stochastic Theory [Sens EPL 2013, Eqn 6]
v_sweep02 = [0 .0001 .001 .01 .1:.1:200]'; % 2x the velocity sweep from FA and AJ friction clutch simulations

for ii = 1:length(v_sweep02)
    v = v_sweep02(ii);
    
    v_beta = koff0*Fb/Klink;
    r_on = kon/koff0;
    
    fun1 = @(f) f.*exp(-1.*exp(f)./(v./v_beta));
    q1 = integral(fun1,0,Inf);

    fun2 = @(f) exp(-1.*f)./f;
    q2 = integral(fun2,1./(v./v_beta),Inf);

    Numerator = Nlink.*Fb.*r_on.*exp(1./(v./v_beta)).*q1;
    Denominator = (v./v_beta) + r_on.*exp(1./(v./v_beta)).*q2;

    F_sweep02(ii,:) = Numerator./Denominator;
end

Fv_Analytical{3} = [v_sweep02, F_sweep02];

%% [Analytical Exp for Sim 4] Slip Bond Case, Inf Rigid External Spring (Sens, EPL, 2013, Eqns 2-6)
% Modified Params
kon = 2;
koff0 = 10/2; % 1/2 x koff0
Fb = 5;
Nlink = 50;
Klink = 5;

v_beta = koff0*Fb/Klink;
r_on = kon/koff0;

% v_star (matches Stochastic Theory) and F_max (overpredicts Stochastic
% Theory), Mean Field Analysis [Sens EPL 2013, Eqn 2-4]
F_max = Nlink * Fb * lambertw(r_on/exp(1));
v_star = v_beta*r_on*(1 + 1/lambertw(r_on/exp(1)) );

% Force-Velocity Relationship, Stochastic Theory [Sens EPL 2013, Eqn 6]
v_sweep02 = [0 .0001 .001 .01 .1:.1:200]'; % 2x the velocity sweep from FA and AJ friction clutch simulations

for ii = 1:length(v_sweep02)
    v = v_sweep02(ii);
    
    v_beta = koff0*Fb/Klink;
    r_on = kon/koff0;
    
    fun1 = @(f) f.*exp(-1.*exp(f)./(v./v_beta));
    q1 = integral(fun1,0,Inf);

    fun2 = @(f) exp(-1.*f)./f;
    q2 = integral(fun2,1./(v./v_beta),Inf);

    Numerator = Nlink.*Fb.*r_on.*exp(1./(v./v_beta)).*q1;
    Denominator = (v./v_beta) + r_on.*exp(1./(v./v_beta)).*q2;

    F_sweep02(ii,:) = Numerator./Denominator;
end

Fv_Analytical{4} = [v_sweep02, F_sweep02];

%% [Analytical Exp for Sim 5] Slip Bond Case, Inf Rigid External Spring (Sens, EPL, 2013, Eqns 2-6)
% Modified Params
kon = 2;
koff0 = 10;
Fb = 10; % 2 x Fb
Nlink = 50;
Klink = 5;

v_beta = koff0*Fb/Klink;
r_on = kon/koff0;

% v_star (matches Stochastic Theory) and F_max (overpredicts Stochastic
% Theory), Mean Field Analysis [Sens EPL 2013, Eqn 2-4]
F_max = Nlink * Fb * lambertw(r_on/exp(1));
v_star = v_beta*r_on*(1 + 1/lambertw(r_on/exp(1)) );

% Force-Velocity Relationship, Stochastic Theory [Sens EPL 2013, Eqn 6]
v_sweep02 = [0 .0001 .001 .01 .1:.1:200]'; % 2x the velocity sweep from FA and AJ friction clutch simulations

for ii = 1:length(v_sweep02)
    v = v_sweep02(ii);
    
    v_beta = koff0*Fb/Klink;
    r_on = kon/koff0;
    
    fun1 = @(f) f.*exp(-1.*exp(f)./(v./v_beta));
    q1 = integral(fun1,0,Inf);

    fun2 = @(f) exp(-1.*f)./f;
    q2 = integral(fun2,1./(v./v_beta),Inf);

    Numerator = Nlink.*Fb.*r_on.*exp(1./(v./v_beta)).*q1;
    Denominator = (v./v_beta) + r_on.*exp(1./(v./v_beta)).*q2;

    F_sweep02(ii,:) = Numerator./Denominator;
end

Fv_Analytical{5} = [v_sweep02, F_sweep02];

%% [SUMMARY PLOT] F/v vs v
set(0,'defaultAxesFontSize',8);
set(0,'defaultAxesFontName','Arial');

f = figure;
subplot(1,2,1); hold on;
scatter(Fv_Simulated{1}(:,1),Fv_Simulated{1}(:,2)./Fv_Simulated{1}(:,1), 'ok'); box on;
scatter(Fv_Simulated{2}(:,1),Fv_Simulated{2}(:,2)./Fv_Simulated{2}(:,1), 'sk');
plot(Fv_Analytical{1}(:,1),Fv_Analytical{1}(:,2)./Fv_Analytical{1}(:,1),'-k');
plot(Fv_Analytical{2}(:,1),Fv_Analytical{2}(:,2)./Fv_Analytical{2}(:,1),'-k');
xlabel('v');
ylabel('F/v');
set(gca, 'YScale', 'log'); ylim([10.^([-1 1.5])]);
set(gca, 'XScale', 'log'); xlim([10.^([-1 2])]);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legend('Ideal koff=10 (Sim)',...
    'Ideal koff=5 (Sim)',...
    'Analytical Solns (Sens 2013)',...
    'Location', 'Best');
subplot(1,2,2); hold on;
scatter(Fv_Simulated{3}(:,1),Fv_Simulated{3}(:,2)./Fv_Simulated{3}(:,1), 'or'); box on;
scatter(Fv_Simulated{4}(:,1),Fv_Simulated{4}(:,2)./Fv_Simulated{4}(:,1), 'sr');
scatter(Fv_Simulated{5}(:,1),Fv_Simulated{5}(:,2)./Fv_Simulated{5}(:,1), '^r');
plot(Fv_Analytical{3}(:,1),Fv_Analytical{3}(:,2)./Fv_Analytical{3}(:,1),'-r');
plot(Fv_Analytical{4}(:,1),Fv_Analytical{4}(:,2)./Fv_Analytical{4}(:,1),'-r');
plot(Fv_Analytical{5}(:,1),Fv_Analytical{5}(:,2)./Fv_Analytical{5}(:,1),'-r');
plot(Fv_Analytical{1}(:,1),Fv_Analytical{1}(:,2)./Fv_Analytical{1}(:,1),'--k');
plot(Fv_Analytical{2}(:,1),Fv_Analytical{2}(:,2)./Fv_Analytical{2}(:,1),'--k');
xlabel('v');
ylabel('F/v');
set(gca, 'YScale', 'log'); ylim([10.^([-1 1.5])]);
set(gca, 'XScale', 'log'); xlim([10.^([-1 2])]);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legend('Slip koff0=10, Fb=5 (Sim)',...
    'Slip koff0=5, Fb=5 (Sim)',...
    'Slip koff0=10, Fb=10 (Sim)',...
    'Analytical Solns (Sens 2013)',...
    'Location', 'Best');

% Save full figure panel in specified file formats
figPosition = [100 100 900 300]; f.Position = figPosition;
Resolution_FigurePanel_png = 300; % Figure Print Resolution for PNG

print(f,fullfile(PlotOutDir,['FrictionClutchValidation_EffFricCoeffvsVPlots.png']),'-dpng',['-r' num2str(Resolution_FigurePanel_png)],'-painters');
print(f,fullfile(PlotOutDir,['FrictionClutchValidation_EffFricCoeffvsVPlots.svg']),'-dsvg','-painters');

%% [SUMMARY PLOT] Plot Linkage Engagement Lifetime vs V
set(0,'defaultAxesFontSize',8);
set(0,'defaultAxesFontName','Arial');

f = figure;
subplot(1,2,1); hold on;
scatter(Fv_Simulated{1}(:,1),Fv_Simulated{1}(:,4), 'ok'); box on;
scatter(Fv_Simulated{2}(:,1),Fv_Simulated{2}(:,4), 'sk');
xlabel('v'); set(gca, 'XScale', 'log'); xlim([10.^([-1 2])]);
ylabel('Eng Lifetime [s]'); ylim([0 .25]);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legend('Ideal koff=10 (Sim)',...
    'Ideal koff=5 (Sim)',...
    'Location', 'Best');
subplot(1,2,2); hold on;
scatter(Fv_Simulated{3}(:,1),Fv_Simulated{3}(:,4), 'or'); box on;
scatter(Fv_Simulated{4}(:,1),Fv_Simulated{4}(:,4), 'sr');
scatter(Fv_Simulated{5}(:,1),Fv_Simulated{5}(:,4), '^r');
xlabel('v'); set(gca, 'XScale', 'log'); xlim([10.^([-1 2])]);
ylabel('Eng Lifetime [s]'); ylim([0 .25]);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
legend('Slip koff0=10, Fb=5 (Sim)',...
    'Slip koff0=5, Fb=5 (Sim)',...
    'Slip koff0=10, Fb=10 (Sim)',...
    'Location', 'Best');

% Save full figure panel in specified file formats
figPosition = [100 100 900 300]; f.Position = figPosition;
Resolution_FigurePanel_png = 300; % Figure Print Resolution for PNG
print(f,fullfile(PlotOutDir,['FrictionClutchValidation_EngageLifetimevsVPlots.png']),'-dpng',['-r' num2str(Resolution_FigurePanel_png)],'-painters');
print(f,fullfile(PlotOutDir,['FrictionClutchValidation_EngageLifetimevsVPlots.svg']),'-dsvg','-painters');

%% [TIMECOURSE PLOTS] Plot Selected Timecourse for Showing Representative Behavior
Toggle_TimescourePlots = 1;

% Ideal Bond for V = [1, 20, 100] @ Kext=1E6
if Toggle_TimescourePlots
simIDsToPlot =  [21,51,91]; % V0 = [1, 20, 100]
plotNumSamples = 10000;
plotTimeInterval = [80 90]; % 10 seconds
Toggle_PlotAllNoneEngagedDataPoints = 0;
YLims_Force = [0 100];
YLims_numEngaged = [0 25];
for simID = simIDsToPlot  
    
    load(fullfile(SimulationDirectory,['Data_simID_' num2str(simID) '.mat'])); 
    
    depVarsToPlot = {'Total_Force','numEngaged'};
    depVarsLabels = {'F', 'N'};
    depVarsYLims = {YLims_Force, YLims_numEngaged}; 

    % Obtain timeseries
    id = find(strcmp('t',AllVariableNames.ClutchEnsemble));
    timevals = DataLog_ClutchEnsemble(:,id);
    datavals = [];
    for ii=1:length(depVarsToPlot)
        id = find(strcmp(depVarsToPlot{ii},AllVariableNames.ClutchEnsemble));
        datavals = [datavals, DataLog_ClutchEnsemble(:,id)];
    end

    % Select timespan
    if ~isempty(plotTimeInterval)
        include = timevals>plotTimeInterval(1) & timevals<plotTimeInterval(2);
        timevals = timevals(include);
        datavals = datavals(include,:);
    end

    % Downsample timeseries
    if length(timevals) > plotNumSamples
        interval_downsample = floor(length(timevals)/plotNumSamples);
        index_downsample = [1:interval_downsample:length(timevals)]';
        % Add back in data points where no linkages are bound ("catastrophic
        % failure")
        if Toggle_PlotAllNoneEngagedDataPoints == 1
            numEngaged = datavals(:,3);   
            index_NoneEngaged = find(numEngaged==0);
            if ~isempty(numEngaged)
                index_downsample = [index_downsample; index_NoneEngaged];
                index_downsample = sort(index_downsample);
                index_downsample = unique(index_downsample);
            end
        end
        time = timevals(index_downsample);
        data = datavals(index_downsample,:);
    else
        time = timevals;
        data = datavals;
    end

    % Plot Timecourse
    f = figure;
    f.Position = [100 100 300 150];
    for ii=1:length(depVarsToPlot)
        subplot(length(depVarsToPlot),1,ii);    
        plot(time,data(:,ii),'k');    
        if ~isempty(depVarsLabels{ii})
            ylabel(depVarsLabels{ii});
        end    
        if ~isempty(depVarsYLims{ii})
            if max(data(:,ii))<=depVarsYLims{ii}(2)
                ylim(depVarsYLims{ii});
            end
        end    
        set(gca,'FontSize',6);
        set(gca,'FontName','Arial');
        set(gca,'XMinorTick','on','YMinorTick','on');
        if ii<length(depVarsToPlot)
            xticklabels({});
        end
    end
    subplot(length(depVarsToPlot),1,1);
    title(['SimID = ' num2str(simID)]);
    subplot(length(depVarsToPlot),1,length(depVarsToPlot));
    xlabel('Time [s]');
    
    set(0,'defaultAxesFontSize',8);
    set(0,'defaultAxesFontName','Arial');
    print(gcf,fullfile(PlotOutDir,['TimeCourse_IdealBond_Lowerkoff' '_simID_' num2str(simID) '.svg']),'-dsvg','-painters');
end
end

% Slip Bond for V = [1, 20, 100] @ Kext=1E6
if Toggle_TimescourePlots
simIDsToPlot =  [21,51,91] + 3; % V0 = [1, 20, 100]
plotNumSamples = 10000;
plotTimeInterval = [80 90];
Toggle_PlotAllNoneEngagedDataPoints = 0;
YLims_Force = [0 100];
YLims_numEngaged = [0 25];
for simID = simIDsToPlot   
    
    load(fullfile(SimulationDirectory,['Data_simID_' num2str(simID) '.mat'])); 
    
    depVarsToPlot = {'Total_Force','numEngaged'};
    depVarsLabels = {'F', 'N'};
    depVarsYLims = {YLims_Force, YLims_numEngaged}; 

    % Obtain timeseries
    id = find(strcmp('t',AllVariableNames.ClutchEnsemble));
    timevals = DataLog_ClutchEnsemble(:,id);
    datavals = [];
    for ii=1:length(depVarsToPlot)
        id = find(strcmp(depVarsToPlot{ii},AllVariableNames.ClutchEnsemble));
        datavals = [datavals, DataLog_ClutchEnsemble(:,id)];
    end

    % Select timespan
    if ~isempty(plotTimeInterval)
        include = timevals>plotTimeInterval(1) & timevals<plotTimeInterval(2);
        timevals = timevals(include);
        datavals = datavals(include,:);
    end

    % Downsample timeseries
    if length(timevals) > plotNumSamples
        interval_downsample = floor(length(timevals)/plotNumSamples);
        index_downsample = [1:interval_downsample:length(timevals)]';
        % Add back in data points where no linkages are bound ("catastrophic
        % failure")
        if Toggle_PlotAllNoneEngagedDataPoints == 1
            numEngaged = datavals(:,3);   
            index_NoneEngaged = find(numEngaged==0);
            if ~isempty(numEngaged)
                index_downsample = [index_downsample; index_NoneEngaged];
                index_downsample = sort(index_downsample);
                index_downsample = unique(index_downsample);
            end
        end
        time = timevals(index_downsample);
        data = datavals(index_downsample,:);
    else
        time = timevals;
        data = datavals;
    end

    % Plot Timecourse
    f = figure;
    f.Position = [100 100 300 150];
    for ii=1:length(depVarsToPlot)
        subplot(length(depVarsToPlot),1,ii);    
        plot(time,data(:,ii),'k');    
        if ~isempty(depVarsLabels{ii})
            ylabel(depVarsLabels{ii});
        end    
        if ~isempty(depVarsYLims{ii})
            if max(data(:,ii))<=depVarsYLims{ii}(2)
                ylim(depVarsYLims{ii});
            end
        end    
        set(gca,'FontSize',6);
        set(gca,'FontName','Arial');
        set(gca,'XMinorTick','on','YMinorTick','on');
        if ii<length(depVarsToPlot)
            xticklabels({});
        end
    end
    subplot(length(depVarsToPlot),1,1);
    title(['SimID = ' num2str(simID)]);
    subplot(length(depVarsToPlot),1,length(depVarsToPlot));
    xlabel('Time [s]');
    
    set(0,'defaultAxesFontSize',8);
    set(0,'defaultAxesFontName','Arial');
    print(gcf,fullfile(PlotOutDir,['TimeCourse_SlipBond_Lowerkoff' '_simID_' num2str(simID) '.svg']),'-dsvg','-painters');
end
end