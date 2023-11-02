% This script generates a minimal set of summary figures. It requires a
% SummaryTable that was generated using
% PostProcess_AnalyzeTimescouse_FrictionClutch.m. Also makes trajectories
% of representative timecourses.

close all; clear all; clc;

%% Figure Output Directory
PlotOutDir = 'Example_SummaryOutputs';
if ~exist(PlotOutDir)
    mkdir(PlotOutDir)
end

%% Summary Table Input Directory
SimulationDirectory = 'G:\AJ_FA_ClutchModels_Demo\FA_FrictionClutch_VelocitySweep';
load(fullfile(SimulationDirectory,['SummaryTable.mat']), 'SummaryTable');
load(fullfile(SimulationDirectory,['AllParamSweep.mat']), 'ParamSweep');

%% Figure Settings
xLabel = 'V [nm/s]';
MarkerSize = 8;
LineSpec = '.-';
set(0,'defaultAxesFontSize',9);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesLineWidth',.75);

% Reference Lines for Range of Speeds in um/hr
Toggle_AddRefLines_umperhr = 1;
umhr_per_nms = 0.277778; % Conversion factor (um/hr)/(nm/s)
RefLine_1 = 1*umhr_per_nms;
RefLine_2 = 30*umhr_per_nms;

%% [SUMMARY FIGURE] FA Friction Clutch Velocity Sweep, 5 Diff Vcl Fractions 
% Plot: 
% Linkage Lifetime
% Fraction Linkages Engaged
% Total Force
% Effective Friction Coefficient: Total Force/Speed

f=figure();
% Recall: Linkage Configurations
%   1: 50 Total Linkages, 100% Type 1
%   2: 100 Total Linakges, 100% Type 1
%   3: 50 Total Linkages, 75% Type 1 / 25% Type 3
%   4: 50 Total Linkages, 50% Type 1 / 50% Type 3
%   5: 50 Total Linkages, 25% Type 1 / 75% Type 3
%   6: 50 Total Linkages, 100% Type 3
%   7: 50 Total Linkages, 100% Type 2
legendEntries = {};
for linkConfigSweepID = [1 3 4 5 6] % --> Link Configs [0%,25%,50%,75%,100%] Vcl Reinforced
% for sensitivitySweepID = unique(SummaryTable.sensitivitySweepID')    
    
    thisSummaryTable = SummaryTable(SummaryTable.linkConfigSweepID==linkConfigSweepID,:);
    
    legendEntries = [legendEntries; ['\rho_{Vcl} = ' num2str(thisSummaryTable(1,:).fractionLinkageType3)] ];
    
    V0 = thisSummaryTable.V0;
    
    % Engagement Lifetime
    subplot(2,3,1);
    z = thisSummaryTable.mean_LinkageEngagedLifetime;
    plot(V0,z,LineSpec,'MarkerSize',MarkerSize); hold on;
    set(gca, 'XScale', 'log');
    xlabel(xLabel); xticks(logspace(-1,2,4));
    ylabel('Eng Lifetime [s]');
    xlim([1E-1 1E2]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    ylim([0 .7]);
    yticks([0:.1:.7]);
    
    % Fraction Engaged
    subplot(2,3,2);
    z = thisSummaryTable.mean_fractionEngaged;
    plot(V0,z,LineSpec,'MarkerSize',MarkerSize); hold on;
    set(gca, 'XScale', 'log');
    xlabel(xLabel); xticks(logspace(-1,2,4));
    ylabel('Frac Engaged');
    xlim([1E-1 1E2]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    ylim([0 .7]);
    yticks([0:.1:.7]);
    
    % Total Force
    subplot(2,3,3);
    z = thisSummaryTable.mean_totalForce;
    plot(V0,z,LineSpec,'MarkerSize',MarkerSize); hold on;
    set(gca, 'XScale', 'log');
    xlabel(xLabel);
    ylabel('Force [pN]');
    xlim([1E-1 1E2]); xticks(logspace(-1,2,4));
    set(gca,'XMinorTick','on','YMinorTick','on');
    ylim([0 650]);
    yticks([0:100:650]);
    
    % Effective Fric Coeff = F/V
    subplot(2,3,4);
    z = thisSummaryTable.mean_totalForce./thisSummaryTable.V0;
    plot(V0,z,LineSpec,'MarkerSize',MarkerSize); hold on;
    set(gca, 'XScale', 'log');
    xlabel(xLabel);
    ylabel('F/V [pN*s/nm]');
    xlim([1E-1 1E2]); xticks(logspace(-1,2,4));
    set(gca, 'YScale', 'log');
    set(gca,'XMinorTick','on','YMinorTick','on');
    ylim([10^(-.5) 10^(2.5)]);
    yticks([1E0 1E1 1E2]);  
    
end

if Toggle_AddRefLines_umperhr
    for ii=1:4
        subplot(2,3,ii);
        plot([RefLine_1 RefLine_1],[1E-12 1E12],'Color',[1 1 1]*.4,'LineStyle',':','LineWidth',0.75);
        plot([RefLine_2 RefLine_2],[1E-12 1E12],'Color',[1 1 1]*.4,'LineStyle',':','LineWidth',0.75);
    end
end

% Add Legend
subplot(2,3,6);
plot([0],[0]); hold on;
plot([0],[0]);
plot([0],[0]);
plot([0],[0]);
plot([0],[0]);
legend(legendEntries);

% Figure Panel Display and Print Size ( [x0,y0,w,h] in inches)
FigPanel_x0_inches = 0; % Figure Panel Display and Print Size ( [x0,y0,w,h] in inches)
FigPanel_y0_inches = 0; % Figure Panel Display and Print Size ( [x0,y0,w,h] in inches)
FigPanel_w_inches = 8; % Figure Panel Display and Print Size ( [x0,y0,w,h] in inches)
FigPanel_h_inches = 4; % Figure Panel Display and Print Size ( [x0,y0,w,h] in inches)
Resolution_FigurePanel_png = 300; % Figure Print Resolution for PNG

% Save full figure panel in specified file formats
f.PaperUnits = 'inches';
f.PaperPosition = [FigPanel_x0_inches,FigPanel_y0_inches,FigPanel_w_inches,FigPanel_h_inches];
print(f,fullfile(PlotOutDir,['SummaryPlots_FAFrictionClutch.png']),'-dpng',['-r' num2str(Resolution_FigurePanel_png)],'-painters');
print(f,fullfile(PlotOutDir,['SummaryPlots_FAFrictionClutch.svg']),'-dsvg','-painters');

%% [TIMECOURSE PLOTS 1] Plot Selected Timecourse for Showing Representative Behavior @ V0 = 1 nm/s
Toggle_TimescourePlots = 1;
if Toggle_TimescourePlots
simIDsToPlot =  [29,34]; % Rho_Vcl = 0% (100% Linkage Type 1) and Rho_Vcl = 100% (100% Linkage Type 3)
Labels{29} = {'simID 29: V0=1, \rho_{Vcl}=0'};
Labels{34} = {'simID 34: V0=1, \rho_{Vcl}=1'};
plotNumSamples = 10000;
plotTimeInterval = [900 930]; % Arbitrary 30 second interval
Toggle_PlotAllNoneEngagedDataPoints = 1;
YLims_Force = [0 300];
YLims_numEngaged = [0 35];
for simID = simIDsToPlot  
    
    load(fullfile(SimulationDirectory,['Data_simID_' num2str(simID) '.mat'])); 
    
    depVarsToPlot = {'Total_Force','numEngaged','sum_Bond_VclFactin_StrongBound_overAllBonds'};
    depVarsLabels = {'Total Force [pN]', '# Engaged Links', '# Vcl:Factin in Strong State'};
    depVarsYLims = {YLims_Force, YLims_numEngaged, YLims_numEngaged, YLims_numEngaged};

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
    f.Position = [100 100 400 700];
    for ii=1:length(depVarsToPlot)
        subplot(4,1,ii);    
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
    subplot(4,1,1);
    % title(['SimID = ' num2str(simID)]);
    title(Labels{simID});
    subplot(4,1,length(depVarsToPlot));
    xlabel('Time [s]');
    
    set(0,'defaultAxesFontSize',8);
    set(0,'defaultAxesFontName','Arial');
    print(gcf,fullfile(PlotOutDir,['TimeCourse_simID_' num2str(simID) '.svg']),'-dsvg','-painters');
end
end

%% [TIMECOURSE PLOTS 2] Plot Selected Timecourse for Showing Representative Behavior @ V0 = 1 nm/s
Toggle_TimescourePlots = 1;
if Toggle_TimescourePlots
simIDsToPlot =  [57,62]; % Rho_Vcl = 0% (100% Linkage Type 1) and Rho_Vcl = 100% (100% Linkage Type 3)
Labels{57} = {'simID 57: V0=10, \rho_{Vcl}=0'};
Labels{62} = {'simID 62: V0=10, \rho_{Vcl}=1'};
plotNumSamples = 10000;
plotTimeInterval = [900 930]; % Arbitrary 30 second interval
Toggle_PlotAllNoneEngagedDataPoints = 1;
YLims_Force = [0 1000];
YLims_numEngaged = [0 35];
for simID = simIDsToPlot  
    
    load(fullfile(SimulationDirectory,['Data_simID_' num2str(simID) '.mat'])); 
    
    depVarsToPlot = {'Total_Force','numEngaged','sum_Bond_VclFactin_StrongBound_overAllBonds'};
    depVarsLabels = {'Total Force [pN]', '# Engaged Links', '# Vcl:Factin in Strong State'};
    depVarsYLims = {YLims_Force, YLims_numEngaged, YLims_numEngaged, YLims_numEngaged};

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
    f.Position = [100 100 400 700];
    for ii=1:length(depVarsToPlot)
        subplot(4,1,ii);    
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
    subplot(4,1,1);
    % title(['SimID = ' num2str(simID)]);
    title(Labels{simID});
    subplot(4,1,length(depVarsToPlot));
    xlabel('Time [s]');
    
    set(0,'defaultAxesFontSize',8);
    set(0,'defaultAxesFontName','Arial');
    print(gcf,fullfile(PlotOutDir,['TimeCourse_simID_' num2str(simID) '.svg']),'-dsvg','-painters');
end
end