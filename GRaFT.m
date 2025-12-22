%% GRaFT – Grain-Resolved Fabric and Texture analysis
% =========================================================================
% Description:
%   GRaFT is the analytical engine of the MAPClean ecosystem. It takes 
%   reconstructed grain data (*_finalGrains.mat) produced by GRaMC and 
%   performs automated batch quantification of microstructural fabrics.
%
%   Key Capabilities:
%     1. Texture Strength: J-index (Texture sharpness) & M-index.
%     2. Fabric Shape: Woodcock/Vollmer indices (BA, BC) for shape fabrics.
%     3. Morphometrics: Aspect Ratio maps, Flinn Plots, and CSDs.
%     4. Comparative Analysis: Dataset-wide Kernel Density Estimation (KDE).
%
% Dependencies:
%   - MTEX Toolbox (Tested on v 6.0.0)
%
% Author: Rahul Subbaraman
% Date: December 2025
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% ================= USER CONFIGURATION =================

%% --- 1. Directories ---
% Assumes "Unified Workspace" structure
dataDir   = fullfile(pwd, 'checkpoints');       % Input: Where GRaMC saved final grains
exportDir = fullfile(pwd, 'exports/Textures');  % Output: Analysis results

% Create output directory if it doesn't exist
if ~exist(exportDir, 'dir'), mkdir(exportDir); end

%% --- 2. Stage Control Flags ---
runPlots      = true;   % Generate individual sample plots (PF, IPF, Maps)
runComparison = true;   % Run final comparative analysis (KDE, Contours)

%% --- 3. Global Parameters ---
halfwidth     = 5*degree;          % Halfwidth for ODF calculation
accessibleMap = parula;            % Colormap for Pole Figures
resolution    = 5*degree;          % Resolution for fabric tensor integration
exportRes     = 300;               % DPI for exported images

%% ================= INITIALISATION =================
fprintf('>> GRaFT Analysis Started.\n');

% File Selection: Looks for the specific output format of GRaMC
fileList = dir(fullfile(dataDir, '*_finalGrains.mat'));

if isempty(fileList)
    error('GRaFT Error: No "*_finalGrains.mat" files found in %s. Run GRaMC first.', dataDir);
end

fprintf('   Found %d grain files to process.\n', numel(fileList));

% Initialize Storage Containers
allStats    = table();
allEqDiams  = {};
sampleNames = {};
allGrains   = {}; 

%% ================= LOOP OVER SAMPLES =================
for fi = 1:numel(fileList)
    % --- File Setup ---
    filePath = fullfile(dataDir, fileList(fi).name);
    [~, rawName, ~] = fileparts(fileList(fi).name);
    sampleName = erase(rawName, '_finalGrains'); % Clean sample name
    
    fprintf('\n--- Processing Sample [%d/%d]: %s ---\n', fi, numel(fileList), sampleName);
    
    % Setup Output Directory for this sample
    outDir = fullfile(exportDir, sampleName);
    if ~exist(outDir, 'dir') && runPlots, mkdir(outDir); end

    % --- Load Data ---
    if ~exist(filePath, 'file'), warning('File not found: %s', filePath); continue; end
    S = load(filePath, 'finalGrains');
    if ~isfield(S, 'finalGrains')
        warning('Variable "finalGrains" not found in %s. Skipping.', fileList(fi).name);
        continue;
    end
    grains = S.finalGrains;
    
    % Validation: Check for Orientation Data
    if ~isprop(grains, 'meanOrientation')
        warning('Grains in %s lack meanOrientation data. Skipping.', sampleName);
        continue;
    end
    
    ori = grains.meanOrientation;
    cs  = grains.CS;
    fprintf('  ✔ Data Loaded: %d grains.\n', grains.length);

    % --- Compute Texture Statistics ---
    % Calculate ODF
    odf = calcDensity(ori, 'halfwidth', halfwidth);
    
    % Texture Strength Indices
    J_index = norm(odf);        % J-index (Texture sharpness)
    M_index = calcMIndex(odf);  % M-index (Misorientation index)
    
    % Fabric Indices (Woodcock/Vollmer derived from PF tensor)
    % Defaulting to principal axes: 100, 010, 001
    h = {Miller(1,0,0,cs), Miller(0,1,0,cs), Miller(0,0,1,cs)};
    [pfJ, pfMax] = calcPFStats(odf, h, resolution);
    [BA, BC]     = calcFabricIndices(odf, h, resolution);
    
    fprintf('  ✔ Indices: J=%.2f | M=%.2f | BA=%.2f | BC=%.2f\n', J_index, M_index, BA, BC);

    % --- Compile Statistics Table ---
    T = table();
    T.Sample          = string(sampleName);
    T.NumGrains       = numel(grains);
    T.MeanECD         = mean(2 * sqrt(grains.area / pi));
    T.MeanGrainArea   = mean(grains.area);
    T.MaxGrainArea    = max(grains.area);
    T.MeanAspectRatio = mean(grains.aspectRatio);
    T.MedAspectRatio  = median(grains.aspectRatio);
    T.J_index         = J_index;
    T.M_index         = M_index;
    T.BA_index        = BA;
    T.BC_index        = BC;
    T.MeanGOS         = mean(grains.GOS);
    
    % Append Pole Figure specific stats
    T.pfJ_100 = pfJ(1); T.pfJ_010 = pfJ(2); T.pfJ_001 = pfJ(3);
    T.pfMax_100 = pfMax(1); T.pfMax_010 = pfMax(2); T.pfMax_001 = pfMax(3);
    
    allStats = [allStats; T];

    % --- Store Data for Comparisons ---
    if runComparison
        allEqDiams{fi}  = 2 * sqrt(grains.area / pi); %#ok<SAGROW>
        allGrains{fi}   = grains;                     %#ok<SAGROW>
        sampleNames{fi} = sampleName;                 %#ok<SAGROW>
    end

    % --- Generate Plots ---
    if runPlots
        plotCPOAnalysis(odf, h, cs, outDir, sampleName, accessibleMap, exportRes);
        plotShapeAnalysis(grains, outDir, sampleName, exportRes);
        plotSPOAnalysis(grains, outDir, sampleName, exportRes);
    end
end

%% ================= COMPARATIVE ANALYSIS =================
if runComparison && ~isempty(allEqDiams)
    fprintf('\n--- Generating Comparative Analysis ---\n');
    
    % 1. Export Consolidated Statistics
    outFile = fullfile(exportDir, 'AllSamples_TextureStats.csv');
    writetable(allStats, outFile);
    fprintf('  ✔ Saved Statistics: %s\n', outFile);

    % 2. Grain Size Distribution (KDE Comparison)
    fComp = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 600]);
    hold on;
    colors = lines(numel(allEqDiams));
    
    validSamples = 0;
    for i = 1:numel(allEqDiams)
        if isempty(allEqDiams{i}), continue; end
        [f, x] = ksdensity(log10(allEqDiams{i}));
        plot(10.^x, f, 'LineWidth', 2.5, 'Color', colors(i,:), 'DisplayName', sampleNames{i});
        validSamples = validSamples + 1;
    end
    
    if validSamples > 0
        set(gca, 'XScale', 'log');
        xlabel('Equivalent Diameter (\mum)', 'FontSize', 12);
        ylabel('Probability Density', 'FontSize', 12);
        title('Comparative Grain Size Distributions');
        legend('Interpreter', 'none', 'Location', 'northeast');
        grid on;
        
        exportgraphics(fComp, fullfile(exportDir, 'Combined_CSD_Comparison.png'), 'Resolution', exportRes);
        fprintf('  ✔ Exported: Combined_CSD_Comparison.png\n');
    end
    close(fComp);

    % 3. AR vs Diameter Contours (Combined)
    fAR = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 600]);
    hold on;
    
    validContours = 0;
    for fi = 1:numel(allGrains)
        if isempty(allGrains{fi}), continue; end
        g = allGrains{fi};
        eqDiam = 2 * sqrt(g.area / pi);
        ar = g.aspectRatio;
        
        % Filter invalid data
        valid = eqDiam > 0 & ar > 0 & isfinite(eqDiam) & isfinite(ar);
        if sum(valid) < 10, continue; end % Skip if too few grains
        
        logx = log10(eqDiam(valid)); 
        ar = ar(valid);
        
        % 2D Grid for KDE
        [xg, yg] = meshgrid(linspace(min(logx), max(logx), 100), ...
                            linspace(min(ar), max(ar), 100));
        f_kde = ksdensity([logx(:), ar(:)], [xg(:), yg(:)]);
        F = reshape(f_kde, size(xg));
        
        % Plot contours (Top 3 levels only for clarity)
        maxVal = max(F(:));
        levels = linspace(maxVal*0.3, maxVal*0.9, 3);
        
        for j = 1:3
            [~, hCont] = contour(10.^xg, yg, F, [levels(j) levels(j)], ...
                'LineColor', colors(fi,:), 'LineWidth', j);
            
            % Manage Legend: Only list the thickest contour per sample
            if j == 3
                hCont.DisplayName = sampleNames{fi};
            else
                hCont.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end
        validContours = validContours + 1;
    end
    
    if validContours > 0
        set(gca, 'XScale', 'log');
        xlim([50 1500]); ylim([1 4]);
        xlabel('Equivalent Diameter (\mum)', 'FontSize', 12);
        ylabel('Aspect Ratio', 'FontSize', 12);
        title('Grain Shape vs. Size Trends');
        legend('Location', 'northeastoutside', 'Interpreter', 'none');
        grid on;
        exportgraphics(fAR, fullfile(exportDir, 'Combined_AR_Contours.png'), 'Resolution', exportRes);
        fprintf('  ✔ Exported: Combined_AR_Contours.png\n');
    end
    close(fAR);
end

fprintf('\n>> Analysis Complete.\n');

%% ================= HELPER FUNCTIONS =================

function plotCPOAnalysis(odf, h, cs, outDir, sampleName, cmap, res)
% PLOTCPOANALYSIS Plots Pole Figures (PF) and Inverse Pole Figures (IPF).
    
    % 1. Pole Figures
    fPF = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 400]);
    maxMRD = 0;
    for i = 1:length(h)
        maxMRD = max(maxMRD, max(calcPDF(odf, h{i})));
    end
    maxMRD = ceil(max(maxMRD, 1)); 
    
    for i = 1:length(h)
        nextAxis;
        plotPDF(odf, h{i}, 'antipodal', 'contourf');
        set(gca, 'CLim', [0 maxMRD]);
        mtexTitle(char(h{i}), 'FontSize', 14);
    end
    
    colormap(min(cmap + 0.3, 1)); 
    mtexColorbar('title', 'mrd');
    exportgraphics(fPF, fullfile(outDir, sprintf('%s_01_PF.png', sampleName)), 'Resolution', res);
    close(fPF);

    % 2. Inverse Pole Figures
    dirs = {vector3d.X, vector3d.Y, vector3d.Z};
    dirLabels = {'X', 'Y', 'Z'};
    
    fIPF = figure('Visible', 'off', 'Position', [100 100 1500 500], 'Color', 'w');
    plotIPDF(odf, [dirs{1}, dirs{2}, dirs{3}], 'antipodal', 'contourf');
    colormap(min(cmap + 0.3, 1));
    mtexColorbar;
    
    hLabels = [Miller(1,0,0,cs), Miller(0,1,0,cs), Miller(0,0,1,cs)];
    for i = 1:3
        nextAxis(i);
        title(dirLabels{i}, 'FontSize', 22, 'FontWeight', 'bold');
        hold on;
        plot(hLabels, 'marker', 'o', 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', 'w', 'markersize', 10);
        annotate(hLabels, 'label', {'100', '010', '001'}, 'color', 'k', ...
             'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'none', ...
             'labeloffset', 0.12);
        hold off;
    end
    exportgraphics(fIPF, fullfile(outDir, sprintf('%s_02_IPF.png', sampleName)), 'Resolution', res);
    close(fIPF);
end

function plotShapeAnalysis(grains, exportPath, sName, res)
% PLOTSHAPEANALYSIS Generates Aspect Ratio maps, Flinn plots, and CSDs.

    eqDiam = 2 * sqrt(grains.area / pi);
    
    % 1. Aspect Ratio Map
    fMap = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
    plot(grains, grains.aspectRatio, 'linewidth', 0.1); 
    colormap(parula); mtexColorbar; title('Aspect Ratio Map');
    exportgraphics(fMap, fullfile(exportPath, sprintf('%s_03_AR_Map.png', sName)), 'Resolution', res); 
    close(fMap);

    % 2. Flinn Plot
    fFl = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 600, 600]);
    scatter(norm(grains.shortAxis), norm(grains.longAxis), 50, 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k'); 
    hold on;
    maxVal = max([max(norm(grains.longAxis)), max(norm(grains.shortAxis))]);
    plot([0 maxVal], [0 maxVal], 'r--', 'LineWidth', 2); hold off;
    xlabel('Short Axis Length'); ylabel('Long Axis Length'); 
    axis equal; grid on; title('Flinn Plot');
    exportgraphics(fFl, fullfile(exportPath, sprintf('%s_04_Flinn.png', sName)), 'Resolution', res); 
    close(fFl);
    
    % 3. Grain Size Distribution (KDE)
    fCS = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 800, 500]);
    [f, x] = ksdensity(log10(eqDiam)); 
    plot(10.^x, f, 'LineWidth', 2, 'Color', [0.8 0.4 0.4]);
    set(gca, 'XScale', 'log'); grid on; 
    xlabel('Equivalent Diameter (\mum)'); ylabel('Probability Density');
    title('Grain Size Distribution');
    exportgraphics(fCS, fullfile(exportPath, sprintf('%s_06_CSD_KDE.png', sName)), 'Resolution', res); 
    close(fCS);
end

function plotSPOAnalysis(grains, exportPath, sName, res)
% PLOTSPOANALYSIS Generates a Shape Preferred Orientation (SPO) Rose Diagram.
    fRo = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 600, 600]);
    w = (grains.aspectRatio - 1);
    histogram(grains.longAxis, 50, 'weights', w, 'FaceColor', [0.2 0.4 0.8]);
    ax = gca; 
    if isa(ax, 'PolarAxes') 
        ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; 
        title('Shape Preferred Orientation (SPO)');
    end
    exportgraphics(fRo, fullfile(exportPath, sprintf('%s_05_SPO_Rose.png', sName)), 'Resolution', res); 
    close(fRo);
end

function [pfJ, pfMax] = calcPFStats(odf, hkl, res)
    r = regularS2Grid('resolution', res); w = calcQuadratureWeights(r);             
    pfJ = zeros(numel(hkl),1); pfMax = zeros(numel(hkl),1);
    for i = 1:numel(hkl)
        pf = calcPDF(odf, hkl{i}, r); 
        pfJ(i) = sum(pf.^2 .* w); pfMax(i) = max(pf); 
    end
end

function [BA, BC] = calcFabricIndices(odf, hkl, res)
    lambda = zeros(3, numel(hkl));           
    for i = 1:numel(hkl)
        T_tensor = orientationTensorFromPoleFigure(odf, hkl{i}, res); 
        lambda(:,i) = sort(eig(T_tensor), 'descend'); 
    end
    P = lambda(1,:) - lambda(3,:); G = 2 * (lambda(2,:) - lambda(3,:));       
    BA = 0.5 * (2 - P(2)/(P(2)+G(2)) - P(1)/(P(1)+G(1))); 
    BC = 0.5 * (2 - P(2)/(P(2)+G(2)) - P(3)/(P(3)+G(3))); 
end

function T = orientationTensorFromPoleFigure(odf, hkl, res)
    r = regularS2Grid('resolution', res); w = calcQuadratureWeights(r); 
    pf = calcPDF(odf, hkl, r); T = zeros(3,3); 
    for i = 1:length(r), v = double(r(i)); v = v(:); T = T + pf(i) * (v * v') * w(i); end
    T = T / sum(pf .* w);                         
end
