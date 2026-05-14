%% GRaFT – Grain-Resolved Fabric and Texture Analysis
% =========================================================================
% Process-agnostic characterisation workflow
%
% Design split:
%   - MATLAB / MTEX = maps, EBSD-derived calculations, checkpointing, and CSV export
%   - Python        = non-map statistical plots, cross-sample comparison, and interpretation
%
% Companion analysis:
%   This MATLAB script is designed to be used with the accompanying GRaFT
%   Python notebook. The MATLAB stage exports grain, contact, pixel, and
%   summary CSV files. The Python notebook reads these exported tables and
%   performs downstream statistical visualisation, population comparison, and
%   interpretative plotting.
%
% Phase logic:
%   - Anorthite  -> core texture calculations + maps 01–11 + grain CSV
%   - Forsterite -> core texture calculations + maps 01–11 + grain CSV
%   - Diopside   -> core texture calculations + maps 01–11 + grain CSV
%
% Maps produced here:
%   Core maps (all phases):
%     01 PF
%     02 IPF
%     03 Aspect ratio
%     04 Size (ECD)
%     05 Long-axis orientation
%     06 GOS
%     07 GROD
%     08 GROD/GOS ratio
%     09 Crystal shape
%     10 KAM
%     11 Ellipse fit
%
% CSVs exported:
%   All phases:
%     *_grains.csv
%
%   Phase-conditional:
%     *_contacts.csv     (Anorthite only)
%     *_pixels.csv       (Anorthite only)
%
%   Summary:
%     AllSamples_TextureStats.csv
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode', 'silent');
warning('off', 'all');

%% ================= USER CONFIGURATION =================
dataDir       = fullfile(pwd, 'checkpoints');
exportDir     = fullfile(pwd, 'GRaFT');
csvDir        = fullfile(exportDir, 'exports');
checkpointDir = fullfile(pwd, 'checkpoints');

if ~exist(exportDir, 'dir'), mkdir(exportDir); end
if ~exist(csvDir, 'dir'), mkdir(csvDir); end
if ~exist(checkpointDir, 'dir'), mkdir(checkpointDir); end

%% ================= STAGE FLAGS =================
run_statsCalc        = true;
run_samplePlots      = true;
run_exportCSVs       = true;
run_contactAnalysis  = true;
run_pixelExport      = true;

%% ================= PARAMETERS =================
global params
params.overwritePlots      = true;
params.halfwidth           = 5 * degree;
params.resolution          = 5 * degree;
params.exportRes           = 300;
params.mapLineWidth        = 0.5;
params.fontName            = 'Arial';
params.axisFontSize        = 11;
params.vesicleMinAreaPx    = 100;
params.contactDilateSize   = 3;
params.minGrainsForODF     = 100;
params.saveIntermediateMAT = true;

%% ================= LOOP OVER SAMPLE-PHASE FILES =================
fileList = dir(fullfile(dataDir, '*_*_finalGrains.mat'));
summaryRows = table();

if isempty(fileList)
    error('GRaFT Error: No *_*_finalGrains.mat files found in %s.', dataDir);
end

for fi = 1:numel(fileList)
    tSample = tic;

    filePath = fullfile(dataDir, fileList(fi).name);
    [~, rawName, ~] = fileparts(fileList(fi).name);
    baseName = erase(rawName, '_finalGrains');

    parts = split(baseName, '_');
    if numel(parts) < 2
        fprintf('Skipping unrecognised file name: %s\n', fileList(fi).name);
        continue;
    end

    sampleName = string(parts{1});
    phaseName  = string(strjoin(parts(2:end), '_'));
    isAnorthite = strcmpi(phaseName, 'Anorthite');

    phaseExportDir = fullfile(exportDir, char(sampleName), char(phaseName));
    if ~exist(phaseExportDir, 'dir'), mkdir(phaseExportDir); end

    diaryFile = fullfile(phaseExportDir, sprintf('%s_%s_logfile.txt', sampleName, phaseName));
    diary off;
    if exist(diaryFile, 'file'), delete(diaryFile); end
    diary(diaryFile); diary on;

    fprintf('\n====================================================\n');
    fprintf('Processing Sample: %s | Phase: %s\n', sampleName, phaseName);
    fprintf('====================================================\n');

    graftStatsFile = fullfile(checkpointDir, sprintf('%s_%s_GRaFTstats.mat', sampleName, phaseName));
    contactMatFile = fullfile(checkpointDir, sprintf('%s_%s_ContactData.mat', sampleName, phaseName));
    vesicleMatFile = fullfile(checkpointDir, sprintf('%s_%s_VesicleData.mat', sampleName, phaseName));

    grainCSVFile   = fullfile(csvDir, sprintf('%s_%s_grains.csv', sampleName, phaseName));
    contactCSVFile = fullfile(csvDir, sprintf('%s_%s_contacts.csv', sampleName, phaseName));
    pixelCSVFile   = fullfile(csvDir, sprintf('%s_%s_pixels.csv', sampleName, phaseName));
    vesRelCSVFile  = fullfile(csvDir, sprintf('%s_%s_vesicle_grain_relationships.csv', sampleName, phaseName));

    %% -------------------------------------------------
    %% STEP 1: CORE SAMPLE STATISTICS
    %% -------------------------------------------------
    tStep = tic;
    fprintf('[STEP 1] Core sample statistics\n');

    doCalc = run_statsCalc || ~exist(graftStatsFile, 'file');

    if doCalc
        fprintf('  -> Calculating / refreshing checkpoint\n');

        S_load = load(filePath, 'finalGrains', 'ebsd_final');
        grains = S_load.finalGrains;
        ebsd   = S_load.ebsd_final;

        ori = grains.meanOrientation;
        cs  = grains.CS;
        useODF = grains.length >= params.minGrainsForODF;
        h = {Miller(1,0,0,cs), Miller(0,1,0,cs), Miller(0,0,1,cs)};

        if useODF
            odf = calcDensity(ori, 'halfwidth', params.halfwidth);
            J_index = norm(odf);
            M_index = calcMIndex(odf);
            [pfJ, pfMax] = calcPFStats(odf, h, params.resolution);
            [BA, BC] = calcFabricIndices(odf, h, params.resolution);
        else
            odf = [];
            J_index = NaN;
            M_index = NaN;
            pfJ = nan(3,1);
            pfMax = nan(3,1);
            BA = NaN;
            BC = NaN;
        end

        eqDiam = 2 * sqrt(grains.area / pi);
        longAxisAngle = atan2(grains.longAxis.y, grains.longAxis.x) / degree;
        longAxisAngle(longAxisAngle < 0) = longAxisAngle(longAxisAngle < 0) + 180;

        T = table();
        T.Sample = string(sampleName);
        T.Phase = string(phaseName);
        T.NGrains = grains.length;
        T.J_index = J_index;
        T.M_index = M_index;
        T.BA_index = BA;
        T.BC_index = BC;
        T.pfJ_100 = pfJ(1);
        T.pfJ_010 = pfJ(2);
        T.pfJ_001 = pfJ(3);
        T.pfMax_100 = pfMax(1);
        T.pfMax_010 = pfMax(2);
        T.pfMax_001 = pfMax(3);
        T.UseODF = useODF;

        save(graftStatsFile, 'sampleName', 'phaseName', 'grains', 'ebsd', 'odf', 'useODF', ...
            'eqDiam', 'h', 'cs', 'longAxisAngle', 'J_index', 'M_index', 'BA', 'BC', 'T');

        fprintf('  -> Saved: %s\n', graftStatsFile);
    else
        fprintf('  -> Loading existing checkpoint\n');
    end

    S = load(graftStatsFile);
    fprintf('✔ STEP 1 complete (%.2f s)\n', toc(tStep));

    %% -------------------------------------------------
    %% STEP 2: MAPS
    %% -------------------------------------------------
    tStep = tic;
    fprintf('[STEP 2] Maps\n');

    if run_samplePlots
        prefix = sprintf('%s_%s', char(S.sampleName), char(S.phaseName));

        filePF      = fullfile(phaseExportDir, [prefix '_01_PF.png']);
        fileIPF     = fullfile(phaseExportDir, [prefix '_02_IPF.png']);
        fileAR      = fullfile(phaseExportDir, [prefix '_03_AR_Map.png']);
        fileSize    = fullfile(phaseExportDir, [prefix '_04_Size_Map.png']);
        fileOri     = fullfile(phaseExportDir, [prefix '_05_LongAxis_Map.png']);
        fileGOS     = fullfile(phaseExportDir, [prefix '_06_GOS_Map.png']);
        fileGROD    = fullfile(phaseExportDir, [prefix '_07_GROD_Map.png']);
        fileHet     = fullfile(phaseExportDir, [prefix '_08_GROD_GOS_Ratio_Map.png']);
        fileCS      = fullfile(phaseExportDir, [prefix '_09_CrystalShape_Map.png']);
        fileKAM     = fullfile(phaseExportDir, [prefix '_10_KAM_Map.png']);
        fileEllipse = fullfile(phaseExportDir, [prefix '_11_EllipseFit_Map.png']);

        if params.overwritePlots || ~exist(filePF, 'file') || ~exist(fileIPF, 'file')
            tSub = tic;
            plotCPOAnalysis(S.odf, S.h, S.cs, S.grains, S.useODF, phaseExportDir, prefix, params);
            fprintf('  -> PF + IPF exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> PF + IPF exist: skipped\n');
        end

        if params.overwritePlots || ~exist(fileAR, 'file')
            tSub = tic;
            fAR = figure('Visible', 'off', 'Color', 'w');
            plot(S.grains, S.grains.aspectRatio, 'linewidth', params.mapLineWidth, 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'linewidth', params.mapLineWidth, 'color', 'k'); hold off;
            colormap(plasma); mtexColorbar;
            set(gca, 'CLim', [1, 4]);
            exportgraphics(fAR, fileAR, 'Resolution', params.exportRes);
            close(fAR);
            fprintf('  -> Aspect Ratio map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> Aspect Ratio map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileSize, 'file')
            tSub = tic;
            fSize = figure('Visible', 'off', 'Color', 'w');
            plot(S.grains, S.eqDiam, 'linewidth', params.mapLineWidth, 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'linewidth', params.mapLineWidth, 'color', 'k'); hold off;
            set(gca, 'ColorScale', 'log'); colormap(viridis);
            mtexColorbar;
            exportgraphics(fSize, fileSize, 'Resolution', params.exportRes);
            close(fSize);
            fprintf('  -> Size map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> Size map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileOri, 'file')
            tSub = tic;
            fOriMap = figure('Visible', 'off', 'Color', 'w');
            plot(S.grains, S.longAxisAngle, 'linewidth', params.mapLineWidth, 'micronbar', 'off');
            hold on;
            plot(S.grains.boundary, 'linewidth', params.mapLineWidth, 'color', 'k');
            quiver(S.grains, S.grains.longAxis, 'Color', 'black');
            hold off;
            colormap(plasma); set(gca, 'CLim', [0 180]);
            mtexColorbar;
            exportgraphics(fOriMap, fileOri, 'Resolution', params.exportRes);
            close(fOriMap);
            fprintf('  -> Long-axis map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> Long-axis map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileGOS, 'file')
            tSub = tic;
            fGOS = figure('Visible', 'off', 'Color', 'w');
            plot(S.grains, S.grains.GOS / degree, 'linewidth', params.mapLineWidth, 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'linewidth', params.mapLineWidth, 'color', 'k'); hold off;
            colormap(plasma); mtexColorbar;
            set(gca, 'CLim', [0, 5]);
            exportgraphics(fGOS, fileGOS, 'Resolution', params.exportRes);
            close(fGOS);
            fprintf('  -> GOS map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> GOS map exists: skipped\n');
        end

        grod_obj = S.ebsd.calcGROD(S.grains);

        if params.overwritePlots || ~exist(fileGROD, 'file')
            tSub = tic;
            fGROD = figure('Visible', 'off', 'Color', 'w');
            plot(S.ebsd, grod_obj.angle / degree, 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'linewidth', 1.0, 'color', 'k'); hold off;
            colormap(viridis); mtexColorbar;
            set(gca, 'CLim', [0, 10]);
            exportgraphics(fGROD, fileGROD, 'Resolution', params.exportRes);
            close(fGROD);
            fprintf('  -> GROD map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> GROD map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileHet, 'file')
            tSub = tic;
            fHet = figure('Visible', 'off', 'Color', 'w');
            grod_ang = grod_obj.angle / degree;
            isIndexed = S.ebsd.grainId > 0;
            gos_pixel = nan(size(grod_ang));
            gos_pixel(isIndexed) = S.grains(S.ebsd.grainId(isIndexed)).GOS / degree;
            rel_strain = grod_ang ./ max(gos_pixel, 0.05);
            plot(S.ebsd(isIndexed), rel_strain(isIndexed), 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'linewidth', 1.0, 'color', 'k'); hold off;
            colormap(viridis); mtexColorbar;
            set(gca, 'CLim', [0, 2.5]);
            exportgraphics(fHet, fileHet, 'Resolution', params.exportRes);
            close(fHet);
            fprintf('  -> GROD/GOS ratio map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> GROD/GOS ratio map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileCS, 'file')
            tSub = tic;
            fCS = figure('Visible', 'off', 'Color', 'w');
            cKey = ipfColorKey(S.grains);
            color_vals = cKey.orientation2color(S.grains.meanOrientation);
            plot(S.grains, color_vals, 'FaceAlpha', 0.3, 'linewidth', params.mapLineWidth, 'micronbar', 'off');
            hold on;
            cS_obj = getCrystalShapeObject(S.ebsd.CS);
            isBig = S.grains.area > 500;
            cSG = S.grains(isBig).meanOrientation * cS_obj * 0.7 * sqrt(S.grains(isBig).area/pi);
            plot(S.grains(isBig).centroid + cSG, 'FaceColor', color_vals(isBig,:), 'FaceAlpha', 0.8);
            hold off;
            exportgraphics(fCS, fileCS, 'Resolution', params.exportRes);
            close(fCS);
            fprintf('  -> Crystal shape map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> Crystal shape map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileKAM, 'file')
            tSub = tic;
            kam = S.ebsd.KAM('threshold', 5*degree);
            fKAM = figure('Visible', 'off', 'Color', 'w');
            plot(S.ebsd, kam ./ degree, 'micronbar', 'off');
            hold on; plot(S.grains.boundary, 'color', 'k', 'linewidth', 1); hold off;
            colormap(viridis); mtexColorbar;
            set(gca, 'CLim', [0 2]);
            exportgraphics(fKAM, fileKAM, 'Resolution', params.exportRes);
            close(fKAM);
            fprintf('  -> KAM map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> KAM map exists: skipped\n');
        end

        if params.overwritePlots || ~exist(fileEllipse, 'file')
            tSub = tic;
            [cEllipse, aEllipse, bEllipse] = fitEllipse(S.grains);
            fEllipse = figure('Visible', 'off', 'Color', 'w');
            plot(S.grains.boundary, 'linewidth', 0.4, 'color', [0.4 0.4 0.4]);
            hold on;
            plotEllipse(cEllipse, aEllipse, bEllipse, 'lineColor', 'r', 'linewidth', 1);
            hold off;
            axis equal tight;
            exportgraphics(fEllipse, fileEllipse, 'Resolution', params.exportRes);
            close(fEllipse);
            fprintf('  -> Ellipse fit map exported (%.2f s)\n', toc(tSub));
        else
            fprintf('  -> Ellipse fit map exists: skipped\n');
        end

        if params.saveIntermediateMAT
            step2Mat = fullfile(checkpointDir, sprintf('%s_%s_STEP2_maps.mat', sampleName, phaseName));
            save(step2Mat, 'sampleName', 'phaseName', 'filePF', 'fileIPF', 'fileAR', 'fileSize', 'fileOri', ...
                'fileGOS', 'fileGROD', 'fileHet', 'fileCS', 'fileKAM', 'fileEllipse');
            fprintf('  -> STEP 2 checkpoint saved\n');
        end
    else
        fprintf('  -> Skipped by flag\n');
    end

    fprintf('✔ STEP 2 complete (%.2f s)\n', toc(tStep));

    %% -------------------------------------------------
    %% STEP 3: GRAIN CSV EXPORT
    %% -------------------------------------------------
    tStep = tic;
    fprintf('[STEP 3] Grain CSV export\n');

    if run_exportCSVs
        grainTable = buildGrainTable(S);
        writetable(grainTable, grainCSVFile);

        if isfield(S, 'T')
            summaryRows = [summaryRows; S.T]; %#ok<AGROW>
            writetable(summaryRows, fullfile(csvDir, 'AllSamples_TextureStats.csv'));
        end

        if params.saveIntermediateMAT
            step3Mat = fullfile(checkpointDir, sprintf('%s_%s_STEP3_grainCSV.mat', sampleName, phaseName));
            save(step3Mat, 'grainTable');
        end
        fprintf('  -> Grain CSV exported / overwritten\n');
    else
        fprintf('  -> Skipped by flag\n');
    end

    fprintf('✔ STEP 3 complete (%.2f s)\n', toc(tStep));

    %% -------------------------------------------------
    %% STEP 4.1: CONTACT / PAIR EXPORT (ANORTHITE ONLY)
    %% -------------------------------------------------
    tStep = tic;
    fprintf('[STEP 4.1] Contact / pair export\n');

    if run_contactAnalysis && isAnorthite
        C = struct();
        C.contactTable = buildPairContactTable(S, sampleName, phaseName);
        save(contactMatFile, '-struct', 'C');

        if run_exportCSVs
            writetable(C.contactTable, contactCSVFile);
            fprintf('  -> Contact CSV exported / overwritten\n');
        end
    else
        fprintf('  -> Skipped by flag or phase\n');
    end

    fprintf('✔ STEP 4.1 complete (%.2f s)\n', toc(tStep));

    %% -------------------------------------------------
    %% STEP 4.2: PIXEL EXPORT
    %% -------------------------------------------------
    tStep = tic;
    fprintf('[STEP 4.2] Pixel export\n');

    if run_exportCSVs && run_pixelExport && isAnorthite
        pixelTable = buildPixelTable(S);
        writetable(pixelTable, pixelCSVFile);
        fprintf('  -> Pixel CSV exported / overwritten\n');
    else
        fprintf('  -> Skipped by flag or phase\n');
    end

    fprintf('✔ STEP 4.2 complete (%.2f s)\n', toc(tStep));

    %% -------------------------------------------------
    fprintf('>> Total Processing: %s | %s (%.2f s)\n', sampleName, phaseName, toc(tSample));
    diary off;
end

fprintf('\n>> GRaFT complete.\n');

%% ================= HELPER FUNCTIONS =================

function grainTable = buildGrainTable(S)
% BUILDGRAINTABLE Constructs a grain-resolved export table for one sample-phase.
%
% Inputs:
%   S - Structure containing GRaFT outputs including:
%       grains, eqDiam, longAxisAngle, cs, sampleName, and phaseName.
%
% Outputs:
%   grainTable - Table containing grain-scale geometric, crystallographic,
%                and ellipse-fit parameters for CSV export.
%
% Notes:
%   Exported properties include:
%     - Grain area and equivalent circular diameter (ECD)
%     - Aspect ratio and long-axis orientation
%     - Grain Orientation Spread (GOS)
%     - Projected crystallographic axis orientations
%     - Ellipse-fit centroid, axis lengths, and axis orientations.
    nGrains = S.grains.length;
    longAxisAngle = S.longAxisAngle(:);
    ori = S.grains.meanOrientation;
    cs  = S.cs;

    a_crystal = vector3d(Miller(1,0,0,cs,'uvw'));
    b_crystal = vector3d(Miller(0,1,0,cs,'uvw'));
    c_crystal = vector3d(Miller(0,0,1,cs,'uvw'));

    a_sample = ori * a_crystal;
    b_sample = ori * b_crystal;
    c_sample = ori * c_crystal;

    aAxisAngle = atan2(a_sample.y, a_sample.x) / degree;
    bAxisAngle = atan2(b_sample.y, b_sample.x) / degree;
    cAxisAngle = atan2(c_sample.y, c_sample.x) / degree;

    aAxisAngle(aAxisAngle < 0) = aAxisAngle(aAxisAngle < 0) + 180;
    bAxisAngle(bAxisAngle < 0) = bAxisAngle(bAxisAngle < 0) + 180;
    cAxisAngle(cAxisAngle < 0) = cAxisAngle(cAxisAngle < 0) + 180;

    [cEllipse, aEllipse, bEllipse] = fitEllipse(S.grains);

    majorAxisLength = 2 * aEllipse.norm;
    minorAxisLength = 2 * bEllipse.norm;

    majorAxisAngle = atan2(aEllipse.y, aEllipse.x) / degree;
    minorAxisAngle = atan2(bEllipse.y, bEllipse.x) / degree;

    majorAxisAngle(majorAxisAngle < 0) = majorAxisAngle(majorAxisAngle < 0) + 180;
    minorAxisAngle(minorAxisAngle < 0) = minorAxisAngle(minorAxisAngle < 0) + 180;

    grainTable = table();
    grainTable.Sample = repmat(string(S.sampleName), nGrains, 1);
    grainTable.Phase = repmat(string(S.phaseName), nGrains, 1);
    grainTable.GrainID = S.grains.id(:);
    grainTable.Area = S.grains.area(:);
    grainTable.ECD = S.eqDiam(:);
    grainTable.AspectRatio = S.grains.aspectRatio(:);
    grainTable.LongAxisAngle = longAxisAngle;
    grainTable.GOS = S.grains.GOS(:) / degree;
    grainTable.aAxisAngle = aAxisAngle(:);
    grainTable.bAxisAngle = bAxisAngle(:);
    grainTable.cAxisAngle = cAxisAngle(:);
    grainTable.EllipseCentroidX = cEllipse.x(:);
    grainTable.EllipseCentroidY = cEllipse.y(:);
    grainTable.MajorAxisLength = majorAxisLength(:);
    grainTable.MinorAxisLength = minorAxisLength(:);
    grainTable.MajorAxisAngle = majorAxisAngle(:);
    grainTable.MinorAxisAngle = minorAxisAngle(:);
end

function contactTable = buildPairContactTable(S, sampleName, phaseName)
% BUILDPAIRCONTACTTABLE Constructs a grain-contact relationship table.
%
% Inputs:
%   S          - Structure containing GRaFT outputs including grains,
%                boundaries, orientations, and long-axis data.
%   sampleName - Current sample name.
%   phaseName  - Current phase name.
%
% Outputs:
%   contactTable - Table containing pairwise grain-contact relationships
%                  and crystallographic comparisons.
%
% Notes:
%   The exported table contains:
%     - Grain-ID contact pairs
%     - Misorientation angles
%     - Long-axis angular differences
%     - Projected crystallographic axis angular differences.
%
%   Only valid grain-boundary contacts between indexed grains are retained.
    gB = S.grains.boundary;
    pairIdsAll = sort(gB.grainId, 2);
    valid = all(pairIdsAll > 0, 2);

    if ~any(valid)
        contactTable = table();
        return;
    end

    pairIds = pairIdsAll(valid, :);
    pairUnique = unique(pairIds, 'rows');

    grainIds = S.grains.id(:);
    nGrains = S.grains.length;
    maxId = max(grainIds);
    id2pos = zeros(maxId, 1);
    id2pos(grainIds) = 1:nGrains;

    pos1 = id2pos(pairUnique(:,1));
    pos2 = id2pos(pairUnique(:,2));

    longAxisAngle = S.longAxisAngle(:);

    ori = S.grains.meanOrientation;
    cs = S.cs;

    a_crystal = vector3d(Miller(1,0,0,cs,'uvw'));
    b_crystal = vector3d(Miller(0,1,0,cs,'uvw'));
    c_crystal = vector3d(Miller(0,0,1,cs,'uvw'));

    a_sample = ori * a_crystal;
    b_sample = ori * b_crystal;
    c_sample = ori * c_crystal;

    aAxisAngle = atan2(a_sample.y, a_sample.x) / degree;
    bAxisAngle = atan2(b_sample.y, b_sample.x) / degree;
    cAxisAngle = atan2(c_sample.y, c_sample.x) / degree;

    aAxisAngle(aAxisAngle < 0) = aAxisAngle(aAxisAngle < 0) + 180;
    bAxisAngle(bAxisAngle < 0) = bAxisAngle(bAxisAngle < 0) + 180;
    cAxisAngle(cAxisAngle < 0) = cAxisAngle(cAxisAngle < 0) + 180;

    az1 = longAxisAngle(pos1);
    az2 = longAxisAngle(pos2);
    dAz = abs(az1 - az2);
    dAz = min(dAz, 180 - dAz);

    da = abs(aAxisAngle(pos1) - aAxisAngle(pos2));
    db = abs(bAxisAngle(pos1) - bAxisAngle(pos2));
    dc = abs(cAxisAngle(pos1) - cAxisAngle(pos2));

    da = min(da, 180 - da);
    db = min(db, 180 - db);
    dc = min(dc, 180 - dc);

    misPair = angle(S.grains(pos1).meanOrientation, S.grains(pos2).meanOrientation) ./ degree;

    contactTable = table();
    contactTable.Sample = repmat(string(sampleName), size(pairUnique,1), 1);
    contactTable.Phase = repmat(string(phaseName), size(pairUnique,1), 1);
    contactTable.GrainID_1 = pairUnique(:,1);
    contactTable.GrainID_2 = pairUnique(:,2);
    contactTable.MisorientationDeg = misPair;
    contactTable.DeltaLongAxisDeg = dAz;
    contactTable.DeltaAAxisDeg = da;
    contactTable.DeltaBAxisDeg = db;
    contactTable.DeltaCAxisDeg = dc;
end

function pixelTable = buildPixelTable(S)
% BUILDPIXELTABLE Constructs a pixel-resolved deformation export table.
%
% Inputs:
%   S - Structure containing GRaFT outputs including EBSD and grain objects.
%
% Outputs:
%   pixelTable - Table containing pixel-scale deformation metrics for CSV export.
%
% Notes:
%   Exported properties include:
%     - Grain IDs
%     - Grain Reference Orientation Deviation (GROD)
%     - Kernel Average Misorientation (KAM)
%     - GROD/GOS ratio.
%
%   GROD/GOS ratios are stabilised using a minimum denominator threshold
%   to avoid division by near-zero GOS values.
    grod_obj = S.ebsd.calcGROD(S.grains);
    grod_deg = grod_obj.angle / degree;
    kam = S.ebsd.KAM('threshold', 5*degree) / degree;

    isIndexed = S.ebsd.grainId > 0;
    gos_pixel = nan(size(grod_deg));
    gos_pixel(isIndexed) = S.grains(S.ebsd.grainId(isIndexed)).GOS / degree;

    grod_gos_ratio = nan(size(grod_deg));
    grod_gos_ratio(isIndexed) = grod_deg(isIndexed) ./ max(gos_pixel(isIndexed), 0.05);

    pixelID = (1:numel(S.ebsd.grainId))';

    pixelTable = table();
    pixelTable.Sample = repmat(string(S.sampleName), numel(pixelID), 1);
    pixelTable.Phase = repmat(string(S.phaseName), numel(pixelID), 1);
    pixelTable.PixelID = pixelID;
    pixelTable.GrainID = S.ebsd.grainId(:);
    pixelTable.GROD = grod_deg(:);
    pixelTable.KAM = kam(:);
    pixelTable.GROD_GOS_Ratio = grod_gos_ratio(:);
end

function plotCPOAnalysis(odf, h, cs, grains, useODF, outDir, prefix, params)
% PLOTCPOANALYSIS Generates pole figure and inverse pole figure maps.
%
% Inputs:
%   odf      - Orientation Distribution Function object.
%   h        - Cell array of crystallographic directions.
%   cs       - Crystal symmetry object.
%   grains   - MTEX grain object.
%   useODF   - Logical flag controlling ODF-based or discrete plotting.
%   outDir   - Output directory for exported figures.
%   prefix   - Filename prefix for exported maps.
%   params   - Global plotting parameter structure.
%
% Outputs:
%   Exports:
%     - Pole figure (PF) map
%     - Inverse pole figure (IPF) map.
%
% Notes:
%   If insufficient grains exist for robust ODF calculation, discrete
%   orientation plots are generated instead of contour-density maps.
    fPF = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 400]);

    if useODF
        plotPDF(odf, h, 'antipodal', 'contourf');
        colormap(min(parula + 0.3, 1));
        mtexColorbar;
    else
        plotPDF(grains.meanOrientation, h, 'antipodal', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k');
    end

    exportgraphics(fPF, fullfile(outDir, [prefix '_01_PF.png']), 'Resolution', params.exportRes);
    close(fPF);

    fIPF = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 500]);

    if useODF
        plotIPDF(odf, [vector3d.X, vector3d.Y, vector3d.Z], 'antipodal', 'contourf', 'noLabel');
        colormap(min(parula + 0.3, 1));
        mtexColorbar;
    else
        plotIPDF(grains.meanOrientation, [vector3d.X, vector3d.Y, vector3d.Z], 'antipodal', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k', ...
            'noLabel');
    end

    hLabels = [Miller(1,0,0,cs), Miller(0,1,0,cs), Miller(0,0,1,cs)];
    axisLabels = {'X', 'Y', 'Z'};

    for i = 1:3
        nextAxis(i);
        hold on;
        title(axisLabels{i}, 'FontWeight', 'bold', 'FontSize', 12);
        plot(hLabels, 'marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', 'markersize', 8);
        annotate(hLabels, 'label', {'100', '010', '001'}, 'color', 'k', 'FontSize', 10, 'FontWeight', 'bold', 'labeloffset', 0.12);
    end

    exportgraphics(fIPF, fullfile(outDir, [prefix '_02_IPF.png']), 'Resolution', params.exportRes);
    close(fIPF);
end

function [pfJ, pfMax] = calcPFStats(odf, hkl, res)
% CALCPFSTATS Calculates pole-figure texture statistics.
%
% Inputs:
%   odf - Orientation Distribution Function object.
%   hkl - Cell array of Miller indices defining pole figures.
%   res - Spherical grid resolution.
%
% Outputs:
%   pfJ   - Pole-figure J-indices for each crystallographic direction.
%   pfMax - Maximum pole density values for each pole figure.
%
% Notes:
%   Pole-figure statistics are calculated on a regular spherical grid
%   using quadrature-weighted integration.
    r = regularS2Grid('resolution', res);
    w = calcQuadratureWeights(r);
    pfJ = zeros(numel(hkl),1);
    pfMax = zeros(numel(hkl),1);

    for i = 1:numel(hkl)
        pf = calcPDF(odf, hkl{i}, r);
        pfJ(i) = sum(pf.^2 .* w);
        pfMax(i) = max(pf);
    end
end

function [BA, BC] = calcFabricIndices(odf, hkl, res)
% CALCFABRICINDICES Calculates BA and BC fabric-shape indices.
%
% Inputs:
%   odf - Orientation Distribution Function object.
%   hkl - Cell array of Miller indices defining analysed pole figures.
%   res - Spherical grid resolution.
%
% Outputs:
%   BA - Fabric symmetry index comparing [100] and [010] fabrics.
%   BC - Fabric symmetry index comparing [001] and [010] fabrics.
%
% Notes:
%   Fabric indices are derived from eigenvalue analysis of orientation
%   tensors calculated from pole-density distributions.
    lambda = zeros(3, numel(hkl));
    r = regularS2Grid('resolution', res);
    w = calcQuadratureWeights(r);

    for i = 1:numel(hkl)
        pf = calcPDF(odf, hkl{i}, r);
        T = zeros(3,3);

        for j = 1:length(r)
            v = double(r(j))';
            T = T + pf(j) * (v * v') * w(j);
        end

        lambda(:,i) = sort(eig(T / sum(pf .* w)), 'descend');
    end

    P = lambda(1,:) - lambda(3,:);
    G = 2 * (lambda(2,:) - lambda(3,:));
    BA = 0.5 * (2 - P(2)/(P(2)+G(2)) - P(1)/(P(1)+G(1)));
    BC = 0.5 * (2 - P(2)/(P(2)+G(2)) - P(3)/(P(3)+G(3)));
end

function cS_obj = getCrystalShapeObject(CS)
% GETCRYSTALSHAPEOBJECT Returns a crystal-shape object for plotting.
%
% Inputs:
%   CS - MTEX crystal symmetry object.
%
% Outputs:
%   cS_obj - MTEX crystalShape object corresponding to the mineral phase.
%
% Notes:
%   Currently supported:
%     - Anorthite
%     - Forsterite
%     - Diopside
%
%   Unsupported phases default to olivine crystal shapes.
    mineralName = char(CS.mineral);

    switch lower(mineralName)
        case {'anorthite'}
            cS_obj = crystalShape.plagioclase(CS);
        case {'forsterite'}
            cS_obj = crystalShape.forsterite;
        case {'diopside'}
            cS_obj = crystalShape.diopside;
        otherwise
            cS_obj = crystalShape.olivine;
    end
end

function cmap = plasma(n)
% PLASMA Returns an interpolated plasma-style colormap.
%
% Inputs:
%   n - Number of colours requested. Defaults to 256 if omitted.
%
% Outputs:
%   cmap - n-by-3 RGB colormap array.
%
% Notes:
%   Colours are linearly interpolated from a reduced plasma anchor palette.
    if nargin < 1, n = 256; end

    base = [ ...
        0.0504 0.0298 0.5280
        0.2546 0.0139 0.6154
        0.4176 0.0006 0.6584
        0.5627 0.0515 0.6415
        0.6928 0.1651 0.5645
        0.7982 0.2802 0.4695
        0.8814 0.3925 0.3832
        0.9440 0.5532 0.2871
        0.9796 0.7048 0.2129
        0.9892 0.8463 0.1400];

    x = linspace(0,1,size(base,1));
    xi = linspace(0,1,n);
    cmap = interp1(x, base, xi, 'linear');
end

function cmap = viridis(n)
% VIRIDIS Returns an interpolated viridis-style colormap.
%
% Inputs:
%   n - Number of colours requested. Defaults to 256 if omitted.
%
% Outputs:
%   cmap - n-by-3 RGB colormap array.
%
% Notes:
%   Colours are linearly interpolated from a reduced viridis anchor palette.
    if nargin < 1, n = 256; end

    base = [ ...
        0.2670 0.0049 0.3294
        0.2823 0.1409 0.4575
        0.2539 0.2653 0.5300
        0.2068 0.3718 0.5531
        0.1636 0.4711 0.5581
        0.1276 0.5669 0.5506
        0.1347 0.6586 0.5176
        0.2669 0.7488 0.4406
        0.4775 0.8214 0.3182
        0.7414 0.8734 0.1496];

    x = linspace(0,1,size(base,1));
    xi = linspace(0,1,n);
    cmap = interp1(x, base, xi, 'linear');
end
