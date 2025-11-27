function fit_capsid_timecourses_from_excel()
% fit_capsid_timecourses_from_excel.m
%
% - Ask user to select an Excel file with timecourses:
%       * Each SHEET = one experiment (keep sheet names)
%       * Each COLUMN = one cell
%       * Each ROW    = one time point (1..45), corresponding to 2..90 min
% - Capsid type is inferred from the filename (19nls, 23nls, 29nls, 38nls),
%   which sets the normalization factor.
% - For each sheet & column (cell):
%       1) Subtract fixed background (BG_VALUE).
%       2) Normalize by capsid-specific factor.
%       3) Detect delay region (flat or dipping).
%       4) Fit rising part to:
%              I(t) = A + IMAX * (1 - exp(-k * (t - delay)))
%          where delay is the onset of monotonic increase.
%       5) Compute A, IMAX, k, J = IMAX*k, delay (min), R^2.
%       6) Save per-cell PNG with trace + fit + R^2.
% - Create a new Excel file with same sheet names:
%       Rows = cells, Columns = [CellID, delay_min, A, IMAX, k, J, R2].
% - Make dot plots for k, IMAX, J, delay, grouped by sheet (experiment).

    clc;
    fprintf('--- fit_capsid_timecourses_from_excel ---\n');

    %% Choose Excel file
    [fileName, filePath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
                                     'Select Excel file with timecourses');
    if isequal(fileName,0)
        error('No Excel file selected.');
    end
    excelFile = fullfile(filePath, fileName);
    [~, baseName, ~] = fileparts(excelFile);
    fprintf('Selected file: %s\n', excelFile);

    %% Determine capsid type & normalization factor from filename
    normFactor = infer_normalization_factor_from_filename(fileName);

    fprintf('Normalization factor inferred from filename = %.3f\n', normFactor);

    %% Fixed background value (you can adjust this)
    BG_VALUE = 108.5;
    fprintf('Using fixed background value: %.2f\n', BG_VALUE);

    %% Get sheet names
    [~, sheetNames] = xlsfinfo(excelFile);
    if isempty(sheetNames)
        error('No sheets found in Excel file.');
    end

    % Output folder for PNGs and for result Excel
    outDir = fullfile(filePath, [baseName '_fit_output']);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    resultExcel = fullfile(outDir, sprintf('%s_fitResults.xlsx', baseName));

    % For grouped dot plots across sheets
    allSheetNames = {};
    allK   = {};
    allI   = {};
    allJ   = {};
    allDel = {};

    %% Process each sheet (experiment)
    for s = 1:numel(sheetNames)
        shName = sheetNames{s};
        fprintf('\n--- Sheet %d/%d: %s ---\n', s, numel(sheetNames), shName);

        % Read numeric data from sheet (ignores header rows/cols)
        [numData, ~, ~] = xlsread(excelFile, shName);
        if isempty(numData)
            fprintf('  Sheet %s: no numeric data, skipping.\n', shName);
            continue;
        end

        [nTimePoints, nCells] = size(numData);
        fprintf('  Found %d time points x %d cells.\n', nTimePoints, nCells);

        % Build time vector: 2,4,...,2*nTimePoints (minutes)
        tVec = (1:nTimePoints)' * 2;

        % Preallocate fit parameter arrays
        delay_min = nan(nCells, 1);
        A_vals    = nan(nCells, 1);
        Imax_vals = nan(nCells, 1);
        k_vals    = nan(nCells, 1);
        J_vals    = nan(nCells, 1);
        R2_vals   = nan(nCells, 1);

        % Sheet-specific output folder for PNGs
        sheetSafeName = sanitize_for_filename(shName);
        sheetDir = fullfile(outDir, sheetSafeName);
        if ~exist(sheetDir, 'dir'), mkdir(sheetDir); end

        %% Per-cell fit
        for c = 1:nCells
            y_raw = numData(:, c);

            % Skip empty or all-NaN columns
            if all(isnan(y_raw)) || all(y_raw == 0)
                fprintf('  Cell %d: empty or all zero, skipping.\n', c);
                continue;
            end

            % Background subtraction
            y_bgsub = y_raw - BG_VALUE;
            y_bgsub(y_bgsub < 0) = 0;

            % Normalize by capsid factor
            y_norm = y_bgsub / normFactor;

            % Detect delay & fit
            fitRes = fit_single_cell_trace(tVec, y_norm);

            delay_min(c) = fitRes.delay;
            A_vals(c)    = fitRes.A;
            Imax_vals(c) = fitRes.Imax;
            k_vals(c)    = fitRes.k;
            J_vals(c)    = fitRes.J;
            R2_vals(c)   = fitRes.R2;

            % Plot trace and fit with R^2 and save PNG
            pngName = sprintf('%s_%s_col%02d.png', baseName, sheetSafeName, c);
            pngPath = fullfile(sheetDir, pngName);
            plot_and_save_cell_fit(tVec, y_norm, fitRes, shName, c, pngPath);
        end

        %% Save fit results for this sheet to Excel
        headers = {'CellID','Delay_min','A','IMAX','k','J','R2'};
        data = [(1:nCells)', delay_min, A_vals, Imax_vals, k_vals, J_vals, R2_vals];
        outCell = [headers; num2cell(data)];
        writecell(outCell, resultExcel, 'Sheet', shName);

        % Store for grouped dot plots
        allSheetNames{end+1} = shName; %#ok<AGROW>
        allK{end+1}   = k_vals;        %#ok<AGROW>
        allI{end+1}   = Imax_vals;     %#ok<AGROW>
        allJ{end+1}   = J_vals;        %#ok<AGROW>
        allDel{end+1} = delay_min;     %#ok<AGROW>
    end

    %% Grouped dot plots across sheets (experiments)
    fprintf('\nCreating grouped dot plots across sheets...\n');

    if ~isempty(allSheetNames)
        make_grouped_dotplot(allSheetNames, allK,   'k (1/min)',      fullfile(outDir, [baseName '_k_dotplot.png']));
        make_grouped_dotplot(allSheetNames, allI,   'IMAX (norm)',    fullfile(outDir, [baseName '_IMAX_dotplot.png']));
        make_grouped_dotplot(allSheetNames, allJ,   'J = IMAX*k',     fullfile(outDir, [baseName '_J_dotplot.png']));
        make_grouped_dotplot(allSheetNames, allDel, 'Delay (min)',    fullfile(outDir, [baseName '_Delay_dotplot.png']));
    end

    fprintf('\nAll done.\nFit results Excel: %s\nPlots in: %s\n', resultExcel, outDir);
end

%% ========================================================================
%% Helper: infer normalization factor from filename
function nf = infer_normalization_factor_from_filename(fname)
    lowerName = lower(fname);

    if contains(lowerName, '19nls')
        nf = 29;
    elseif contains(lowerName, '23nls')
        nf = 27.7;
    elseif contains(lowerName, '29nls')
        nf = 26.8;
    elseif contains(lowerName, '38nls')
        nf = 38.2;
    else
        warning('Could not infer NLS from filename "%s". Using normalization factor = 1.', fname);
        nf = 1;
    end
end

%% ========================================================================
%% Helper: sanitize string to be safe as filename
function s = sanitize_for_filename(str)
    s = str;
    % Replace problematic characters with underscore
    s = regexprep(s, '[/\\:*?"<>|]', '_');
    s = strrep(s, ' ', '_');
end

%% ========================================================================
%% Helper: fit a single cell trace with delay detection
function fitRes = fit_single_cell_trace(tVec, y_norm)
    % Returns struct with fields:
    %   A, Imax, k, J, delay, R2

    fitRes = struct('A', NaN, 'Imax', NaN, 'k', NaN, ...
                    'J', NaN, 'delay', NaN, 'R2', NaN);

    % Remove NaNs
    valid = ~isnan(y_norm);
    t = tVec(valid);
    y = y_norm(valid);

    n = numel(t);
    if n < 5
        return;
    end

    % If basically flat and tiny amplitude => no meaningful fit
    if (max(y) - min(y)) < 0.01 * max(1, max(y))
        return;
    end

    % Smooth for delay detection
    ys = smooth_trace(y, 3);
    diffs = diff(ys);

    % Delay detection: look for first index where
    %  - local slope is non-negative (within tolerance)
    %  - from there to the end we gain enough amplitude
    delayIdx = 1;
    yMax = max(ys);
    smallTol = 0.01 * max(yMax, 1);   % tolerance for small negative noise
    minAmp = 0.05 * yMax;             % need at least 5% of max amplitude after delay

    for k = 1:(n-3)
        windowDiffs = diffs(k:k+2);
        if all(windowDiffs >= -smallTol) && (ys(end) - ys(k)) > minAmp
            delayIdx = k;
            break;
        end
    end

    delayTime = t(delayIdx);  % in minutes

    tFit = t(delayIdx:end) - t(delayIdx);  % start at 0
    yFit = y(delayIdx:end);

    if numel(tFit) < 4
        fitRes.delay = delayTime;
        return;
    end

    % Initial parameter guesses
    A0 = yFit(1);                          % baseline offset at onset
    Imax0 = max(yFit) - A0;
    if Imax0 <= 0
        Imax0 = max(yFit);
    end
    if Imax0 <= 0
        % still flat; give up
        fitRes.delay = delayTime;
        return;
    end
    k0 = 1 / max(tFit);                    % rough guess

    p0 = [A0, Imax0, k0];

    % Model: I(t) = A + Imax * (1 - exp(-k * t))
    modelFun = @(p, tt) p(1) + p(2) .* (1 - exp(-p(3) .* tt));
    errFun   = @(p) sum((modelFun(p, tFit) - yFit).^2);

    opts = optimset('Display','off');
    [pFit, ~, exitflag] = fminsearch(errFun, p0, opts);
    if exitflag <= 0
        pFit = p0;   % if optimization failed, fall back to initial
    end

    A     = pFit(1);
    Imax  = max(pFit(2), 0);
    k     = max(pFit(3), 0);
    yHat  = modelFun([A, Imax, k], tFit);

    % R^2
    SSres = sum((yFit - yHat).^2);
    SStot = sum((yFit - mean(yFit)).^2);
    if SStot > 0
        R2 = 1 - SSres / SStot;
    else
        R2 = NaN;
    end

    J = Imax * k;

    fitRes.A     = A;
    fitRes.Imax  = Imax;
    fitRes.k     = k;
    fitRes.J     = J;
    fitRes.delay = delayTime;
    fitRes.R2    = R2;
end

%% ========================================================================
%% Helper: simple moving average smoother
function ys = smooth_trace(y, win)
    if nargin < 2
        win = 3;
    end
    n = numel(y);
    ys = zeros(size(y));
    half = floor(win / 2);
    for i = 1:n
        i1 = max(1, i-half);
        i2 = min(n, i+half);
        ys(i) = mean(y(i1:i2));
    end
end

%% ========================================================================
%% Helper: plot and save per-cell fit
function plot_and_save_cell_fit(tVec, y_norm, fitRes, sheetName, cellID, pngPath)
    % full fitted curve across all time points
    yFitFull = nan(size(y_norm));
    if ~isnan(fitRes.A) && ~isnan(fitRes.k)
        % Reconstruct fit over full time vector
        t0 = fitRes.delay;
        t_shift = max(tVec - t0, 0);  % 0 before delay
        yFitFull = fitRes.A + fitRes.Imax .* (1 - exp(-fitRes.k .* t_shift));
    end

    fig = figure('Visible','off','Color','w');
    hold on;
    plot(tVec, y_norm, 'ko-', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName','Data');
    if ~all(isnan(yFitFull))
        plot(tVec, yFitFull, 'r-', 'LineWidth', 1.5, 'DisplayName','Fit');
    end
    xlabel('Time (min)');
    ylabel('Normalized intensity (bg-subtracted / normFactor)');
    title(sprintf('%s | %s | Cell %d | R^2 = %.3f', ...
          sheetName, pngPathShort(pngPath), cellID, fitRes.R2), ...
          'Interpreter','none');
    legend('Location','best');
    grid on;
    hold off;

    % Save PNG
    try
        exportgraphics(fig, pngPath, 'Resolution', 200);
    catch
        % fallback
        saveas(fig, pngPath);
    end
    close(fig);
end

%% small helper for cleaner title path
function s = pngPathShort(p)
    [~, name, ext] = fileparts(p);
    s = [name ext];
end

%% ========================================================================
%% Helper: grouped dot plots by sheet
function make_grouped_dotplot(sheetNames, paramCells, yLabel, outPng)
    nSheets = numel(sheetNames);
    if nSheets == 0
        return;
    end

    fig = figure('Visible','off','Color','w');
    hold on;

    for s = 1:nSheets
        vals = paramCells{s};
        vals = vals(~isnan(vals));
        if isempty(vals)
            continue;
        end
        x = s + 0.15 * randn(numel(vals), 1);  % jitter around sheet index
        scatter(x, vals, 'filled');
    end

    xlim([0.5, nSheets + 0.5]);
    set(gca, 'XTick', 1:nSheets, 'XTickLabel', sheetNames, 'XTickLabelRotation', 45);
    ylabel(yLabel);
    title(sprintf('%s grouped by sheet (experiment)', yLabel), 'Interpreter','none');
    grid on;
    hold off;

    try
        exportgraphics(fig, outPng, 'Resolution', 200);
    catch
        saveas(fig, outPng);
    end
    close(fig);
end
