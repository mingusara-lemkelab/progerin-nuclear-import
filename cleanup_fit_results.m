function cleanup_fit_results()
% cleanup_fit_results.m
%
% Post-process fitResults Excel files produced by the fitting pipeline.
%
% For each selected Excel file:
%   - Each sheet is one experiment, with columns including:
%       CellID, Delay_min, A, IMAX, k, J, R2_full  (names may vary in case)
%   - Steps per sheet:
%       1) Remove rows with R2_full < 0.9.
%       2) On the remaining rows, compute 90% "CI" for IMAX and k as:
%              mean ± 1.645 * std
%          and remove rows where IMAX or k are outside their respective CI.
%   - Save cleaned data into <basename>_cleaned.xlsx (same sheet names).
%   - Compute mean and std of IMAX, k, J for the cleaned data per sheet and
%     write them into <basename>_summary_cleaned.xlsx:
%         SheetName | Imax_mean | Imax_std | k_mean | k_std | J_mean | J_std
%
% You can select one or multiple Excel files in the dialog.

    clc;
    fprintf('--- cleanup_fit_results.m ---\n');

    %% Ask user for Excel files
    [fnames, fpath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
                                'Select fitResults Excel file(s)', ...
                                'MultiSelect','on');
    if isequal(fnames,0)
        fprintf('No files selected. Aborting.\n');
        return;
    end
    if ischar(fnames)
        fnames = {fnames};
    end

    % Loop over selected files
    for kFile = 1:numel(fnames)
        srcName = fnames{kFile};
        srcPath = fullfile(fpath, srcName);
        [baseName, ~, ~] = fileparts(srcName);

        fprintf('\nProcessing file: %s\n', srcPath);

        % New filenames
        cleanedFile  = fullfile(fpath, sprintf('%s_cleaned.xlsx', baseName));
        summaryFile  = fullfile(fpath, sprintf('%s_summary_cleaned.xlsx', baseName));

        % Get sheet names
        [~, sheetNames] = xlsfinfo(srcPath);
        if isempty(sheetNames)
            fprintf('  No sheets found in %s, skipping.\n', srcPath);
            continue;
        end

        % Prepare summary container
        summaryRows = {};
        summaryHeader = {'SheetName','Imax_mean','Imax_std','k_mean','k_std','J_mean','J_std'};

        % Process each sheet
        for s = 1:numel(sheetNames)
            shName = sheetNames{s};
            fprintf('  Sheet %d/%d: %s\n', s, numel(sheetNames), shName);

            T = readtable(srcPath, 'Sheet', shName);

            if isempty(T) || height(T) == 0
                fprintf('    Empty sheet, skipping.\n');
                continue;
            end

            % Identify relevant columns by (case-insensitive) name
            vnames = T.Properties.VariableNames;

            idxImax = find(strcmpi(vnames, 'IMAX') | strcmpi(vnames, 'Imax') | strcmpi(vnames, 'I_max'), 1);
            idxK    = find(strcmpi(vnames, 'k'), 1);
            idxJ    = find(strcmpi(vnames, 'J'), 1);
            idxR2   = find(strcmpi(vnames, 'R2') | strcmpi(vnames, 'R2_full') | strcmpi(vnames, 'Rsq'), 1);

            if isempty(idxImax) || isempty(idxK) || isempty(idxJ) || isempty(idxR2)
                warning('    Could not find IMAX, k, J, or R2 columns in sheet %s. Skipping this sheet.', shName);
                continue;
            end

            Imax = T{:, idxImax};
            kval = T{:, idxK};
            Jval = T{:, idxJ};
            R2   = T{:, idxR2};

            % Initial good-fit mask: R2 >= 0.9
            goodFitMask = (R2 >= 0.9) & ~isnan(Imax) & ~isnan(kval);

            if nnz(goodFitMask) < 3
                fprintf('    Fewer than 3 good-fit rows (R2>=0.9). Keeping only those (no outlier removal).\n');
                finalMask = goodFitMask;
            else
                % Compute 90% CI for Imax and k on good-fit data:
                % mean ± 1.645 * std
                z90 = 1.645;

                I_good = Imax(goodFitMask);
                k_good = kval(goodFitMask);

                muI = mean(I_good, 'omitnan');
                sdI = std(I_good,  'omitnan');
                muk = mean(k_good, 'omitnan');
                sdk = std(k_good,  'omitnan');

                if sdI == 0
                    inCI_I = true(size(Imax));
                else
                    inCI_I = (Imax >= (muI - z90*sdI)) & (Imax <= (muI + z90*sdI));
                end

                if sdk == 0
                    inCI_k = true(size(kval));
                else
                    inCI_k = (kval >= (muk - z90*sdk)) & (kval <= (muk + z90*sdk));
                end

                finalMask = goodFitMask & inCI_I & inCI_k;
            end

            Tclean = T(finalMask, :);

            % Write cleaned table to cleanedFile (overwrite/create per sheet)
            if ~isempty(Tclean)
                writetable(Tclean, cleanedFile, 'Sheet', shName, 'WriteMode', 'overwrite');
            else
                % If no rows left, write an empty sheet with headers
                writetable(T(1,:), cleanedFile, 'Sheet', shName, 'WriteMode', 'overwrite', ...
                           'WriteVariableNames', true);
                % then immediately delete data rows by writing empty after header
                writetable(T([],:), cleanedFile, 'Sheet', shName, 'WriteMode', 'overwrite', ...
                           'WriteVariableNames', true);
            end

            % Compute summary stats on cleaned data
            if ~isempty(Tclean)
                I_clean = Tclean{:, idxImax};
                k_clean = Tclean{:, idxK};
                J_clean = Tclean{:, idxJ};

                Imean = mean(I_clean, 'omitnan');
                Isd   = std(I_clean,  'omitnan');
                kmean = mean(k_clean, 'omitnan');
                ksd   = std(k_clean,  'omitnan');
                Jmean = mean(J_clean, 'omitnan');
                Jsd   = std(J_clean,  'omitnan');
            else
                Imean = NaN; Isd = NaN;
                kmean = NaN; ksd = NaN;
                Jmean = NaN; Jsd = NaN;
            end

            summaryRows(end+1, :) = {shName, Imean, Isd, kmean, ksd, Jmean, Jsd}; %#ok<AGROW>
        end

        % Write summary Excel
        if ~isempty(summaryRows)
            summaryCell = [summaryHeader; summaryRows];
            writecell(summaryCell, summaryFile, 'Sheet', 'Summary');
            fprintf('  Saved cleaned data to:   %s\n', cleanedFile);
            fprintf('  Saved summary stats to:  %s\n', summaryFile);
        else
            fprintf('  No summary rows generated for %s.\n', srcPath);
        end
    end

    fprintf('\nDone.\n');
end
