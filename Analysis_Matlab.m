
clear all;

% Prompt user to select a .csv file
[filename, filepath] = uigetfile('*.csv', 'Select a CSV file');

% Check if user clicked cancel
if isequal(filename, 0) || isequal(filepath, 0)
    disp('User canceled the operation.');
    return;
end

% Construct the full file path
fullFilePath = fullfile(filepath, filename);

% Read the file as a table
T = readtable(fullFilePath, 'Delimiter', ',');


% Extract the 4th column as a numerical vector
column4 = table2array(T(:, 4));


% Check if the number of elements is divisible by 3
isDivisibleBy3 = mod(numel(column4), 3) == 0;

% make new variable with data from separate channels
if isDivisibleBy3
    numRows = numel(column4) / 3;
    % Reshape the column vector into a matrix with three columns and arranged row-wise
    threeColumnMatrix = reshape(column4, 3, numRows)';
else
    disp('The number of elements is not divisible by 3.');
end

%Store red channel data in red channel variable:
redchannel=threeColumnMatrix(:, 1);

% Check if the number of elements is divisible by 45, which is the number
% of time frames in the file
isDivisibleBy45 = mod(numel(redchannel), 45) == 0;

% make new variable with data from separate cells
if isDivisibleBy45
    % Reshape the variable to have 45 rows and as many columns as needed
    redchannel_table = reshape(redchannel, 45, []);
else
    disp('The number of elements is not divisible by 3.');
end

% Now I want to make a figure to show all the traces simultaneously. Loop through each column and plot

f=figure;
hold on;
for i = 1:size(redchannel_table, 2)
    % Extract the data from the current column
    thiscell = redchannel_table(:, i);
    
    % Create the plot with a semi-transparent line
    plot(thiscell, 'LineWidth', 2, 'Color', [0 0.4470 0.7410 0.5]); % RGBA color with alpha for transparency
end
hold off;

%Prepare for saving the figure with all plots;
figname=filename;
figname=matlab.lang.makeValidName(figname);
figname=strcat(figname,'_all_traces.png');
base_path = 'L:\LemkeLab\SaraMingu\Progeria Transport Assay Project\Analysis\Results from Analysis_Matlab';
figpath = fullfile(base_path, figname);
saveas(f, figpath);
close(f);
% 
% x = (2:2:90)';
% 
% %Make a new matrix to store fit results
% redchannel_table_fit_results=zeros(9,size(redchannel_table, 2));
% %Make a new matrix to store best fit lines; here i just copy the original
% %one and will over-write later
% redchannel_table_bestfits=redchannel_table;
% 
% 
% for i=1:size(redchannel_table, 2)
%     y = redchannel_table(:, i);
%    
%     % Ensure there are no NaN values in y
%     if any(isnan(y))
%         disp(['Skipping column ', num2str(i), ' due to NaN values.']);
%         continue;
%     end
%    
%     % Fit data
%     f = fittype('a + imax*(1-exp(-tau*x))', 'independent', 'x', 'coefficients', {'a', 'imax', 'tau'});
%     initial_guess = [y(1), max(y) - y(1), 0.1];
%    
%     try
%         [cfun, rsquare] = fit(x, y, f, 'StartPoint', initial_guess, 'Lower', [0, 0, 0], 'Upper', [Inf, Inf, Inf]);
%         
%         % Optionally, store or process the fit results
%     catch ME
%         disp(['Skipping column ', num2str(i), ' due to fitting error: ', ME.message]);
%         continue;
%     end
%     
%     %store results in a separate variable redchannel_table_fit_results
%     redchannel_table_fit_results(1,i)=cfun.a;
%     redchannel_table_fit_results(2,i)=cfun.imax;
%     redchannel_table_fit_results(3,i)=cfun.tau;
%     redchannel_table_fit_results(4,i)=rsquare.sse;
%     redchannel_table_fit_results(5,i)=rsquare.rsquare;
%     redchannel_table_fit_results(6,i)=rsquare.dfe;
%     redchannel_table_fit_results(7,i)=rsquare.adjrsquare;
%     redchannel_table_fit_results(8,i)=rsquare.rmse;
%     min_y=min(y);
%     max_y=max(y);
%     import_fraction = (max_y-min_y)/min_y;
%     redchannel_table_fit_results(9,i)=import_fraction;
% 
%     f2=figure;
%     yi=cfun(x);
%     
%     %Add the best fit lines to new table
%     redchannel_table_bestfits(:,i)=yi;
%     
%     plot(x,y,'r*',x,yi,'b-');
%     %Save figure
%     figname=filename;
%     figname=matlab.lang.makeValidName(figname);
%     %add information of i
%     thiscell_nr=num2str(i); 
%     figname=strcat(figname, '_', thiscell_nr, '.png');
%     figpath = fullfile(base_path, figname);
%     saveas(f2, figpath);
%     close(f2);
% end 

% export files

base_path = 'L:\LemkeLab\SaraMingu\Progeria Transport Assay Project\Analysis\Results from Analysis_Matlab';

filename_traces=matlab.lang.makeValidName(filename);
filename_traces=strcat(filename_traces,'.xlsx');
redchannel_table=array2table(redchannel_table);
writetable(redchannel_table, fullfile(base_path,filename_traces));

% filename_bestfits=matlab.lang.makeValidName(filename);
% filename_bestfits=strcat(filename_bestfits,'_bestfits.xlsx');
% redchannel_table_bestfits=array2table(redchannel_table_bestfits);
% writetable(redchannel_table_bestfits, fullfile(base_path,filename_bestfits));
% 
% filename_2=matlab.lang.makeValidName(filename);
% filename_2=strcat(filename_2,'_fittingresults.xlsx');
% redchannel_table_fit_results=array2table(redchannel_table_fit_results);
% %add the column with labels
% writetable(redchannel_table_fit_results, fullfile(base_path,filename_2));