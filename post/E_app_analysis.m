function E_app_analysis()
%%% ====================================================================== %%%
% This function will...
% Take calculated cell modulus data from AFM_POST_JON (as .mat)
% Filter out unwanted data
% Check/eliminate if cell data is incomplete
% Append apparent cell modulus to cell file (.mat)
% NOTE: Working directory should containly only .mat files analyzed by one method (e.g. raw_pointwise=0 from AFM_INPUT)
%
%   Appends:
%       1. E_apparent_med (initial median apparent modulus)
%       2. E_apparent_stats (initial descriptive statistics of E_Matrix)
%       3. E_filt_Matrix
%       4. E_apparent_filt_med (median apparent modulus after filtering) **USE THIS VALUE FOR CELL**
%       5. % of NaNs
%       6. FileName
%
%   OPTIONS:
%       Add: ?? rsq_cutoff, outlier_filters, E_app_stats ??
%
% SAVE AS: Appended variables to same .mat file
% 
%%%% Robert J. Wiener (c) Jan 2023 %%%%
%=========================================================================%

%%% INITIALIZE DIRECTORY FILE LOAD %%%
% for-loop through all files in directory
%mat = dir('*.mat'); for q = 1:length(mat) cont = load(mat(q).name); end
Files=dir('*.mat*');
for k=1:length(Files)
   FileName=Files(k).name


% load file
load(FileName)


% NaN zeros
E_Matrix(E_Matrix==0) = NaN;


%%% E_apparent_stats %%%
% Calculate initial stats on E_Matrix (gives SD, range, etc. for deciding how to move forward)
E_apparent_stats{1,1} = 'median'; E_apparent_stats{2,1} = nanmedian(E_Matrix, 'all');          % median
E_apparent_stats{1,2} = 'mean'; E_apparent_stats{2,2} = nanmean(E_Matrix, 'all');              % mean
E_apparent_stats{1,3} = 'min'; E_apparent_stats{2,3} = nanmin(E_Matrix, [], 'all');            % min
E_apparent_stats{1,4} = 'max'; E_apparent_stats{2,4} = nanmax(E_Matrix, [], 'all');            % max
E_apparent_stats{1,5} = 'range'; E_apparent_stats{2,5} = range(E_Matrix, 'all');               % range
E_apparent_stats{1,6} = 'SD'; E_apparent_stats{2,6} = nanstd(E_Matrix, 0, 'all');              % SD
E_apparent_stats{1,7} = 'outliers'; E_apparent_stats{2,7} = isoutlier(E_Matrix);               % outliers


%%% E_apparent_med %%%
% calculate initial median apparent modulus
E_apparent_med = nanmedian(E_Matrix, 'all');



%%% E_apparent_filt_med %%%

%%% Do this if linear fit (has rsq_Matrix)
% 1. Filter based on rsq values
% 2. Filter outliers ( E_apparent_stats{1,7} )
% 3. NaN "E_apparent_filt_med" if over half is bad
E_filt_Matrix = E_Matrix; % initialize filter matrix
if exist('rsq_Matrix','var')

	% 1. NaN filter E values with under 95% R2
	E_filt_Matrix(rsq_Matrix<0.95) = NaN;

	% 2. NaN filter E values for outliers ( E_apparent_stats{1,7} )
	E_filt_Matrix(isoutlier(E_Matrix)) = NaN;

		% 3. NaN "E_apparent_filt_med" if over half is bad
		if sum(isnan(E_filt_Matrix),'all') > (numel(E_filt_Matrix)/2)
			E_apparent_filt_med = NaN;
        else
            E_apparent_filt_med = nanmedian(E_filt_Matrix, 'all');
		end


else
%%% Do this if no rsq_Matrix exists
% 1. NaN filter E values for outliers ( E_apparent_stats{1,7} )
E_filt_Matrix(isoutlier(E_Matrix)) = NaN;
		
	% 2. NaN "E_apparent_filt_med" if over half is bad
		if sum(isnan(E_filt_Matrix),'all') > (numel(E_filt_Matrix)/2)
			E_apparent_filt_med = NaN;
        else
            E_apparent_filt_med = nanmedian(E_filt_Matrix, 'all');
		end


end


%%% Percentage_nan %%%
% Calculate final percentage of NaNs for force map
Percentage_nan = sum(isnan(E_filt_Matrix),'all') / numel(E_filt_Matrix);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SAVING %%%%%%%%%%%%%
save(FileName, 'E_apparent_med', 'E_apparent_stats', 'E_filt_Matrix', 'E_apparent_filt_med', 'Percentage_nan', 'FileName', '-append')
%   Appends:
%       1. E_apparent_med (initial median apparent modulus)
%       2. E_apparent_stats (initial descriptive statistics of E_Matrix)
%       3. E_filt_Matrix
%       4. E_apparent_filt_med (median apparent modulus after filtering)
%       5. % of NaNs
%       6. FileName
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%{
%%%%% NOTE FOR EXISTING LIN_HERTZ DELETE "E_appt_med" and "E_filt"
% Use this for deleting specific vars from existing .mat files
Files=dir('*.mat*');
for k=1:length(Files)
   FileName=Files(k).name
% load file
load(FileName)
clear E_appt_med E_filt name
save(FileName, '-regexp', '^(?!(Files)$).')
%%%%%
end
%}

