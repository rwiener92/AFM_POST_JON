%%%%%%
% This script takes all the cell data and aggregates their important
% information into a cell array "**_Matrix" and saves
%%%%%%

clear; clc
% Change all AngTh
Files=dir('AngTh*.mat');

%initialize matrix labels
AngTh_Matrix{1,1} = 'FileName'; 
AngTh_Matrix{2,1} = 'E_apparent_stats'; 
AngTh_Matrix{3,1} = 'E_apparent_med'; 
AngTh_Matrix{4,1} = 'E_filt_Matrix'; 
AngTh_Matrix{5,1} = 'E_apparent_filt_med'; 
AngTh_Matrix{6,1} = 'Percentage_nan'; 

for k=1:length(Files)
   FileName=Files(k).name;
% load file
load(FileName)

%store each cells data in aggregated matrix
AngTh_Matrix{1,k+1} = FileName;
AngTh_Matrix{2,k+1} = E_apparent_stats;
AngTh_Matrix{3,k+1} = E_apparent_med;
AngTh_Matrix{4,k+1} = E_filt_Matrix;
AngTh_Matrix{5,k+1} = E_apparent_filt_med;
AngTh_Matrix{6,k+1} = Percentage_nan;


%%%%%
save AngTh_Matrix
%%%%%
end

% average apparent modulus of group
mean(cell2mat(AngTh_Matrix(5,2:end)))
