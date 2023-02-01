%%%%%%
% This script takes all the cell data and aggregates their important
% information into a cell array "**_Matrix" and saves
%%%%%%

clear; clc
% Change all d1_WT09
Files=dir('d1_WT09*.mat');

%initialize matrix labels
d1_WT09_Matrix{1,1} = 'FileName'; 
d1_WT09_Matrix{2,1} = 'E_apparent_stats'; 
d1_WT09_Matrix{3,1} = 'E_apparent_med'; 
d1_WT09_Matrix{4,1} = 'E_filt_Matrix'; 
d1_WT09_Matrix{5,1} = 'E_apparent_filt_med'; 
d1_WT09_Matrix{6,1} = 'Percentage_nan'; 

for k=1:length(Files)
   FileName=Files(k).name;
% load file
load(FileName)

%store each cells data in aggregated matrix
d1_WT09_Matrix{1,k+1} = FileName;
d1_WT09_Matrix{2,k+1} = E_apparent_stats;
d1_WT09_Matrix{3,k+1} = E_apparent_med;
d1_WT09_Matrix{4,k+1} = E_filt_Matrix;
d1_WT09_Matrix{5,k+1} = E_apparent_filt_med;
d1_WT09_Matrix{6,k+1} = Percentage_nan;


%%%%%
%save d1_WT09_Matrix
save('d1_WT09_Matrix', '-regexp', '^(?!(j|k)$).') %change "Files" to detele another var
%%%%%
end

% reshape all calculated E indentations into single array for group
load('d1_WT09_Matrix.mat')
d1_WT09_Matrix{7,1} = 'E_filt_individual';
d1_WT09_Matrix{7,2} = cell2mat(d1_WT09_Matrix(4,2:end)); 
d1_WT09_Matrix{7,2} = d1_WT09_Matrix{7,2}(:);
save('d1_WT09_Matrix')

% average apparent modulus of group
nanmean(cell2mat(d1_WT09_Matrix(5,2:end)))
