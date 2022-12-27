%%%% E_filt Loop %%%%
% After analyzing raw data from AFM_POST_JON and getting modulus arrays
% Use this loop to get median apparent_modulus for cell
% This also gives an option to filter based on rsq values & valid half-min
clear
for i = 0:4   %loop through all cells
    
name=sprintf('d1_WT09_0%s.mat',num2str(i))  %name extension from directory
load(name);
   
    
E_appt_med = median(E_Matrix(:)) %total apparent modulus (median)
if sum(sum((rsq_Matrix>0.90)))>2 E_filt=median(E_Matrix(rsq_Matrix>0.95)), else  E_filt=NaN, end
%here is a 90 percent R2 filter


clear ans
clear i
save(name);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% After running E_filt Loop
%put all E_filt values into one array
j=1
for i = 5:6   %loop through all cells

    
name=sprintf('AngTh_0%s.mat',num2str(i))  %name extension from directory
load(name);
   
AngTh_E_filt_arry(j) = E_filt

j=j+1
    

end

