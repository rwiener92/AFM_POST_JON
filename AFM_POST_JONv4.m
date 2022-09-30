% This script processes the h5 file made by using Asylum's .ARDF converter
%
% J.H. 2020


clear
clc
close all
warning off
%% File locations, reading data, setting up strings
AFM_INPUT


%Load the h5 into matlab and read it.
info = h5info(h5_file_loc);
info2 = h5info(h5_file_loc,'/ForceMap');
%Get the strings labels for the force curves.
mdata = info2.Groups.Datasets;
%Combine them into a cell arra of size N+1 curves x 1 (the plus one is for
%a dataset called segment, which contains index of where the data changes
%from extension to retraction)
FCnames = {mdata.Name}';

%Convert the string array FCnames into numbers
FCnames_str = string(FCnames(1:end-1));
sndx = strfind(FCnames_str,':');
sndx = [sndx{:}]';
rndx = zeros(length(FCnames_str),1);
cndx = rndx;
for j = 1:length(FCnames_str)
    this_string = char(FCnames_str(j));
    rndx(j) = str2num(this_string(1:sndx(j)-1));
    cndx(j) =  str2num(this_string(sndx(j)+1:end));
end

% read in spring constant from metadeta

%read in the spring constant from metadata
mdata = info.Attributes.Value;
mdata = strsplit(mdata,'\n');
mdata = string(mdata);
TF = contains(mdata,'SpringConstant: 0.');
mdata = mdata(TF);
string_w_sc = char(mdata(1));
pndx = strfind(string_w_sc,'.');
sc = string_w_sc(pndx-1:end);
spring_constant = str2num(sc); %N/m or nN/nM

%-----------------------------------------------------------------
%% Looping over force curves
% Gather the index for which the curves switches from ext to ret to away
%the segment dataset is at the end.
segmentdata = h5read(h5_file_loc,['/ForceMap/0/' FCnames{end}]);

%initialize a matrix to hold all the info for the forcemap and contact points.
nRows = max(rndx)+1;
nCols = max(cndx)+1;
E_Matrix = zeros(nRows,nCols);
CP_Matrix = E_Matrix; % contact points
if MODEL_QUADRATIC_FIT == 2
    E_Matrix_raw = E_Matrix;
elseif MODEL_QUADRATIC_FIT == 3
    rsq_Matrix = E_Matrix;
end

if SAVE_OPT == 1
    F_Matrix = cell(nRows,nCols);
    D_Matrix = F_Matrix;
    RoV_Matrix = F_Matrix; % ratio of variance
    RoVZ_Matrix = F_Matrix; %ratio of variance Z sensor
    Ext_Matrix = F_Matrix;%extension matrix used for fitting
    ExtDefl_Matrix = F_Matrix;% z sensor of extension matrix
    if MODEL_QUADRATIC_FIT == 0
        PWE_Matrix = F_Matrix;
    end
end
%length-1 comes from not wanting to loop over Segments dataset at the end.
tic

% read in spring constant from the meta-deta

%CHANGE THIS FOR TO PARFOR TO RUN FASTER
%CAN ONLY SET BREAKPOINTS WHEN USING FOR
parfor i = 1:nCols
    warning off % needed when using parfor, redundant when not
    thiscol = zeros(nRows,1);
    thiscol_CP = thiscol;
    %temp variables for parfor depending on model
    if MODEL_QUADRATIC_FIT == 2
        thiscol_raw = thiscol;
    elseif MODEL_QUADRATIC_FIT == 3
        thiscol_rsq = thiscol;
    end
    %temp variables for parfor if saving
    if SAVE_OPT == 1
        thiscol_F = cell(nRows,1);
        thiscol_D = thiscol_F;
        thiscol_RoV = thiscol_F;
        thiscol_RoVZ = thiscol_F;
        thiscol_Ext = thiscol_F;
        thiscol_ExtDefl = thiscol_F;
        if MODEL_QUADRATIC_FIT == 0
            thiscol_PWE = thiscol_F;
        end
    end
    for i2 = 1:nRows
        close all
        %for some unknown reason, sometimes a specific AFM indentation is
        %missing from the data. make sure it's there before proceeding with
        %analyis
        %
        this = string([num2str(i2-1) ':' num2str(i-1)]);
        if any(strcmp(FCnames_str,this))
            data = h5read(h5_file_loc,['/ForceMap/0/' num2str(i2-1) ':' num2str(i-1)]);
        else
            continue
        end
        % data is a npoints x 3 matrix
        % Channels (columns): Raw (s) , Defl (m) , ZSnsr (m)
        % Each channel contains values corresponding to the entire indentation
        % and retraction. The index corresponding to the end of each respective
        % phase is stored in 3 x n x m matrix in the the variable
        % segmentdata. However, I calculate it on my own, as I have found
        % it to not be reliable.
        raw = data(:,1);
        defl = data(:,2)/1e-9; %converted to nm
        Zsnsr = data(:,3)/1e-9; %converted to nm
        
        % Extract the two indices e.g., '2' and '3' from dataset '2' '3'
        
        ndx_1 = i2; %row ndx
        ndx_2 = i; % col ndx

        
        extndx = segmentdata(1,ndx_2,ndx_1);
        retndx = segmentdata(2,ndx_2,ndx_1);
        if PLOT_OPT == 1
            figure
            plot(Zsnsr(1:extndx),defl(1:extndx));
            title('Raw Data: Extension vs Deflection')
            xlabel('Extension (nm)')
            ylabel('Deflection (nm)')
        end
        

        
        start_idx = max(1,extndx-NUM_PTS_CONSIDERED+1);
        fit_ext = Zsnsr(start_idx:extndx); % +1 for off by one error
        fit_defl = defl(start_idx:extndx);
        
        
        
        if any( diff(fit_ext)<0 )
            error('The deflection curve contains some of the retraction portion.')
            %         else
            %             extndx = find( diff(Zsnsr)<0, 1);
            %             fit_ext = Zsnsr(extndx-NUM_PTS_CONSIDERED+1:extndx); % +1 for off by one error
            %             fit_defl = defl(extndx-NUM_PTS_CONSIDERED+1:extndx);
        end
        %% -------------------- Find the point of contact, xc -----------------
        if CONTACT_METHOD_OPT == 1
            
            
            baseline_defl = mean(fit_defl(1:NUM_PTS_TO_AVERAGE));
            baseline_std = std(fit_defl(1:NUM_PTS_TO_AVERAGE));
            
            
            if baseline_std > MAX_STD_RAISE_ERROR
                error(['The standard deviation of the baseline deflection data set' ...
                    ' is higher than specified maximum.'])
            end
            
            % Find the first place where deflection is MAX_DEFL_FIT more than this baseline
            fit_defl_temp = fit_defl-baseline_defl;
            cutoff_ndx=find(fit_defl_temp>MAX_DEFL_FIT,1,'first');
            if isempty(cutoff_ndx)
                deflstring = sprintf('Objective deflection not obtained for row %i col %i', ndx_1,ndx_2);
                disp(deflstring)
                cutoff_ndx = length(fit_ext);
                
            end
            fit_ext = fit_ext(1:cutoff_ndx);
            
            
            fit_defl = fit_defl(1:cutoff_ndx);
            
            % Create the index of the crosspoints to be considered and some
            % initializations
            xcs_fit = 1:length(fit_ext);
            errvec = zeros(length(xcs_fit),1);
            fits = cell(length(xcs_fit),1);
            %empty vector to store fits
            
            for j = 1:length(xcs_fit)
                %note that extension is in um
                [errvec(j),soln, fits{j} ] = lsqfitFC_Jon(fit_ext,fit_defl,xcs_fit(j));
            end
            
            % Get the minimum error xc
            [minerr,minxcndx_fit] = min(errvec);
            %convert the index back into context of entire extnesion vector
            %minxcndx = minxcndx_fit + extndx-NUM_PTS_CONSIDERED;
            %for the following line of code, if everything from the start
            %is kept, then extndx-NUM_PTS_CONSIDERED will fuck up indexing
            minxcndx = minxcndx_fit + start_idx-1;
            
            %Plotting check
            if PLOT_OPT == 1
                figure
                plot(fit_ext,fit_defl,'-*');
                hold all
                plot(fit_ext,fits{minxcndx_fit},'-');
                hold all
                plot(fit_ext(minxcndx_fit),fit_defl(minxcndx_fit),'ko','markersize',10)
                xlabel('Z sensor / Extension (nm)')
                title('Deflection Curve')
                ylabel('Deflection (nm)')
                legend('Raw','Fit','location','best')
            end
            
        elseif CONTACT_METHOD_OPT == 2
            % RATIO OF VARIANCE
            leftpts = zeros(ROV_INTERVAL_N,length(fit_defl));
            rightpts = leftpts;
            % The following loop uses circshift functions to allow for
            % vectorization when calculating the variances.
            for j = 1:ROV_INTERVAL_N
                leftpts(j,:) = circshift(fit_defl,j);
                rightpts(j,:) = circshift(fit_defl,-j);
            end
            leftvar = var(leftpts(:,ROV_INTERVAL_N+1:end-ROV_INTERVAL_N));
            rightvar = var(rightpts(:,ROV_INTERVAL_N+1:end-ROV_INTERVAL_N));
            ROV = rightvar ./ leftvar;
            
            [maxROV,minxcndx_fit] = max(ROV);
            minxcndx_fit = minxcndx_fit + ROV_INTERVAL_N; %account for offset
            minxcndx = minxcndx_fit + extndx-NUM_PTS_CONSIDERED;
            
            %Plotting check
            extplot = fit_ext(ROV_INTERVAL_N+1:end-ROV_INTERVAL_N);
            if length(extplot) < 2
                continue
            end
            if PLOT_OPT == 1
                figure
                subplot(2,1,1)
                plot(fit_ext,fit_defl,'-*');
                title('Deflection Curve')
                ylabel('Deflection (nm)')
                hold all
                plot(fit_ext(minxcndx_fit),fit_defl(minxcndx_fit),'ko','markersize',10)
                legend('Raw data','Contact Point','location','best')
                xlim([extplot(1) fit_ext(end)])
                set(gca,'fontsize',FontSize)
                
                subplot(2,1,2)
                plot(extplot,ROV)
                ylabel('ROV')
                xlim([extplot(1) fit_ext(end)]);
                xlabel('Z sensor / Extension (nm)')
                title('Ratio of Variance')
                set(gca,'fontsize',FontSize)
                
%                 %% this plot is a check to see that ROV is being calculated correctly
%                 ROVvec = zeros(size(fit_defl));
%                 for p = ROV_INTERVAL_N+1:length(ROVvec)-ROV_INTERVAL_N
%                    right = var(fit_defl(p+1:p+1+ROV_INTERVAL_N-1));
%                    left = var(fit_defl(p-ROV_INTERVAL_N:p-1));
%                    ROVvec(p) = right / left;
%                 end
%                 subplot(3,1,3)
%                 plot(fit_ext,ROVvec)
%                 ylabel('ROV')
%                 xlim([extplot(1) fit_ext(end)]);
%                 xlabel('Z sensor / Extension (nm)')
%                 title('Verification of ROV')
%                 xlabel('Z sensor / Extension (nm)')
            end
            
            
        end
        thiscol_CP(i2) = Zsnsr(minxcndx);
        
        %% ---------------------- Calculate the modulus ----------------------
        
        %calculate the index where the indentation starts and create Depth and
        %deflection vector, both in nm.
        minxcndx = minxcndx;
        
        
        
        def =  defl(minxcndx:extndx) - defl(minxcndx);  % h-h0
        D = Zsnsr(minxcndx:extndx) - Zsnsr(minxcndx); %z- z0
        D = D - def; % D = (z-z0) - (h-h0)
        F = def.*spring_constant; % force vector using Hooke's law
        [E_app,regimeChange] = calc_E_app(D,F,R,th,b,'pointwise',PLOT_OPT);
        
        if PLOT_OPT == 0
            %Plot depth vs deflection
            figure
            plot(D,def,'*-')
            xlabel('Depth (nm)')
            ylabel('Deflection (nm)')
            set(gca,'fontsize',FontSize)
        end
        
        if MODEL_QUADRATIC_FIT == 1 || MODEL_QUADRATIC_FIT == 2
            Draw = D;
            def_raw = def;
            Fraw = F;
            E_app_raw = E_app;
            regimeChange_raw = regimeChange;
            
            p = polyfit(D,def,2);
            D = linspace(D(1),max(D));
            def = polyval(p,D);
            F = def.*spring_constant; % force vector using Hooke's law
            [E_app,regimeChange] = calc_E_app(D,F,R,th,b,'pointwise',PLOT_OPT);
        end
        
        
        if MODEL_QUADRATIC_FIT == 1 || MODEL_QUADRATIC_FIT == 2
            [E_app,regimeChange] = calc_E_app(D,F,R,th,b,'pointwise');
        elseif MODEL_QUADRATIC_FIT == 3
            %             Draw = D;
            %             def_raw = def;
            %             Fraw = F;
            %             E_app_raw = E_app;
            %             p = polyfit(D,def,2);
            %             D = linspace(D(1),max(D));
            %             def = polyval(p,D);
            %F = def.*spring_constant; % force vector using Hooke's law
            [E_app,~,rsq] = calc_E_app(D,F,R,th,b,'Hertz',PLOT_OPT);
        end
        
        %% ---------------- PLOTTING E ---------------------------
        % E_app is in units of nN / nM^2
        E_app = E_app * 1e18 * 1e-9 ; %units of N/m or Pa
        E_app = E_app /1000; % units of kPa
        E = E_app .* 2 .* (1-v^2);
        if MODEL_QUADRATIC_FIT == 2
            E_app_raw = E_app_raw * 1e18 * 1e-9 ; %units of N/m or Pa
            E_app_raw = E_app_raw /1000; % units of kPa
            E_raw = E_app_raw .* 2 .* (1-v^2);
        end
        
        
        
        if PLOT_OPT == 1
            if MODEL_QUADRATIC_FIT == 0 || MODEL_QUADRATIC_FIT == 1
                figure
                subplot(2,1,1)
                plot(D,E,'*-')
                xlabel('Depth (nm)'); ylabel('E (kPa)')
                
                hold all
                if regimeChange ~= 0
                    plot(D(regimeChange),E(regimeChange),'o','markersize',10,'linewidth',2)
                    axis auto
                end
                set(gca,'fontsize',FontSize)
                title('Pointwise Elastic Modulus')
                legend('E','Conical Geometry Transition','location','best')
                ylim([0 inf])
                subplot(2,1,2)
                if MODEL_QUADRATIC_FIT == 1
                    plot(D(2:end),def(2:end),'r-')
                    hold on
                    plot(Draw,def_raw,'*');
                    legend('Fit','Raw','location','best')
                else
                    plot(D(2:end),def(2:end),'*-')
                end
                xlabel('Depth (nm)'); ylabel('Deflection (nm)')
                set(gca,'fontsize',FontSize)
                if regimeChange ~= 0
                    ylim([0 max(E(regimeChange:end))])
                end
                axis auto
                title('Force Curve (post-contact)')
            elseif MODEL_QUADRATIC_FIT == 2
                figure
                subplot(2,1,1)
                plot(D,E,'-')
                hold on; plot(Draw,E_raw,'r*');
                xlabel('Depth (nm)'); ylabel('E (kPa)')
                hold all
                if regimeChange ~= 0
                    plot(D(regimeChange),E(regimeChange),'o','markersize',10,'linewidth',2)
                    plot(Draw(regimeChange_raw),E_raw(regimeChange_raw),'o','markersize',10,'linewidth',2)
                    axis auto
                end
                set(gca,'fontsize',FontSize)
                title('Pointwise Elastic Modulus')
                legend('E_f_i_t','E_r_a_w','Conical Transition (fit)','Conical Transition (raw)','location','best')
                
                subplot(2,1,2)
                plot(D(2:end),def(2:end),'-')
                hold on
                plot(Draw,def_raw,'r*');
                legend('Fit','Raw','location','best')
                xlabel('Depth (nm)'); ylabel('Deflection (nm)')
                set(gca,'fontsize',FontSize)
                if regimeChange ~= 0
                    ylim([0 max(E(regimeChange:end))])
                end
                axis auto
                title('Force Curve (post-contact)')
            end
            % If any of the data is to be saved
        end
        
        thiscol(i2) = E(end);
        
        if MODEL_QUADRATIC_FIT == 2
            thiscol_raw(i2) = E_raw(end);
        elseif MODEL_QUADRATIC_FIT == 3
            thiscol_rsq(i2) = rsq;
        end
        
        if SAVE_OPT == 1
            thiscol_F{i2} = F;
            thiscol_D{i2} = D;
            thiscol_Ext{i2} = fit_ext;
            thiscol_ExtDefl{i2} = fit_defl;
            if MODEL_QUADRATIC_FIT == 0
                thiscol_PWE{i2} = real(E);
            end
            if CONTACT_METHOD_OPT == 2
                thiscol_RoV{i2} = ROV;
                thiscol_RoVZ{i2} = extplot;
            end
        end
        
    end
    E_Matrix(:,i) = thiscol;
    CP_Matrix(:,i) = thiscol_CP;
    if MODEL_QUADRATIC_FIT == 2
        E_Matrix_raw(:,i) = thiscol_raw;
    elseif MODEL_QUADRATIC_FIT == 3
        rsq_Matrix(:,i) = thiscol_rsq;
    end
    
    if SAVE_OPT == 1
        F_Matrix(:,i) = thiscol_F;
        D_Matrix(:,i) = thiscol_D;
        RoV_Matrix(:,i) = thiscol_RoV;
        RoVZ_Matrix(:,i) = thiscol_RoVZ;
        Ext_Matrix(:,i) = thiscol_Ext;
        ExtDefl_Matrix(:,i) = thiscol_ExtDefl;
        
        
        if MODEL_QUADRATIC_FIT == 0
            PWE_Matrix(:,i) = thiscol_PWE;
        end
    end
    
    s = '%.2f percent left to process.';
    s = sprintf(s,i/nCols*100);
    disp(s)
    
end
toc

%% Plotting maps

% Occasionally floating point errors will lead to numbers being considered
% as having an imaginary part.

E_Matrix = real(E_Matrix);
if MODEL_QUADRATIC_FIT == 0 || MODEL_QUADRATIC_FIT == 1
    figure
    imagesc(E_Matrix);
    title('Elastic Modulus')
    c = colorbar; c.Label.String = 'kPa';
    set(gca,'fontsize',FontSize)
    set(gca,'YDir','normal')
elseif MODEL_QUADRATIC_FIT == 2
    figure
    subplot(1,2,1)
    imagesc(E_Matrix)
    title('Elastic Modulus (fit)')
    c = colorbar; c.Label.String = 'kPa';
    set(gca,'fontsize',FontSize)
    set(gca,'YDir','normal')
    subplot(1,2,2)
    imagesc(E_Matrix_raw)
    title('Elastic Modulus (raw)')
    c = colorbar; c.Label.String = 'kPa';
    set(gca,'fontsize',FontSize)
    set(gca,'YDir','normal')
elseif MODEL_QUADRATIC_FIT == 3
    figure
    imagesc(E_Matrix);
    title('Elastic Modulus')
    c = colorbar; c.Label.String = 'kPa';
    set(gca,'fontsize',FontSize)
    set(gca,'YDir','normal')
    
    figure
    imagesc(rsq_Matrix);
    title('R^2')
    c = colorbar;
    set(gca,'fontsize',FontSize)
    set(gca,'YDir','normal')
    
end


% For Surface Height
Height_Matrix = max(CP_Matrix(:)) - CP_Matrix;
figure; imagesc(Height_Matrix);
title('Relative Height (CP estimate)')
c = colorbar; c.Label.String = 'nm';
set(gca,'fontsize',FontSize)
set(gca,'YDir','normal')

%% Saving
if SAVE_OPT == 1
    save(SAVE_NAME,'E_Matrix','F_Matrix','D_Matrix','CP_Matrix', ...
        'h5_file_loc','CONTACT_METHOD_OPT','R','th','b', ...
        'RoV_Matrix','RoVZ_Matrix','Ext_Matrix','ExtDefl_Matrix','Height_Matrix', ...
        'spring_constant','v');
    if MODEL_QUADRATIC_FIT == 0
        save(SAVE_NAME,'PWE_Matrix','-append');
    elseif MODEL_QUADRATIC_FIT == 3
        rsq_Matrix = real(rsq_Matrix);
        save(SAVE_NAME,'rsq_Matrix','-append')
    end
end

%%

warning on
