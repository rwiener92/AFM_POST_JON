%% Runtime parameters related to File I/O, graphics, etc

h5_file_loc ="D:\Microsopy\AFM\Patterns_Smiti_April_2021\210423\Dten00.h5";
%h5_file_loc = "D:\Microsopy\AFM\Patterns_Smiti_April_2021\cell_00200.h5";
PLOT_OPT =  0; % 1 makes plot, 0 doesn't.
FontSize = 10;
% If the data is saved, it will save all force curves, all maps, r^2
% values, as either single/double array. The force curves are stored into a
% cell array.
SAVE_OPT = 1 ; % 1 saves, 0 doesn't
SAVE_NAME = 'Dten00.mat';

% A note on the saved results
% F_Matrix (cell array) : contains the Force of deflection of cantilever (To make Force
% vs Depth)
% D_Matrix (cell array) : contains the indentation depth vectors for each indentation
% E_Matrix (double array) : contains last value (deepest value) of the point wise modulus 
% Ext_Matrix (cell array) : vector of the raw extension values for each
% indentation
% ExtDefl_Matrix (cell array) : vector of the deflection values for each
% indentation
%CP_Matrix (double array) : contact points
% Height_Matrix : visualization of relative heights based on contact points
% PWE_Matrix (cell array) : contains vector of pointwise modulus for each
% indentationj




%% Parameters for finding the contact point.

CONTACT_METHOD_OPT = 1;
%1: Least square fit to linear-quadratic piecewise function
%2: Uses ratio of variance method

%%%%%%%%%%%%%---------- PARAMETERS FOR LINEAR?QUADRATIC FIT METHOD --------
% Sets the max number of points to be considered in fitting. The data
% to be fit will start from NUM_PTS_CONSIDERED before the end of the
% extension.
NUM_PTS_CONSIDERED = 3000;
% Sets the maximum amount of deflection to be considered in nm. The
% goal is to have enough to capture the first part of the indentation.
MAX_DEFL_FIT = 20; %nm
% Approximate baseline deflection. Use the first ten points
% of the data used for fitting, where there should be no contact.
NUM_PTS_TO_AVERAGE = 10;
% If the standard deviation in the data used to calculate the baseline
% deflection is above this threshold, an error will be thrown.
MAX_STD_RAISE_ERROR = 1; %nm

%%%%%%%%%%%%----------- PARAMETERS FOR RATIO OF VARIANCE METHOD
%Rov is based on the ratio of variance for an interval before and after a
%given point. This sets the number of samples used in that interval.
ROV_INTERVAL_N = 10;

%% Parameters for AFM tip - Note no need to input spring constant. It is read from the experimentally determined value.
% https://www.nanoandmore.com/AFM-Probe-hq-xsc11-hard-al-bs - pattern
% cell. These probes are long enough to touch the cell.
R = 20; % Radius of curvature (nm) https://www.nanoandmore.com/AFM-Probe-PNP-TR
th = 20*pi/180; % Cone semi-angle
b = R*cos(th); % Cylindrical radius %%double check
v = .5; % Poisson ratio

%% Parameters for modulus calculation


MODEL_QUADRATIC_FIT = 0;
% This parameter allows you to fit the Force-Depth curve.
% 0: Pointwise, uses raw data
% 1: Pointwise, Quadratic fit
% 2: Pointwise and quadratic fit (compared plots)
% 3: Hertz contact (single E value), no fit