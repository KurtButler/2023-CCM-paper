%% Matlab main file
% Title:   On Causal Discovery with Convergent Cross Mapping
% Author:   Kurt Butler
% Date:     June 17, 2023
% Description:  Reproduces all listed Figures and Tables from the paper, and configures the MATLAB path to spit out files.

% Notes:
% - Figures 1 and 3 were deemed too unimportant to include in this repository, even though we indeed used code in their preparation. 
% - I didnt record which random seed that we used when doing simulations for the paper, so this code may produce results that look mildly different from what appears in the final paper. Nonetheless, the results were similar enough from realization to realization that we felt this was O.K.


%% Set up the path
addpath('./Functions')
addpath('./Results')
addpath('./Scripts')



%% Turn this guy off
% Sometimes, the surrogate signal analyses could trigger warnings from Matlab while training a GP for autoprediction with a random signal. I turned off warnings here to suppress these messages from flooding the log file.
warning('off','all')


%% Table I
% This will return the table data as text in the 'output' file
disp('      <> <> <> TABLE I <> <> <> ')
System_L1
System_L2
System_B1
System_B2
System_K1
System_K2
System_K3
System_K4
% This system requires the Electricity Loads data set. See INSTRUCTIONS.md in the data folder.
if exist('LD2011_2014.txt','file')
  System_E1;
else
  disp('Electricity data set skipped.')
end
disp('Table I done.')


%% Figures
% Figure 2: Varying the tau parameter in SSR
Figure_2
disp('Figure 2 done.')

% Figure 4: Basic demonstration of CCM.
Figure_4
disp('Figure 4 done.')

% Figure 5: Comparison of the surrogate data test for a deterministic signal vs a random one.
Figure_5
disp('Figure 5 done.')

% Figure 6: Recurrence matrices.
Figure_6
disp('Figure 6 done.')

% Figure 7: Demonstration of the recurrence principle
Figure_7
disp('Figure 7 done.')

% Figure 8: Histogram of CCM convergence coefficients for a Rossler-Lorenz system.
Figure_8
disp('Figure 8 done.')





disp('Everything is done.')

