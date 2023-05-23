%% Matlab main file
% Title:   On Causal Discovery with Convergent Cross Mapping
% Author:   Kurt Butler
% Date:     May 23, 2022
% Description:  Reproduces all listed Figures and Tables from the paper, 
% and configures the MATLAB path to spit out files.


%% Set up the path
addpath('./Functions')
addpath('./Results')
addpath('./Scripts')


%% Figures
% Figure 4: Basic demonstration of CCM.

% Figure 5: Comparison of the surrogate data test for a deterministic signal vs a random one.

% Figure 6: Comparison of distance, covariance, correlation and recurrence matrices.

% Figure 7: Demonstration of the recurrence principle with a pair of sinusoids.

% Figure 8: Histogram of CCM convergence coefficients for a Rossler-Lorenz system.
Figure_8


%% Table I
% This will return the table data as text in the command window.
System_L1
System_L2
System_B1
System_B2
System_K1
System_K2
