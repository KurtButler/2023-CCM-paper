# On Causal Discovery with Convergent Cross Mapping
In this repository, we provide code to reproduce the figures from our paper "[On Causal Discovery with Convergent Cross Mapping](https://doi.org/10.1109/TSP.2023.3286529)," published to IEEE Transactions on Signal Processing. 

> **Abstract:** Convergent cross mapping (CCM) is a principled causal discovery technique for signals, but its efficacy depends on a number of assumptions about the systems that generated the signals. We present a self-contained introduction to the theory of causality in state-space models, Takens’ theorem, and cross maps, and we propose conditions to check if a signal is appropriate for cross mapping. Further, we propose simple analyses based on Gaussian processes to test for these conditions in data. We show that our proposed techniques detect when convergent cross mapping may conclude erroneous results using several examples from the literature, and we comment on other considerations that are important when applying methods such as CCM.

This repo contains several functions that implement useful algorithms for studying nonlinear systems.
When the main file is run, it creates several files corresponding to tables and figures from the paper:
- Figure 2: Illustration of the effect of the parameter tau on state-space reconstruction (SSR)
- Figure 4: Basic demonstration of CCM
- Figure 5: Surrogate signal test with a deterministic and a random signal
- Figure 6: Comparison of the distance, covariance, correlation and recurrence matrices
- Figure 7: Demonstration of the recurrence principle
- Figure 8: Histogram of CCM convergence coefficients for a Rossler-Lorenz system.
- Table I: Comparison of CCM results and the proposed test statistics for several test cases.

## Download
You can download the code using git:
```
git clone https://github.com/KurtButler/2023-CCM-paper.git
```
Alternatively, you can view this code as a capsule on [Code Ocean](https://codeocean.com/capsule/8338092/tree/v1).

## Instructions
To generate all figures (as .png files), you just need to run `main.m`. The code should run with no issues using Matlab 2022a or later. All generated figures and tables will be saved to the results folder. 

If you wish to reproduce our results that used real data, you will need to download the data from the source. See `INSTRUCTIONS.md` in the data folder to see what you need.

## Dataset
In one of our experiments, we used the [ElectricityLoadDiagrams20112014](https://doi.org/10.24432/C58C86) data set from the UCI Machine Learning Repository. The data was donated by Artur Trindade at the University of Porto, Portugal.

## Citation
If you use any code or results from this project, please cite the orignal paper:
```
@article{butler2023causal,
  title={On Causal Discovery with Convergent Cross Mapping},
  author={Butler, Kurt and Feng, Guanchao and Djuri{\'c}, Petar M},
  journal={IEEE Transactions on Signal Processing},
  year={2023},
  volume={71},
  number={},
  pages={2595-2607},
  doi={10.1109/TSP.2023.3286529},
  publisher={IEEE}
}
```
