# On Causal Discovery with Convergent Cross Mapping
In this repo, we provide code to reproduce the figures from our paper "On Causal Discovery with Convergent Cross Mapping," published to IEEE Transactions on Signal Processing. 

> **Abstract:** Convergent cross mapping is a principled causal discovery technique for signals, but its efficacy depends on a number of assumptions about the systems that generated the signals. We present a self-contained introduction to the theory of causality in state-space models, Takens’ theorem, and cross maps, and we propose conditions to check if a signal is appropriate for cross mapping. Further, we propose simple analyses based on Gaussian processes to test for these conditions in data. We show that our proposed techniques detect when convergent cross mapping may conclude erroneous results using several examples from the literature, and we comment on other considerations that are important when applying methods such as CCM.

This repo contains several functions that implement useful algorithms for studying nonlinear systems.
When the main file is run, it creates several files corresponding to tables and figures from the paper:
- Figure 4: Basic demonstration of CCM.
- Figure 5: Comparison of the surrogate data test for a deterministic signal vs a random one.
- Figure 8: Histogram of CCM convergence coefficients for a Rossler-Lorenz system.
- Table I: Comparison of CCM results and the proposed test statistics for several test cases.

## Instructions
To generate all figures (as .png files), you just need to run `main.m`. Also all the subfiles can be run independently.

We wrote the code with Matlab 2022a.

## Citation
If you use this project, please cite:
```
@article{butler2023causal,
  title={On Causal Discovery with Convergent Cross Mapping},
  author={Butler, Kurt and Feng, Guanchao and Djuri{\'c}, Petar M},
  journal={IEEE Transactions on Signal Processing},
  year={2023},
  publisher={IEEE}
}
```
## Code
You can download the code using git:
```
git clone https://github.com/KurtButler/2023-CCM-paper.git
```
