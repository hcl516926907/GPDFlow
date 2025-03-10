# GPDFlow: Generative Multivariate Threshold Exceedance Modeling via Normalizing Flows

This repository contains the code and resources necessary to reproduce the results presented in the paper:

**"GPDFlow: Generative Multivariate Threshold Exceedance Modeling via Normalizing Flows"**

GPDFlow is a generative modeling approach designed for threshold exceedance modeling by integrating the Multivariate Generalized Pareto Distribution (mGPD) with normalizing flows. It is coded in Python by Pytorch and the normflows package.

## Repository Structure

### Directories

- **`Data_and_Model`**:
  - Contains the dataset used in the application section of the paper.
  - Contains the fitted weights for the GPDFlow model.

- **`KRSW`**:
  - Contains R code for fitting the classic parametric Multivariate Generalized Pareto Distribution (mGPD).
  - For detailed methodology, refer to the supporting materials of the paper *"Peaks Over Thresholds Modeling With Multivariate Generalized Pareto Distributions"*.

### Files

- **`Application.ipynb`**:
  - Jupyter notebook containing code to reproduce results from the application section of the paper.

- **`GPDFlow_Simulation_Scenario_1.ipynb`**:
  - Jupyter notebook for reproducing results of Simulation Scenario 1 described in the paper.

- **`GPDFlow_Simulation_Scenario_2.ipynb`**:
  - Jupyter notebook for reproducing results of Simulation Scenario 2 described in the paper.

- **`GPDFlow.py`**:
  - Python script providing the implementation and setup of the GPDFlow model.

- **`mGPD_Estimation.R`**:
  - R script containing the code to fit the classic Multivariate Generalized Pareto Distribution (mGPD).

## How to Use

Clone the repository and navigate through the notebooks and scripts to replicate experiments, simulations, and applications as described in the paper.

Please ensure you have Python, R, and required libraries installed to run the notebooks and scripts successfully.


