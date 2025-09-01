# Freely Moving Hunger and Thirst (FMHT)

[![Creative Commons License](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by/4.0/)
[![DOI](https://img.shields.io/badge/DOI-Coming%20Soon-orange)](https://github.com/)

This repository contains code and data for analyzing neural recordings from mice engaged in naturalistic foraging behavior, investigating the neural representations underlying drive states (hunger/thirst) and actions.

## Repository Structure
### [`example_data/`](./example_data/)
this contains behavioral and electrophysiological data from a representative session of simultaneous foraging and neuropixels recording along a trajectory spanning hippocampus, thalamus, and hypothalamus saved as a .mat file

### 📊 [`primary_analysis/`](./primary_analysis/)
MATLAB code for analyzing and visualizing behavioral and electrophysiological data collected during the freely moving hunger and thirst paradigm.

### 🔄 [`unitmatch_analysis/`](./unitmatch_analysis/)
MATLAB implementation extending the [UnitMatch framework](https://github.com/EnnyvanBeest/UnitMatch) to track individual neurons across multiple recording days.

### ⚡ [`ephys_preprocessing/`](./ephys_preprocessing/)
Python utilities for preprocessing electrophysiological data from neuropixels recordings

### 🧠 [`anatomy/`](./anatomy/)
Python code for registering neural recording sites to the Allen Brain Atlas

## Citation

If you use this code or data in your research, please cite our paper: insert citation here future Lucas
