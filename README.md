This ROOT/RooFit code performs a statistical analysis of vector boson fusion (VBF) data from CMS experiments in 2017 and 2018. It fits background and signal models to datasets, computes significance and p-value, and visualizes the results.

Key Steps:

Read Data: Loads data from 2017CMSVBF.txt and 2018CMSVBF.txt.

Model Creation: Defines background (exponential) and signal (Gaussian) models for each year.

Simultaneous Fit: Combines both years' data into a simultaneous model using RooSimultaneous.

Fit and Plot: Performs the fit, plots the data and models, and calculates the significance and p-value.

Results: Saves the significance and p-value for further analysis.

Output: Generates plots and saves results in a text file (results/mass.txt).
