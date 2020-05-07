1. h_neg_point_1 is a vector of h's for all neurons and times.

2. Brain_K_pair is a two column array of pairwise K's.  First column true.  Second Column Est.

3. Brain_K_self is the self couplings.  Same format as above.

4. Brain_time_series is a three column array of the dynamics of our model.  You should
exclude the first 100 time points, as the system needed time to equilibrate

Column 1 is simulated from the true LM, column 2 from the MF (no linear response), and 
column 3 is after linear response.


5.  Brain_Temp_mean and var show how the mean and variance of synchrony (the time-dependent magnetization)
change with the inverse temperature, beta.  Each of these files has 4 columns
First, a vector of the beta values used.  Second, using the true LM.  Third, using Linear Response,
and Fourth, using just mean-field.

6.  MCMC.m simulates the Ising model for a particular set of Lagrange multipliers.

7.  MCMC_MaxCal_Ising.m contains the script to replicate our findings.
	******The last section of this code takes hours to run.  As an alternative,
	Brain_Temp_mean.txt and Brain_Temp_var.txt provide the output of this last 
	section only.