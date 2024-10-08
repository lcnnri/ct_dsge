# [ct_dsge](https://lcnnri.github.io/ct_dsge/)
Replication package for the paper "Estimation of continuous-time linear DSGE models from discrete-time measurements" by *B.J. Christensen, L. Neri and J.C. Parra-Alvarez*

`README.md`: This file contains a description of the replication files for the empirical application in Section 6 of “Estimation of Continuous-time Linear DSGE Models from Discrete-time Measurements”  by *Bent Jesper Christensen, Luca Neri, and Juan Carlos Parra-Alvarez* (2024), Journal of Econometrics, forthcoming. 

Copyright ©, 2024, Bent Jesper Christensen, Luca Neri, and Juan Carlos Parra-Alvarez. All Rights Reserved.

The code is provided free of charge. If you use the code, please cite the publication. 			

All files have been tested in MATLAB Release: R2024a.

This version: September, 2024. 

For updates, please check the authors' webpages: 
- [Luca Neri](https://www.lneri.com)
- [Juan Carlos Parra-Alvarez](https://jcparra-alvarez.weebly.com)

For additional information, please contact the corresponding author: *Bent Jesper Christensen*, Department of Economics and Business Economics, Aarhus University, Fuglesangs Alle 4, 8210 Aarhus V, Denmark, [bjchristensen@econ.au.dk](mailto:bjchristensen@econ.au.dk).

***
## Programs
1) [`main.m`](matlab/main.m): replicates Tables 4, 5, and 7 in Section 6 using the configurations described in the paper in one single run.

2) main_table_*.m: replicates the specific table, Table *, in Section 6. 
		user's choice:
			Panel = 'a' for replicating panel a) (left)
			Panel = 'b' for replicating panel b) (right)
	- [`main_table_4.m`](matlab/main_table_4.m)
	- [`main_table_5.m`](matlab/main_table_5.m)
 	- [`main_table_7.m`](matlab/main_table_7.m)	


3) [`check_identification_condition.m`](matlab/check_identification_condition.m): checks if the rank and order conditions in Proposition 2 of Section 3 are satisfied. The file calls the function [isIdentified.m](utils/isIdentified.m) to check if θ is locally identified.
		user's choice:
```matlab
choose_model = 1` % {1 2 3} = {F-SSR, S-SSR, EM-SSR}
choose_condition = 'a' % {'a','b','c'}, i.e., the conditions in Proposition 2
``` 
### Example:

```matlab
	[idFlag, OCflag] = ... % identification, cf. Proposition 2
        		isIdentified(theta_hat,cfg,y,choose_condition);
```
 - where `theta_hat` is the vector of estimated and calibrated parameters, `cfg` is the configuration settings, and `choose_condition` selects the condition (a, b, or c) of Proposition 2.


***
## Folders
The replication files in main_table_*.m and check_identification_condition.m automatically add the following required folders to Matlab's path:

1) [utils/](matlab/utils/): set of routines for estimation of the model. 

2) [Linearization/](matlab/Linearization/): set of routines to compute the solution of the linear DSGE model. These routines were written by SeHyoun Ahn, Greg Kaplan, Benjamin Moll, Thomas Winberry, and Christian Wolf for the article "When Inequality Matters for Macro and Macro Matters for Inequality," NBER Macroeconomics Annual 2017, volume 32. University of Chicago Press, 2017.

3) [model_loglin/](matlab/model_loglin/): folder that contains the files related to the model specification and solution. 

	3.1) [`paramstructbase.m`](matlab/model_loglin/paramstructbase.m): set the model's parameter values.

	3.2) [`equilibrium_conditions_loglin.m`](matlab/model_loglin/equilibrium_conditions_loglin.m): set the model's equilibrium conditions.

	3.3) [`model_solution_loglin.m`](matlab/model_loglin/model_solution_loglin.m): calls the routines in the folder Linearisation to compute the log-linear approximation to the model's solution, 	and writes the solution in the state space form of Section 2 (see Equations (2.1) and (2.2)).

4) [data/](matlab/data/): folder containing the time series used in the estimation routines. Please see Section 6 for a description of the variables used. 

## Instructions
1. Clone this repository.
2. Run the scripts in the **/matlab** folder to reproduce the results.

## License
This work is licensed under the [MIT License](LICENSE)
 
