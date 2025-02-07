# Network-regression-and-control-chart-analysis-for-Trenord
Code associated to the thesis "Network regression and centrality estimation for TRENORD's OD matrices" authored by Giuseppe Cicerone and supervised by Francesca Ieva and Piercesare Secchi.
## Overview
This repository contains the code and resources associated with the research paper titled **Network regression and centrality estimation for TRENORD's OD matrices** authored by Giuseppe Cicerone and supervised by Piercesare Secchi and Francesca Ieva.

In particular:
* Script 1: [supercent_method.R](https://github.com/Giuseppecicerone12/Network-regression-and-control-chart-analysis-for-Trenord/blob/main/supercent_method.R) applies the supercent method to the OD matrices obtained using the code of Adrien Pirin in his report: **ESTIMATION OF DYNAMIC ORIGIN-DESTINATION MATRICES INTEGRATING PARTIAL AUTOMATED PASSENGER COUNTING AND TICKETS SALES DATA**.
* Script 2: [control_charts.R](https://github.com/Giuseppecicerone12/Network-regression-and-control-chart-analysis-for-Trenord/blob/main/control_charts.R) analyzes the OD matrices obtained with Adrien Pirin's code(https://github.com/adri-pi/od-estimation) using statistical quality control tools.
* Function: [dynamic_network_analysis.R](https://github.com/Giuseppecicerone12/Network-regression-and-control-chart-analysis-for-Trenord/blob/main/dynamic_network_analysis.R)) analyzes the OD matrices obtained with Adrien Pirin's code using network analysis and functional data analysis tools. 
The code is developed in R (version 4.4.1), is it mandatory to run Adrien Pirin's code(https://github.com/adri-pi/od-estimation) before running these files.
