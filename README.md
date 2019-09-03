# hawkesinf

Current development version of Hawkes process inference.

Directory aaai18 contains the code to run the experiments in 

Christian R. Shelton, Zhen Qin, and Chandini Shetty
Hawkes Process Inference with Missing Data
AAAI 2018
http://www.cs.ucr.edu/~cshelton/papers/index.cgi?SheQinShe18


This code allows for inference over hidden labels:  data in which for specified periods of time and labels, the events are unobserved (via MCMC sampling).

The code *also* provides maximum likelihood estimation of parameters (with missing data, via EM).
