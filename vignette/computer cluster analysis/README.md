## Analysis on computer cluster

The files included here can be used to run a BayesTwin analysis on a computer cluster and save the workspace. 

The accompanying unix script can be directly used to run the analysis on the gentic computer cluster of LISA surfsara (https://ctg.cncr.nl/research/the_genetic_cluster_computer), but you might need to slightly 
adjust it when you want to run the analysis on another system. 

In order to get the script working, you first have to make sure that BayesTwin and all dependencies are installed on the computer cluster 
as well as is JAGS. In the accompanying unix script, also the path to the folder where R libraries are installed is included (you might need to adjust this to the folder you are using). 

Furthermore, you need to use the same name for the R script and the folder the R script is included (here both are named Bayestwin).
If you use another name, you need to adjust the unix script (e.g., instead of using file=bayestwin use file=**the name you are using**)

