# DecisionsOnNetworks
Simulation software to study consensus collective decisions on random networks with synchronous and asynchronous decision models.

This code allows you to generate all data and results presented in the article

Andreagiovanni Reina, Thomas Bose, Vaibhav Srivastava, and James A. R. Marshall, "Asynchrony rescues statistically-optimal group decisions from information cascades through emergent leaders" _The Royal Society Open Science_ (under review)


## Dependencies

The Python code has been run with Python 3.6.8 and has the following dependencies from external libraries:
* configparser
* json
* math
* numpy
* matplotlib
* scipy
* networkx
* itertools

## Generate the data


#### Figure 2

In order to reproduce the results reported in Figure 2 of the paper, that is, synchronous collective signal detection, you need to run the script
```
DecisionsOnNetworks/scripts/submitSynch.sh
```
which will exectute a large number of simulations for all conisidered conditions.
On line [160 of `submitSynch.sh`](https://github.com/DiODeProject/DecisionsOnNetworks/blob/published_code/scripts/submitSynch.sh#L160) you can change the command to execute the Python code using the method that you prefer (default it is by submitting SLURM jobs).

If you wish to submit a single simulation run, for a specific condition, for the synchronous collective signal detection, you can do it by running the command:
```
cd DecisionsOnNetworks/src
python3 DecNet/DecisionProcess.py ../conf/DecNet.config
```
You can modify any parameters in the configuration file `DecisionsOnNetworks/conf/DecNet.config`.
Once the process had terminated, you can find the results in the folder `DecisionsOnNetworks/data`, which will include the text file `out.txt` (with the data used to generate the paper's Figures), and a graphical representation of each iteration step of each simulation as pdf files. 


#### Figure 3

In order to reproduce the results reported in Figure 3 of the paper, that is, asynchronous collective sequential sampling, you need to run the script
```
DecisionsOnNetworks/scripts/submitAsynch.sh
```
which will exectute a large number of simulations for all conisidered conditions.
On line [139 of `submitAsynch.sh`](https://github.com/DiODeProject/DecisionsOnNetworks/blob/published_code/scripts/submitAsynch.sh#L139) you can change the command to execute the Python code using the method that you prefer (default it is by submitting SLURM jobs).

If you wish to submit a single simulation run, for a specific condition, for the asynchronous collective sequential sampling, you can do it by running the command:
```
cd DecisionsOnNetworks/src
python3 AsynchKicks/DecisionProcess.py ../conf/AsynchK.config
```
You can modify any parameters in the configuration file `DecisionsOnNetworks/conf/AsynchK.config`.
Once the process had terminated, you can find the results in the folder `DecisionsOnNetworks/data`, which will include text files `out.txt` with the information of the decision of the agents at the end of each simulation, the text file `out_cas.txt` with the cascade size of each simulation, and a graphical representation of each iteration step for each simulation as pdf files. 

## Plot the data

#### Figure 1

The graphics of Figure 1 have been generated using the Mathematica notebook `DecisionsOnNetworks/plots/paperPlots.nb`.

#### Figure 2 and 3 and Supplementary Figures SF1, SF2, SF3, and SF6

The data generated with the Python simulator have been plotted using the Mathematica notebook `DecisionsOnNetworks/plots/paperPlots.nb`.

#### Figure 3(b) and Supplementary Figures SF4 and SF5

For Figure 3(b) and Supplementary Figures SF4 and SF5, the data generated with the Python simulator have been plotted using the R script `DecisionsOnNetworks/plots/plotKicks.R`, calling the function `plotCascadeOnRank()`.


