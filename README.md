Supplementary information / reproducible research files for the manuscript 
Title: "Simultaneous Inference for Multiple Comparisons in
Overdispersed Multinomial Data"

Authors: Budig, S., Vogel, C., Schaarschmidt, F.
Code was written by Budig, S. and Vogel, C.
In case of questions or comments please contact budig@cell.uni-hannover.de!

The code was written/evaluated/run in R with the following software versions:
R version 4.3.3 (2024-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

The folder "Code_and_Data" contains the following data and files that can be used to reproduce 
all analyses and figures of the manuscript.
It contains three subfolders with the following files:

./Code/:
This folder contains the Code for the Simulation Study and for the evaluation of the examples 
presented in the manuscript.

        ./Example_Analysis/:
        This folder contains the R file "Example_Analysis.R" to evaulate the three examples.
        The subfolders contain the corresponding datasets

        ./Simulation_Study/:
        This folder contains all r files which are necessary to reproduce the results of the simulation study
        The mult_sim_main.R file is used to run the main simulation study. The Code_and_Data folder must be set as the path. 
        The number of simulations should be adjusted accordingly. The Dirichlet-multinomial model is 
        relatively computationally intensive and can drastically increase the time required.
        The mult_sim_addzero.R file is used to run the additional simulation study in which a one is added 
        to a random cluster within a group. The Code_and_Data folder must be set as the path. 
        The number of simulations should be adjusted accordingly.
        After the simulation has been run, the file mult_sim_data_preprocessing.R is used to 
        process the results and summarise them in an .rds file
        The mult_sim_graphs.R can then be used to reproduce the graphics from the manuscript 
        with the final result file.

./Figures/:
This folder contains all the graphics shown in the manuscript and produced by the given R-code

./Results/:
This folder contains the processed results from the intermediate results of the simulation file. 
Due to the size of the intermediate results, they were not uploaded to Github

    ./main/:
    This subfolder contains the processed results from the main simulation study

    ./zero_handling/:
    This subfolder contains the processed results of the simulation with the added ones




