# XB-sol-cycle-var-manuscript

A fixed repository with all the codes and data needed to create Figures 1-7 of the manuscript

# Requirements

You will need either R or RStudio and certain libraries to run the codes.  We recommend that you install RStudio
which will promt you to install the required libraries.

After cloning the repository:

%cd code

Start RStudio and now you can generate each figure using:

>source('make_figX.R') where X =1,2,..,7

Please note that some of the scripts (make_fig3.R, make_fig4.R and make_fig5.R) take a few minutes to run. 
The runtime of these scripts is determined by the number of bootstraps which is set to 500.  
If you would like to speed-up the scripts you can change this number to 200 (or 100) but note that this does
mean that you will see more variation in the results between runs. 
 
For comments or questions please email us at: pete@predsci.com, mbennun@predsci.com 
