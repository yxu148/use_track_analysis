# use_track_analysis
This folder is to use the Matlab-Track-Analysis, or analyze data after collecting.
Here is a brief introduction to some useful codes, for more detailed information, please refer to the source codes, which have detailed annotation inside.
All the codes below should be ran from your folder. For example in my case, run from \user specific\Yiming.

## processMMF.m 
processMMF.m will process the MMF videos and output BIN files.
You need to edit the directories (startdir1, outputdir1, and logpath) so the codes know where is the MMF video, and where to save
the extracted data and logs.

## processBIN.m
processBIN.m will process the BIN files and output MAT files.
You need to edit the directories (startdir1, outputdir1, and logpath) so the codes know where is the BIN files, and where to save
the extracted data and logs.

## load_data.m, load_multi_data.m
Analyze MAT files. I use load_data.m to try analyzing one experiment and use load_multi_data.m to analyze multiple experiments.
load_multi_data.m will make a data structure data.mat, which saves important features for each larva and will be needed for the other plots in the furture.
**Together with these codes, all the codes below need to adjust the base directory according to the directory of your extracted data.
All of these codes with, for example '_gr21a' in the end is specifically for the data of Gr21a larvae.**

## plot_from_data.m
Add more properties to data.mat and plot figures based on data.mat

## info_from_pturn.m
Plot figures related to turn probability based on data.mat

## track_trace.m
Modified version of Isabel's code to generate video of larvae moving with track ID labeled.

## Other files are not useful for you to run.