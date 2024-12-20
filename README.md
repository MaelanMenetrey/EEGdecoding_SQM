# EEGdecoding_SQM
In this repository we present the code and data to reproduce the results shown in Menétrey, Herzog & Pascucci (under review). Neural correlates of unconscious and conscious visual processing.
Preprint: https://www.biorxiv.org/content/10.1101/2024.11.15.623853v1

# 1. System requirements
Software, toolbox and functions required:
- MATLAB (MathWorks Inc., Natick, MA, USA): The scripts were executed using MATLAB R2022b.
- EEGLAB: The scripts were executed using EEGLAB v2021.1, which can be downloaded from the official EEGLAB website (https://sccn.ucsd.edu/eeglab/download.php).
- Custom-made functions written in MATLAB, which can be downloaded in this repository (EEGdecoding_SQM/Functions/...)

# 2. Installation guide
To reproduce the results shown in the manuscript, you need to download all the files of this repository, keeping the same folder structure. In the Data folder, you need to add the original data (preprocessed EEG data and behavioral data) that can be downloaded here: https://osf.io/d83vs/
Then, you need to change accordingly the paths in each analysis file to indicate the location of the Data folder (main) and add to the path EEGLAB, and the Function folder.

- 1.PrepareDataForAnalyses/PrepareEEGdata.m: Prepare the data for the decoding analyses. Keep onyl the valid trial (removing bad EEG epochs, and trials where the RT was too short). EEG data are also resampled, re-epoched (-0.2 to 1s) and zscored. The data that are obtained with this script can also directly be downloaded in (EEGdecoding_SQM/Data/...)

  

