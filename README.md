# OTtracking
This is the repository contains the code for the paper Exploring Relevant Features for EEG-Based Investigation of Sound Perception in Naturalistic Soundscapes

the main.m script contains the order in which the scripts have to be executed. Some scripts require manual definition of paths and feature specifications. The corresponding data set can be found here:
https://zenodo.org/records/7147701

## Abstract
A comprehensive analysis of everyday sound perception can be achieved using Electroencephalography (EEG) with the concurrent acquisition of information about the environment.
While extensive research has been dedicated to speech perception, the complexities of auditory perception within everyday environments, specifically the types of information and the key features to extract, remain less explored. Our study aims to systematically investigate the relevance of different feature categories: discrete sound-identity markers, general cognitive state information, and acoustic representations, including discrete sound onset, the envelope, and mel-spectrogram. 

Using continuous data analysis, we contrast different methods in terms of their predictive power for unseen data, their distinct contributions to explaining neural data. We also evaluate the results considering the impact of sound context, here the density of acoustic events. For this, we analyse data from a complex audio-visual motor task using a naturalistic soundscape.

The results demonstrated that model prediction is increased using more acoustically detailed features in conjunction with a comprehensive description of the sound identity of the soundscape. Crucially, the outcome hinged on excluding periods devoid of sound onsets in the analysis in the case of the discrete features. Furthermore, we showed that considering the sound event density was crucial when using discrete acoustic onsets. 

Our study highlights the importance of a comprehensive description of the soundscape, using acoustic and non-acoustic aspects, to fully understand the dynamics of sound perception in complex situations. 
This approach can serve as a foundation for future studies aiming to investigate sound perception in natural settings.

#### Toolbox dependencies

mTRF toolbox (https://github.com/mickcrosse/mTRF-Toolbox)

eeglab - https://sccn.ucsd.edu/eeglab/index.php

field trip - https://www.fieldtriptoolbox.org/

violinplot - https://github.com/bastibe/Violinplot-Matlab

sigstar - https://de.mathworks.com/matlabcentral/fileexchange/39696-raacampbell-sigstar

BBCI - https://github.com/bbci/bbci_public


