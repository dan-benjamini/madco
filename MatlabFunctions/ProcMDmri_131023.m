clc
clear 
close all

%% initialize

im_first_ind= [firstIndSample1, firstIndSample2]; % the experiment number of the first 1D-T2 experiment

pathnamesDATA = {'PathToData/Sample1','PathToData/Sample2'}; % a path to the experimental data. e.g., /Users/benjaminidh/Documents/RESEARCH/DATA/300_NICHD/NICHD_Bruker7T/RB1604_0120.WZ1

titles={'NameOfSample1','NameOfSample2'}; %names/numbers of the samples. e.g., "Sarm345"

threshIm=[threshSample1, threshSample2]; %threshold for mask. Should be a threshold value based on ExpData.ImageRef. e.g., threshIm=1e5

%% read the MD-MRI data
[ExpPars, ExpData] = readMDmri(im_first_ind,pathnamesDATA);


%% filter and vectorize the MD-MRI data
[ExpDataF,ExpDataFvect] = filterMDmri(ExpData,threshIm);

%add titles
ExpDataFvect.titles=titles;


%%
% sampleNum contains the samples to process, 
sampleNum=[1, 2];

% sliceNum is a cell array that contains vectors with the slice numbers to process
% within each sample, refers to slice "z" in Image(x,y,z). 
sliceNum={20:30, 23:30};

% pathToSave is a path to a local directory on your computer to save the
% proc data. It will automatically create the directory, and
% sub-directories with the names in the "title" variable 
pathToSave= '/PathToProc';

MultiDprocFun(sampleNum,sliceNum,pathToSave,ExpDataFvect,ExpPars);




