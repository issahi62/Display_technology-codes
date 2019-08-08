%***********************
%% Display Technologies 2019
%% Issah Ibrahim 

%***********************
%% Run these codes to attain  all the results required. 
%% All comments are in the screen_type function, the color_gamut codes is
%% for ploting the chromaticity diagram of all the displays in one plot. 
%% The plot_with_offset or without_offset function is for the choice of the channel chromatic constancy.
%% The combinegraph is for combing the individual graphs of the displays into one graph. 
%% The Issah_Ibrahim_main function was the manipulation of data in to require cell arrangement. 
%% The compared_graphs folder contains the images of all the compared graphs and the other folders 
%% contain their respective images.
%% The ImplementedCodes are functions created by Piotr Bartczak which I used in my data manipulations. 


%
% The screen_type function has two arguments
% 1. is the screen type you will want to analyse 
% 2. is wether you want offset or not 
% choose 1 for offset and 0 for without offset in the main_for_running
% command. 
%
% CREATE the required folders for saving your screen_type. 
% And kindly change the folder names for the required files.
% The codes below the screen_type(screen_name, offset) helps to save all
% the plots at their respective folders.

% The combinegraph function also combines the five display graphs into 1
% figure to compare the results. 
% This was done manually by editing the figures to the right screen_names. 
% codes for running the screen_type function
clear 
close all
clc
%% 
  % put here your own path
path='/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah';
addpath(genpath(path)) %adding subfolders and all folders. 

load( 'IN_LABORATORY_2019.mat', 'CRT', 'Dell24', 'Dell_Konica', 'EIZO_2012', 'Projector_DLP_BENQ', 'xyz31_1nm');
%screen_type(Dell24,1);
screen_type(CRT,0);
FolderName = '/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah/CRT_Figs';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end
%saveas(yaw, fullfile(fname, sprintf('phase_an.png')));

close all 

screen_type(Dell24,0);
FolderName = '/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah/Dell24figs';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end

close all 

screen_type(Dell_Konica,0);
FolderName = '/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah/Dell_konicafigs';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end
close all 

screen_type(EIZO_2012,0);
FolderName = '/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah/Eizo_2012';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end

close all 

screen_type(Projector_DLP_BENQ,0);
FolderName = '/Users/kobbyTilly/Desktop/DT_Ibrahim_Issah/Projectfigs';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end

figure; 
color_gamut; % plot the color gamut of displays. 



%comment%
%%
%A physical model for a CRT display involves two steps: 
%a nonlinear transformation between drive voltages in the digital-to-analog
%converter (DAC- values) and corresponding phosphor radiance levels, and
%a linear transformation between phosphor radiance level and corresponding 
%XYZ tristimulus values? 
%In the case of CRT displays, the Gain-Offset- Gamma model (GOG-model) 
%is normally used to describe the nonlinear part of the transformation \
%The coefficients kg, k0 and ? represent the model?s g
%ain offset and nonlinearity, respectively. N
%is the DAC quantization depth. dr , dg and db are DAC-values for R, G and B
%%

 %This difference is due to the short term instability of the Eizo monitor, 
 %thus, the maximum luminance of the phosphors change within a spam of a 
 %single day causing differences between the training data set and test data set
