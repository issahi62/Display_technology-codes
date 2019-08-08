% Display Tehnology 

% Ibrahim Issah 

%******************************
% First add all the folders and subfolders to the path
% the lab_2019.mat file contains the values obtained from the
% spectraphotometer.
% Issah_Ibrahim_main codes was used to assign these values to the required
% cells 
% main_for_running codes are used for the analysis using the screen_type
% command in the implemented codes developed by me. 
% The screen_type function basically uses most of the Piotr function with
% some changes I made in some of the functions. 
%******************************

%%
%*****************
% Codes for assigning the data into the required cells
%*****************
  
clc; 
load( 'lab_2019.mat')


lambda = 380:2:780; % required wavelength range
cmf = xyz31_1nm(21:2:421,2:4);  %usage of the color matching function xyz of the selected wavelengths.
newwarmup = interp1(Warmup(:, 1), Warmup(:, 2:end), lambda); %using interpolation command to get the required value range
newwarmup = newwarmup'; % transpose for the data. 

 
XYZ  = CalculateXYZ(newwarmup, cmf); % used Piotr CalculateXYZ function. 
time = 0:2:30; % timestep

Warmup_data = [time', XYZ(:, 2)]; % time versus the luminosity of the warmup section. 

Dell24{1,1} = Warmup_data; % placing it in the required cell

%plot(time, XYZ(:,2));

% section 2 for the uniformity
% thus how the luminous is uniformed on the screen

U = interp1(Uniformity(:,1), Uniformity(:, 2:end), lambda); % same interpolation command 
XYZUniformity = CalculateXYZ(U', cmf); %using the calculateXYZ and the color matching function.
xy = XYZtoxy(XYZUniformity); %using XYZ to xy function to get the chromaticity coordinates 
Uniformity = [XYZUniformity(:, 2), xy(:, 1), xy(:,2)]; % storing in the required cells 
Dell24{1,2} = Uniformity; % sectioning in the cell

% viewing angle
% same procedure
V = interp1(Viewingangle(:,1), Viewingangle(:, 2:end), lambda); 
XYZViewangle = CalculateXYZ(V', cmf); 
xy1 = XYZtoxy(XYZViewangle);
viewangle = [XYZViewangle(:, 2), xy1(:, 1), xy1(:,2)];
Dell24{1,3} = viewangle; 

%RGB_ramps section using the interpolation command

RGB_Ramps = [Red(:, 1), Red(:, 2:end), Green(:, 2:end), Blue(:, 2:end)];
RGB_a = interp1(RGB_Ramps(:,1), RGB_Ramps(:, 2:end), lambda);
RGB_a = RGB_a'; 
Dell24{1,4} = RGB_a; 

% color patches of the first 6 colors 

Color_a = interp1(RGB(:,1), RGB(:, 2:7), lambda); 
ColorPatches = Color_a'; 
Dell24{1,5} = ColorPatches; 

%%
%%**************************
%** Remaining codes for the analysis is 
%** found in the screen_type function which I have added to the implemented
%codes of Piotr's codes. 
%**** I basically used most of Piotr's function and added my function to
%inherit from Piotr's functions. 


% The main_for_running file is for running the various screen_types using
% the screen type_function.
%%***************************


%%


