function screen_type(screen_name, offset)

%% Display Technology 2019 
%% Ibrahim Issah
clc; 
load( 'IN_LABORATORY_2019.mat', 'CRT', 'Dell24', 'Dell_Konica', 'EIZO_2012', 'Projector_DLP_BENQ', 'xyz31_1nm');

w_range = 380:2:780;
cmf = xyz31_1nm(21:2:421,2:4); 
%*******************
%% Exercise 1
%*******************
figure(1);
plot(screen_name{1}(:,1),screen_name{1}(:,2))  % plot time VS Luminance
title('stabilization curve')
xlabel('Time [t]')
ylabel('Luminosity [cd]');
hold on; 
% plot tolerence 10 -+ 1 
plot(screen_name{1}([1,end],1),[screen_name{1}(end,2)*0.01+ screen_name{1}(end,2),screen_name{1}(end,2)*0.01+screen_name{1}(end,2)],'r'); 
plot(screen_name{1}([1,end],1),[screen_name{1}(end,2)-screen_name{1}(end,2)*0.01,screen_name{1}(end,2)-screen_name{1}(end,2)*0.01],'r');
hold off;
%%
%**************************
% COMMENTS
%**************************
% The above graphs is a plot of teh luminance of the varied displays as a
% function of time. 1% tolerance was added relative to the last luminance
% values. 
%***************************
% The CRT display showed a well stabilized after 25mins with an initial
% high illuminance which drops as a function of time depicting the effect of
% heat and stresses due to the high metal halide lamp at background. 
%***************************
% The Dell24 display had a sharp increase in the luminance around 0.5mins
% and remain constance for about 30mins depicting the stability of the
% luminance as a function of time. 
% ***************************
% The Dell Konica also is relatively similar to the Dell24 display.
% However, there was a sharp increase in luminance around 0.25mins which
% was seemingly over saturated yet started decreasing as a function of time
% and remain constant around 25mins. 
%****************************
%The EIZ0 display is relatively simular to the CRT display yet the
%luminance descreased systematically as a function of time. The EIZO is
%quite unstable as a function of time
%***************************
% The projector also showed a gradual decrease in luminance as a function
% of time with its stabilization in between the +1% and -1% tolerance
% relative to the luminance values.The display was quite stabilized around
% 6mins till 25mins and showed to vary gradually.
% *************************
% These variation in luminances relative to the display depicts how each
% display backlight in relation to the phospors utilized is affected by
% heat or stress variations in each display. 
%%

%**********
%% EXERCISE 2
%**********
%%
% color ramps 
figure(2);
% 0:15:255 for each Primary 18 levels
%***************
%Red primary
plot(w_range,screen_name{4}(1:18,:),'r'), hold on  
% Green primary
plot(w_range,screen_name{4}(19:36,:),'g'), hold on 
% Blue primary
plot(380:2:780,screen_name{4}(37:54,:),'b'),
xlabel('Wavelength [nm]')
% Depends on the device which was used while measurements
ylabel('Radiance [W m^{-2} sr^{-1}]')
title('screen_name primary color ramps')
%**************************
%COMMENTS 
%**************************
%The CRT display has good color ramps for the red, blue and green. The blue
%has relatively low radiance of around 0.002W/m2/sr yet higher than the green
%color ramps at a step size of 0:15:255. However, the red had higher peaks
%around 650nm and 700nm respectively.There is some small flourescence of
%the green primary mixed with the red primary around 600nm. 
%*****************************
%The Dell24 display showed relatively higher peak of radiance comparable to
%the other screens. The blue ramp showed a peak around 480nm with a
%narrow-band flourescence due to perharps its flourescent backlights.
%*****************************
%The Dell Konica display showed relatively low radiance for the blue
%ramps comparable to the Dell24. In addition to the high radiance of the
%green primary and the red primary there is some small flourescence around
%600nm. However, this display showed less peak of the red primary around
%650nm to 700nm.This flourescence visualization also depict the probability
%of the display have flourescent backlight. 
%*****************************
%The EIZ0 display showed relatively similar radiance to the blue ramps of
%the Dell Konica display. However, the green color ramp showed lower
%radiance between 500nm to 550nm comparable to the Dell Konica display with
%a higher peak respectively. 
%******************************
%The projector also shows good color ramps for the varied primaries with
%relatively gradual increase in radiance as the specific primary level
%increases from 0:15:255.They have relatively broad band for the blue,
%green and red color ramps comparable to the other displays. Even there is
%some small dependencies in the ramps between 480nm to 500nm and 580nm to
%600nm the ramps are quite spatially dependent.
%*****************************
%The Dell Konica and Dell24 displays have higher peak radiance in 
%all three channels than the conventional display due to probably, 
%the use of the narrow-band fluorescent backlights.

%%
%****************************
%% Exercise 3
%****************************
%color RGB patches
figure(3); 
plot(w_range,screen_name{5}(1,:),'r'),hold on
%
white=screen_name{4}(18,:)+screen_name{4}(36,:)+screen_name{4}(54,:)-2*screen_name{4}(1,:);
plot(w_range,white,'b')
xlabel('\lambda [nm]');
ylabel('Radiance [W m^{-2} sr^{-1}]');
title('screen name Channel independence')

%*************************
%COMMENTS 
%*************************
%Testing additivity is the prediction of the white color as a sum of the
%individually measured primarier which is compared with the meaure white. 
% The CRT and EIZ0 showed relatively good additivity testing which showed relatively
%good monochrome of each primaries. 
%**********************
%The small failure of additivity for the Dell24 and Dell Konica display might
%well be due to a small increase in flare at the high luminance levels.
%There may also be circuitry on board the display which increases power 
%sent to the electron guns to compensate for their increased load. 
%******************
%The projector showed large differences between white and R+G+B which is quite disturbing.
%A possible cause may be the strong angular dependency of this display
% and color shift errors.
%******************
%Notably, the degree of additivity is sufficient to justify 
%the use of a 3x3 primary matrix transform.
%****************************
%%

%*************************
%% EXERCISE 4 - color matching function
%*************************
xyz31_1nm(1,1); 
xyz31_1nm(21,1); 
xyz31_1nm(421,1); 
ckobby = xyz31_1nm(21:2:421, 2:4); 
figure(4); 
plot(w_range, ckobby, 'LineWidth', 2.5); 
xlabel('wavelength$(\lambda)$', 'Interpreter', 'latex'); 
ylabel("$V'(\lambda)$", 'Interpreter', 'latex'); 
%************************
%COMMENT
% **********************
% The 2 degrees V(lamda) curve. showing the color matching function
% relative to photopic vision. 
%%

%*************************
%% EXERCISE 5 - VECTOR ALGEBRA
%*************************

A = [1,2,3]; 
B = [-2,4,3]';
[~]=A*B; 
[~]=B*A;
[~]=A.*A; 



%%
% Chromacity diagram
%%
data = eye(201); % spd of monochromatic sources
data(202,:) = data(1,:); % for creating a closed curve
%%
%*****************************
%Exercise 6 - tistimulus values
%****************************
XYZ = CalculateXYZ(data,cmf);% XYZ of the locus points
%***********************
%COMMENT
%***********************
%CalculateXYZ is a function for calculating the XYZ tristimuls values for
%its argument. Thus, the data and the color matching functions. 
%%

%***************************
%% EXERCISE 7 - Contrast ratio
%***************************
% contrast calculation
black=CalculateXYZ(screen_name{4}(1,:),cmf);
white=CalculateXYZ(screen_name{5}(1,:),cmf);
contrast =white(2)./black(2);
disp("contrast value");
disp(contrast); 
%*************
%COMMENT
%*************
% This is the ratio of the luminance of the brightest color to the
% luminance of the darkest color. 
% 1. Dell24 - 682:1
% 2. EIZ0 - 256:1
% 3. Dell Konica - 97:1
% 4. CRT - 66.13:1
% 5. Projector - 37:1
% The results show that projector has the lowest constrat which is relative
% due to broad screen and effect of ambient illumination. Additionally the
% DMD mirrors also accepts and rejects some of the illumination. 
%*******************
%%

%********************
%% EXERCISE 8 - CIE CHROMATICITY COORDINATES
%********************
xy = XYZtoxy(XYZ); % xy of the locus points, calculate points representing the spectral locus

%*****************
% The XYZtoxy function calculates the XYZ tristimulus values calculated to
% chromaticity coordinates x,y. 
%*****************
%%
%************************
%% Exercise 9 - chromaticity diagram
%************************
figure(5); 
plot(xy(:,1),xy(:,2),'g');
axis equal; grid on;
xlabel('x');
ylabel('y');
title('Chromaticity Diagram');
%*********************
%COMMENTS
%*********************
% The graph uses the chromaticity coordinates to plot the chroma diagram.
%I have created another codes known as color_gamut which plots all the
%displays on one chromaticity diagram. 
%From the results
% EIZ0 has the highest color gamut which is comparable to adobe RGB
% chromaticity coordinates. Dell24 and Dell Konica have relatively the same
% color gamut which quite comparable to sRGB color gamut. However the CRT
% and Projector has the lowest color gamut compared to the other displays. 
% The CRT however is relatively beetter than the projector in terms of the
% amount of colors each display can produce. 
%%
hold on

%Standard sRGB and aRGB 
sRGB =[0.64 0.33;0.3 0.6;0.15 0.06];
sRGB(4,:) = sRGB(1,:); % close the curve
aRGB = [0.64 0.33; 0.21 0.71; 0.15 0.06];
aRGB(4,:) = aRGB(1,:); % close the curve
plot(sRGB(:,1),sRGB(:,2),'b-o');
plot(aRGB(:,1),aRGB(:,2),'r-x');

XYZ=CalculateXYZ(screen_name{4}([18 36 54 18] ,:),cmf); % maximum level of each primary-color is at Row nr. 18 36 54 in the cell nr 4. 
% We use [18 36 54 18] in order to create closed curve
xy=XYZtoxy(XYZ);
% subplot(1,3,1);
plot(xy(:,1),xy(:,2),'k'); 

%*************************
%% Exercise 11 -channel chromatic constancy
%*************************

if offset ==1
    plot_with_offset(screen_name, xyz31_1nm);
else
    plot_without_offset(screen_name, xyz31_1nm);
end 
%*******************
%COMMENTS 
%*******************
% The chromaticity coordinates of different excitation levels were plotted
% on the chromaticity. The two functions plot_with_offset and plot_without
% offset. From these two functions it was visualized that after substrating
% the offset values the chromaticity coordinates points which was scattered
% on the graph shifted to the a point on their primariies which is red,
% green and blue. This shows how offset values affect the additivity
% testing and thereby should be removed to form the approximately the same
% color of white when you combine the primaries.
% Although the CRT and EIZO displays showed this shift of chroma points,
% the Dell24, Dell Konica and Projector showed some small offset values
% still reappering on the graph depicting the effect of failure in the
% channel indepency or additive testing of these displays. 
%%
%%
%****************************
%% Exercise 12 -CIELAB
%****************************


tmp=screen_name{2}; % screen_name 2 order is Yxy
tmp(:,4)=tmp(:,1);   
tmp=tmp(:,2:4);
screen_name_lab = zeros(25,3); 
for i=1:25 
%13= center measurement 
screen_name_lab(i,:)=XYZtoLab(xyYtoXYZ(tmp(i,:)),CalculateXYZ(screen_name{5}(1,:),cmf)); 
% In this case we use different transformation to Lab space where L corresponds to Lightness
end
%**********************
%Comments
%********************
%XYZtoLab function converts the tristimulus values to CIELAB color
%coordinates to determine the color difference between two color stimuli. 
%%

%%
% Function dEcalc - calculates color differences E*ab between the reference
% and other measurement
% compared with reference patch.
dE_screen_name=zeros(1,25);
for i=1:25 
dE_screen_name(i,:)=dEcalc(screen_name_lab(i,:),screen_name_lab(13,:));  
end
%*******************
%Comments 
%*********************
%The dECalc function calculates the color differnece of two color stimuli
%Thus, the color difference between the reference middle point and the rest
%%
%*************
%EXERCISE 13- SPATIAL UNIFORMITY
%**************
% show result in 5x5 array

Screen=[dE_screen_name(1:5);dE_screen_name(6:10);dE_screen_name(11:15);dE_screen_name(16:20);dE_screen_name(21:25)];  
% Now you see the difference between middle and other points (higher value means higher difference)
 figure(6); 
 surf(Screen);
 colormap(hot);
 colorbar('location','southoutside');
%********************************
%COMMENTS 
%*******************************
% From the 3D graphs it could be inferred that color difference at
% different positions on the screen varies relative to the display type. 
% The CRT colors are the edges are quite different from the reference white
% depicting less spatial uniformity in terms of how we percieve colors at
% different positions on the screen. CRT has an average color difference of
% 18.5 which is relatively good compared to the reference stimuli. 
% The Dell24 and Dell Konica display showed relatively closer
%color difference from the reference color with Dell24 having a less color
%difference relative to the Dell Konica display.Therefore in terms of
%spatial uniformity Dell24 will be quite better than Dell Konica. 
%However, the EIZO showed (NaN) for the color difference
% which could inferred it has good spatial uniformity.
% The projector also showed the highest average sum of color differences
%from the reference white which is comparable to the Dell Konica display.
% However, CRT is having resonably good spatial uniformity.
%%
%show sum
Kabun =mean(sum(dE_screen_name));
disp('The value color difference in screen_type is :');
disp(Kabun); 
%*****************
%Comment 
%*****************
% The average of these color difference for each display was calculated 
% EIZO - NaN
% CRT - 18.5
% Dell24 - 48.5
% Dell Konica - 71.5
% Projector - 91.9
%%
%**************************
%% EXERCISE 14
%**************************
% Effect of viewing angle should be done in CIE Lab coordinates

% The y AND Y order are different in our CRT structure therefore we
% have to rearrange them for correct calculations

tmp2=screen_name{3}; % CRT 2 order is Yxy
tmp2(:,4)=tmp2(:,1);
% Now the order is xyY to match input in the xyYtoXYZ function
tmp2=tmp2(:,2:4); 
angle_view = zeros(11,3); 
for i=1:11   
angle_view(i,:)=XYZtoLab(xyYtoXYZ(tmp2(i,:)),CalculateXYZ(screen_name{5}(1,:),cmf)); 
% In this case we use different transformation to Lab space where L corresponds to Lightness
end


% Function dEcalc - calculates color differences E*ab between the
% reference and other measurement


% Compare results with the reference white at 0 degree
angle_diff = zeros(1,11); 
for i=1:11
angle_diff(i)=dEcalc(angle_view(i,:),angle_view(6,:));
end
kobbybryan = -75:15:75; 
figure(7)
plot(kobbybryan,angle_diff)
legend('screentype')
xlabel('Angle')
ylabel('Angular Difference')
%*********************
%COMMENTS 
%*********************
% From the graphs it could be inferred that the CRT display has quite
%fluctuations in angle difference relative to the viewing angle of the
%display. Which means that display may have fluctuations in displaying
%colors per solid angle of the viewer. 
%********************
% The Dell24 and Dell Konica display showed similar quite similar viewing
% angle difference with a a nonlinear relations of the neutral gray sample
% at different solid angles. As the angle varies from the 0 degrees the
% differences in color a varies exponentially.
% From the graphs the fluctuations maybe due to human error or the incorrect
% positioning of the measuring device.
%***************
% However the projector also showed a huge variations in the percieved
% colors based on the solid angle. The graph was quite similar band stop
% filter and becomes constant in some relative angle of view. 
%The EIZ0 showed to having constant difference in angle from -30 to 30
%before there was a sharp increase in color difference from the neutral
%grey. 
% This depicts that EIZ0 is like an AMOLED with good viewing angle and high
% contrast. 

%% Nonlinear relationship between luminance and DAC values
% prepare R G and B ramps
r=screen_name{4}(1:18,:);  % assign all red ramps to r variable
g=screen_name{4}(19:36,:);
b=screen_name{4}(37:54,:);
xx= 0:15:255;
% use equation (3)
figure(8)
SumR=sum(r,2)./sum(r(18,:));    % sum(r(1,:))  % Scale each ramp by its maximum digit (its primary), then in all cases maximum Channel values is 1
plot(xx,SumR,'*')
xlabel('d r');ylabel('R')
hold on; 
figure(9)
SumG=sum(g,2)./sum(g(18,:));
plot(xx,SumG,'+')
xlabel('d g');ylabel('G')

figure(10)
SumB=sum(b,2)./sum(b(18,:));
plot(xx,SumB,'o')
xlabel('d b');ylabel('B')

figure(11) % Plot all channels in the same Figure
plot(xx,SumR,'r*'),hold on
plot(xx,SumG,'g+')
plot(xx,SumB,'bo')
xlabel('d r g b');ylabel('R G B')
title('Nonlinear relationship between luminance and DAC values')
%%
%**********************
%COMMENTS 
%**********************
% The varied displays showed the non-linearity between the DAC values
% selected for the colors and its luminance values. This nonlinearity leads
% effect varied color differences errors. 
% The GOG model which is used to model the tone curve characteristics of
% CRT monitors was used to fit the other displays. This basically considers
% the output luminance of a display with a certain input digital value
% respectively. 
% The GoG model is also used to fit relatively small amount of measuremet
% comparing to the interpolation fitting using the Look Up table approach. 
N=8; % because we have 8 bits 0-255
GOy=fitting_GOG(xx,SumR',[1,1,1]);     % SumR=sum(r,2)./sum(r(18,:)); 
tmp=(GOy(1)*(xx/(2^N-1))+GOy(2)).^GOy(3);
figure(12)
subplot(1,3,1)
plot(xx,real(SumR),'o',xx,real(tmp),'r');
title('Red channel')
% % % for the Red channel    (model's Gain, offset, nonlinearity (Gamma))
% % %  ROy =  0.5274    0.4752    3.8792    

% Blue

N=8;
GOy=fitting_GOG(xx,SumB',[2,.055,2]);
tmp=(GOy(1)*(xx/(2^N-1))+GOy(2)).^GOy(3);
subplot(1,3,2)
plot(xx,SumB,'o',xx,tmp,'b');
title('Blue channel');
% Green

N=8;
GOy=fitting_GOG(0:15:255,SumG',[2,.055,2]);
tmp=(GOy(1)*(xx/(2^N-1))+GOy(2)).^GOy(3);
subplot(1,3,3)
plot(xx,SumG,'o',xx,tmp,'g');
title('Green channel');
%%

% The interpolation of 1D-lookup tables (LUTs)
offset=screen_name{4}(1,:);
r=screen_name{4}(1:18,:);
g=screen_name{4}(19:36,:);
b=screen_name{4}(37:54,:);

% Usually the linear part of the transformation includes a separate
% offset term and hence offset can  be subtracted from the nonlinear part

r=r-ones(18,1)*offset;  % substract offset from each ramp  (18 x offset for 18 patches)
g=g-ones(18,1)*offset;
b=b-ones(18,1)*offset;

R=sum(r,2)./sum(r(18,:)); %Normalization
G=sum(g,2)./sum(g(18,:));
B=sum(b,2)./sum(b(18,:));
figure(13)
plot(xx,R,'ro',xx,G,'go',xx,B,'bo')
title("RGBramps");
%***********
% COMMENTS 
%***********
% This depicts the non-linearity between the DAC values after substracting
% the flare or offset from it to reduce the errors. 


%%

% Calculate DAC-values dr , dg , db that correspond to scalar values R = 0.2, G = 0.7,
% and B = 0.5.

% Calculate DAC-values dr
figure(14)
R_dr=0.2;
y0=interp1(R,xx,R_dr);

plot(0:15:255,R,'ro')
hold on
plot(y0,R_dr,'bx')
title('Red channel- Interpolation')
xlabel('Digital count')
hold off
%% comment 
% This section uses lookup table to interpolate the DAC values of some
% specific RGB values to determine if we could get the same luminance with
% those specific DAC values obtained form the fitting. 
figure(15)
% Calculate DAC-values dg
G_dg=0.7;
y0=interp1(G,xx,G_dg,'spline');  % Finding the corresponding value by using interpolation method
% y0=interp1(G,xx,G_dg)
plot(xx,G,'go')
hold on
plot(y0,G_dg,'rx')
xlabel('Digital count')
title('Green channel- Interpolation')
hold off

figure(16)
% Calculate DAC-values db
B_db=0.5;
y0=interp1(B,xx,B_db);
plot(xx,B,'bo')
hold on
plot(y0,B_db,'rx')
xlabel('Digital count')
title('Blue channel- Interpolation')
hold off

%%
% compare the arbitrary colors
test=[[255,255,255];[69,69,69];[24,124,224];[200,0,200];[50,50,100];[18,31,129]]; %Selected Color Patches allocated in the cell nr 5
% Color matching functions  380:2:780   of course it depends on the device used while measurements

%Primaries and offset for screen_name display
offset=screen_name{4}(1,:);

%Calculate XYZ from the spectral data !!! Results are obtained from our measurements
%testXYZ=CalculateXYZ(screen_name{5}(6,:),cmf);  % 4.7721    3.0587   19.0303  XYZ

% As an input we use DAC values [18,31,129]
%Calculate XYZ from the DAC-values using LUT
test_XYZ_LUT=DigittoXYZ(r,g,b,offset,test(6,:),cmf)';
disp('test_XYZ_LUT model value');%  observe what are the XYZ values calculated from the LUT and GOG model
disp(abs(test_XYZ_LUT));
%Calculate XYZ from the DAC-values using GOG model
test_XYZ_GOG=DigittoXYZ_GOG(r,g,b,offset,test(6,:),cmf)';
disp('test_XYZ_GOG_model value');
disp(abs(test_XYZ_GOG));
%Calculate LAB values from the LUT
Real_Predicted_LabDifference_LUT =zeros(1,6);
for i = 1:6
    lab=XYZtoLab(DigittoXYZ(r,g,b,offset,test(i,:),cmf)', CalculateXYZ(screen_name{5}(1,:),cmf));   % XYZ values from the DAC values
    
    real1=XYZtoLab(CalculateXYZ(screen_name{5}(i,:),cmf), CalculateXYZ(screen_name{5}(1,:),cmf));
    
    Real_Predicted_LabDifference_LUT(i)=dEcalc(lab,real1);
end
Real_Predicted_LabDifference_GOG = zeros(1,6);
%Calculate LAB values from the GOG
for i = 1:6
    lab=XYZtoLab(DigittoXYZ_GOG(r,g,b,offset,test(i,:),cmf)',CalculateXYZ(screen_name{5}(1,:),cmf)); % XYZ values from the DAC values
    
    real1=XYZtoLab(CalculateXYZ(screen_name{5}(i,:),cmf),CalculateXYZ(screen_name{5}(1,:),cmf));
    
    Real_Predicted_LabDifference_GOG(i)=dEcalc(lab,real1);
end

% Sum all Errors
AveragedError_LUT=sum(Real_Predicted_LabDifference_LUT)/6;
disp('Absolute value of AverageError_LookUptable');
disp(abs(AveragedError_LUT));
AveragedError_GOG=sum(Real_Predicted_LabDifference_GOG)/6;
disp('Absolute value of AverageError_GoG model');
disp(abs(AveragedError_GOG));

%**********************
% COMMENTS 
%**********************
% The absolute value of average error of the GoG model was relatively
% better than the LUT table for the selected DAC values. This depicts that
% the CRT is quite better of with GOG model which helps in optimizing the linear transfromation matrix for better accuracy in the uniform CIELAB sapce. 
%With these characterization errors there exist also quantization errors.
% Dell24 and Dell Konica display have a higher characterization error due 
% quantization error. Additionally, there may also be temporal instability
% of the display which may incur other errors. 
% In addition, the GOG model for Dell Konica was quite better compared with the LUT
% interpolation. Error of about 15.14 difference is color was quite not
% really good but was better than LUT interpolation. 
% Moreover, the LUT interpolation and GOG model showed appreciably almost
% the same error results for the the tone curve characterization. 

%There is a significant difference between the calibration performance of 
%the Dell monitor and Eizo monitor. This difference is due to the short term 
%instability of the Eizo monitor, i.e. the maximum luminance of the phosphors change within a 
%spam of a time causind differences in tristimulus values.

% The EIZO display showed relative good fitting for both LUT and GOG models
% which depicts a small difference between the DAC values and the
% tristimulus values. 
% The projector showed a very inaccurate results in the GOG model which
% depicts that this fitting model is relatively not good for projectors. 
