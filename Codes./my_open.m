% clear all
% close all
% 
% %fileName = 'D:\Matlab ENVI\read_ENVI\2018-09-20_004\capture\2018-09-20_004';
% fileName = 'D:\Temp\t1\2019-01-16_005\capture\2019-01-16_005';
% 
% 
% 
% [imageStack,wavelengths] = read_ENVI(fileName);
% 
% 
% %imshow(imageStack(:,:,100)/1000)
% 
% imagesc(imageStack(:,:,100)/1000)
% plot(squeeze(imageStack(200,200,:)/1))
% %%


clear all
close all

fileName = 'D:\Matlab ENVI\read_ENVI\Cartilage 1 2018-11-06_003\capture\2018-11-06_003';

[imageStack,wavelengths] = read_ENVI(fileName);


imshow(imageStack(:,:,100)/1000)

plot(squeeze(imageStack(200,200,:)/1))