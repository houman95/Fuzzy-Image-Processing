clc;
clear;
close all;

%% Load Image Data

Filter={'*.jpg;*.jpeg;*.png'};

[FileName, FilePath]=uigetfile(Filter);

if FileName==0
    return;
end

FullFileName=[FilePath FileName];

img=imread(FullFileName);
img = rgb2gray(img);
%img=im2double(img);
img2 = NAFSM(img);
subplot(1,2,1);
imshow(img);
title('Original Image');

subplot(1,2,2);
imshow(img2);