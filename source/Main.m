clc;
clear all;
close all;
warning off;
[file,path]=uigetfile('.jpg','Pick An Image');
im=imread([path file]);
figure;
imshow(im);
title('Input Image');
[R,C,col]=size(im);
if col==3
    In_Img=rgb2gray(im);
    figure;
imshow(im);
title('Gray Image');
else
    In_Img=im;
end
In_Img=imresize(In_Img,[1210 2756]);
[H,W]=size(In_Img);



In_Img=imcrop(In_Img,[1757.45 780.745652173914 721.952173913043 257.626086956522]);
figure;
imshow(In_Img);
title('ROI Extracted Image');
% %Filter-Preprocessing

In_fil=medfilt2(double(In_Img));
figure;imshow(uint8(In_fil));title('PreprocessedImage');

%FeatureExtraction
% originalImage=im;
corners=detectHarrisFeatures(In_Img);

figure;
imshow(In_Img);hold on;
plot(corners.selectStrongest(1000));
title('InputFeaturesImage');
I_thresh=In_Img>80 & In_Img<170;


cc=bwconncomp(I_thresh,26);
props=regionprops(I_thresh,'Area');
L=labelmatrix(cc);
I_thresh=ismember(L,find([props.Area]<=10000 & [props.Area]>=100));
figure;
imshow(I_thresh);title('ThresholdSegmentation');
points=detectSURFFeatures(In_Img);
figure;
imshow(I_thresh);hold on;
plot(points.selectStrongest(25));
title('Detected BRISK Points');
[feat,vp]=extractFeatures(In_Img,points);
testfea=feat(1,:);
load trainfea 
load label;
[result] = multisvm(trainfea,label',testfea);

if result==0
    msgbox('Signature is Original');
else
    msgbox('Signature is Forged');
end

