% clc;
% clear all;
% close all;
ij=1;
for ii=19:36
im=imread(['E:\Signature Verification Code\forged signature\',num2str(ij),'.jpg']);
ij=ij+1;
[R,C,col]=size(im);
if col==3
    In_Img=rgb2gray(im);
    
else
    In_Img=im;
end
In_Img=imresize(In_Img,[1210 2756]);
[H,W]=size(In_Img);



In_Img=imcrop(In_Img,[1757.45 780.745652173914 721.952173913043 257.626086956522]);
 
% %Filter-Preprocessing

In_fil=medfilt2(double(In_Img));

%FeatureExtraction
% originalImage=im;
corners=detectHarrisFeatures(In_Img);


I_thresh=In_Img>80 & In_Img<170;


cc=bwconncomp(I_thresh,26);
props=regionprops(I_thresh,'Area');
L=labelmatrix(cc);
I_thresh=ismember(L,find([props.Area]<=10000 & [props.Area]>=100));
points=detectSURFFeatures(In_Img);
[feat,vp]=extractFeatures(In_Img,points);
trainfea(ii,:)=feat(1,:);
end
save trainfea trainfea