%--------------------------------------------------------------------------
%SIFT Based Digit-RecognitionS
%Date: 10/30/18
%--------------------------------------------------------------------------

clc
close all
clear

run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose

I1 = single(rgb2gray(imread('6.jpg')));
I2 = single(rgb2gray(imread('6','jpg')));
I3 = single(rgb2gray(imread('9','jpg')));

[M, N] = size(I1);

I1 = imgaussfilt(I1,20);
I1 = wiener2(I1,[10 10]);
I1 = single(edge(I1,'Canny'));

se = strel('disk',5);
I1 = imclose(I1,se);

figure(1)
subplot(1,2,1)
imshow(I1,[]);

[f,d] = vl_sift(I1);

perm = randperm(size(f,2)) ;
P = length(perm);
sel = perm(1:P);
h1 = vl_plotframe(f(:,sel)) ;
h2 = vl_plotframe(f(:,sel)) ;
set(h1,'color','k','linewidth',2) ;
set(h2,'color','y','linewidth',1) ;
h3 = vl_plotsiftdescriptor(d(:,sel),f(:,sel)) ;
set(h3,'color','g') ;

I2 = imgaussfilt(I2,20);
I2 = wiener2(I2,[10 10]);
I2 = single(edge(I2,'Canny'));

I2 = imclose(I2,se);

subplot(1,2,2)
imshow(I2,[]);

[f1,d1] = vl_sift(I2);
perm = randperm(size(f1,2)) ;
P = length(perm);
sel = perm(1:P);
h1 = vl_plotframe(f1(:,sel)) ;
h2 = vl_plotframe(f1(:,sel)) ;
set(h1,'color','k','linewidth',2) ;
set(h2,'color','y','linewidth',1) ;
h3 = vl_plotsiftdescriptor(d1(:,sel),f1(:,sel)) ;
set(h3,'color','g') ;

[matches, scores] = vl_ubcmatch(d, d1) ;