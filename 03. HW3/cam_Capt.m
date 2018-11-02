clc, clear, close all

filter=@(n)1/(n^2)*ones(n);

%% webcam setup
% webcamlist_conn = webcamlist;
% webcam_connected = webcam(webcamlist_conn{1});
% webcam_connected.Resolution ='640x480';
% preview(webcam_connected);
% prompt = 'Press enter to capture a frame';
% x = input(prompt);
% I = snapshot(webcam_connected);
% I = rgb2gray(I);


I = imread('7.jpg');
LPF = filter(13);

%I1 = preprocess(I)
I1 = conv2(I,LPF,'same');
%I1 = imgaussfilt(I,10);
% I2 = wiener2(I1,[5 5]);
% I2 = single(edge(I2,'Canny'));
%% matlab nonuniform illumination
[I2, I3, I4] = matlab_nonuni_illu_remov(I,6, 15);
se = strel('disk',20);

I5 = imclose(I4,se);

I_cell = {I1, I2, I3, I4, I5};

figure
sub_plots(I_cell, length(I_cell));

colormap gray;

function sub_plots(I_cell,n)
    for i=1:n
        subplot(2,3,i)
        imshow(I_cell{i},[]);
    end
end

function [I2, I3, I4] = matlab_nonuni_illu_remov(I,disk_size, bin_thresh)
    se = strel('disk',disk_size);
    background = imopen(I,se);
    %imshow(background);
    I2 = I - background;
    %imshow(I2,[]);
    I3 = imadjust(I2);
    %imshow(I3,[]);
    bw = imbinarize(I3);
    I4 = bwareaopen(bw,bin_thresh);
end
function p = preprocess(I)


    


    

    HPF = 1 - LPF;
    %BPF = 1 - (LPF + HPF);       %Band-pass Filter
    
    HP_content = conv2(I,HPF,'same');
    
%     LP_content = LPF .* fftI;         %relevant frequency regions
%     HP_content = HPF .* fftI;
    
    LP_content = fftshift(fft2(LP_content));
    HP_content
    [M, N] = size(I);

    %fftI=fftshift(fft2(I));

    %LPF = mat2gray(fspecial('gaussian', [M N], 15));    %Low-pass Filters
    %LPF = filter2(fspecial('average',15),I)./255;
     
    %LPF=LPF~=0;
%     LPF2 = mat2gray(fspecial('gaussian', [M N], 5));
%     LPF2=LPF2~=0;
%     HPF = 1 - LPF2;                  %High-pass Filter
%     BPF = 1 - (LPF + HPF);           %Band-pass Filter
%     %LP_content = LPF .* fftI;        %Partition spectrum into
%     BP_content = BPF .* fftI;         %relevant frequency regions
%     HP_content = HPF .* fftI;
%     HP_content = HP_content .* 0.01;       %Reduce noise by dampening
    
    backgrnd = uint8(abs(ifft2(LP_content)));
    backgrnd(:) = mean(backgrnd(:));
    LP_content = fftshift(fft2(backgrnd));
    
    %high frequency content
    newSpectrum = LP_content + HP_content;
    p = uint8(abs(ifft2(newSpectrum)));
end