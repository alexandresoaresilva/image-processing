close all; clc, clear

filt_n = @(n)1/n^2*ones(n);
adjust_lg_abs_I_fft =...
    @(lg_abs_I_fft)uint8(round(lg_abs_I_fft./max(max(lg_abs_I_fft))*255));
%% original
I = imread('Proj2.tif');
figure
subplot(231);
imshow(I);
title('original img');
%% fft
threshold = 12;
I = double(I) - conv2(I,norm_filt(3),'same');
I_fft = fft2(I);
abs_I_fft = abs(fftshift(I_fft));
lg_abs_I_fft = log(1+abs_I_fft);
background_i_lo = find(lg_abs_I_fft  <= threshold );
object_i_lo = find(lg_abs_I_fft  > threshold );
%lg_abs_I_fft(background_i_lo ) = 0;

x = -5:-1;
y = x;
lg_abs_I_fft(205,273 + x) = 0;
lg_abs_I_fft(205 +y,273) = 0;
lg_abs_I_fft(205+y,273 + x) = 0;
x = 1:5;
y = x;
lg_abs_I_fft(205,273 + x) = 0;
lg_abs_I_fft(205 +y,273) = 0;
lg_abs_I_fft(205+y,273 + x) = 0;

x = 1:5;
y = -5:-1;
lg_abs_I_fft(205+y,273 + x) = 0;
y = 1:5;
x = -5:-1;
lg_abs_I_fft(205+y,273 + x) = 0;

% SE = strel('diamond',1)
% I_er = imerode(lg_abs_I_fft ,SE);

subplot(232);
imagesc(lg_abs_I_fft);
title('fft log-abs');
%% ifft
I_fft(background_i_lo) = 0;
[m,n] = size(I_fft);
%m = round(m/2);
%n = round(n/2);

x = -5:-1;
y = x;
I_fft(205,273 + x) = 0;
I_fft(205 +y,273) = 0;
I_fft(205+y,273 + x) = 0;
x = 1:5;
y = x;
I_fft(205,273 + x) = 0;
I_fft(205 +y,273) = 0;
I_fft(205+y,273 + x) = 0;

x = 1:5;
y = -5:-1;
I_fft(205+y,273 + x) = 0;
y = 1:5;
x = -5:-1;
I_fft(205+y,273 + x) = 0;


I_g = I_fft;
I_ifft = fftshift(I_g);

%I_ifft(object_i_lo) = 
%I_ifft(background_i_lo) = 0;
I_ifft = abs(ifft2(fftshift(I_ifft)));
%I_ifft_8bit = I_ifft;

subplot(233);
colormap gray
imagesc(I_ifft);
title('ifft');




function S = correct_negatives(r) %corrects negative pixels
    [m, n] = size(r);
    r = double(r);
    r_reshaped = reshape(r, [1,m*n]);
    
    min_r = min(r_reshaped);
    max_r = max(r_reshaped);
    S=uint8(255.*(r- min_r)./(max_r-min_r));
end

function binImage =  binI(threshold, image)
    binImage = image;
    binImage(image<(threshold)) = 0;
    binImage(image>(threshold)) = 255;
end

function filt = norm_filt(n)
    guass_filt =@(x,y,std)exp(-(x.^2+y^2)/(2*std^2))*1/(2*pi*std^2);

    x = -n:1:n;
    y= x;
    len = length(x);
    for i=1:len
        for j=1:len
            filt(i,j) = guass_filt(x(i),y(j),4);
        end
    end
end