close all; clc, clear

filt_n = @(n)1/n^2*ones(n);
adjust_lg_abs_I_fft =...
    @(lg_abs_I_fft)uint8(round(lg_abs_I_fft./max(max(lg_abs_I_fft))*255));
I_read = imread('Proj2.tif');
I = I_read;
%I_conv = conv2(I,filt_n(9),'same');
%I_conv = I_read - uint8(I_conv);
%I = correct_negatives(I_conv)

I_fft = fft2(I);
abs_I_fft = abs(fftshift(I_fft));
lg_abs_I_fft = log(abs_I_fft);
figure
%lg_abs_I_fft = adjust_lg_abs_I_fft(lg_abs_I_fft);
%I_conv = conv2(lg_abs_I_fft,filt_n(3),'same');
%I_conv_th = I_conv;
% I_conv_th =  binI(163, I_conv);
% I_conv_th =  binI(170, lg_abs_I_fft);

pattern_i = find(lg_abs_I_fft > 11 );

background_i_lo = find(lg_abs_I_fft <= 11.5);
% background_i_hi = find(lg_abs_I_fft > 13.4);
% (background_i_lo == background_i_hi)
% pattern_i = find(lg_abs_I_fft > 168 );
% background_i = find(lg_abs_I_fft <= 168);

% [m,n] = size(I_conv_th);
% m = m/2;
% n = n/2;
I_g = fftshift(I_fft); 
I_g(background_i_lo ) = 0;
% I_g(background_i_hi) = 0;

%I_ifft = ifft2(fftshift(I_g));
I_ifft = ifft2(fftshift(I_g));
I_ifft = abs(I_ifft);
% a = uint8(round(I_ifft/max(max(I_ifft))*255));
a = uint8(round(I_ifft));
[m,n] = size(lg_abs_I_fft);

a = histeq(a);
%imshow([I, lg_abs_I_fft, lg_abs_I_fft; lg_abs_I_fft, I_g, a]);

imagesc([I, lg_abs_I_fft, lg_abs_I_fft; lg_abs_I_fft, I_g, a]);
title('original, lg_abs_I_fft, I_conv');
colormap gray

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
    filt_n =...
        @(x,mean,std_dev)1./sqrt(2*pi.*std_dev^2).*exp(-(x-mean)^2./(2*std_dev^2));
    
    for i=-n:n
        filt(
    end
end