%% Clear Workspace
clc; clear;
close all;
%% Read Image
im = imread('Proj2.tif');
subplot(2,3,1);
imagesc(im);
title('Orignal Image');
colormap gray
axis image;
%% Correct Illumination
correct_illi_im = correctIlli(im);
subplot(2,3,2);
imagesc(correct_illi_im);
title('Illumination Corrected Image');
colormap gray
axis image;

%% Extract Periodic Pattern
[orig_fft, fft_im, filt_fft_im, periodic_pattern_normal] = extract_mesh(correct_illi_im, 11.3);

subplot(2,3,3);
imagesc(fft_im);
title('FFT');
colormap gray
axis image;

subplot(2,3,4);
imagesc(filt_fft_im);
title('Low Pass Filtered FFT')
colormap gray
axis image;

subplot(2,3,5);
imagesc(periodic_pattern_normal);
title('Extracted Periodic Pattern')
colormap gray
axis image;

% %% Matlab
% background = imopen(im,strel('disk',15));
% mat_corrected_illi = im-background;
% subplot(2,3,6);
% imagesc(mat_corrected_illi);
% title('Matlab Illi Corrected')
% colormap gray
% axis image;
% 
% [mat_orig_fft, mat_fft_im, mat_filt_fft_im, mat_periodic_pattern_normal] = extract_mesh(mat_corrected_illi, 11.3);
% 
% 
% %% Freq Image Illi Corr 2
% fft_i2 = fft_im;
% fft_i2 = (fft_i2+2*ones(size(fft_i2)));
% fft_i2(fft_i2<0) = 0;
% 
% figure
% subplot(3,1,1)
% imagesc(fft_im)
% colormap gray
% 
% subplot(3,1,2)
% imagesc(mat_fft_im)
% colormap gray
% [filt_g, blur_im] = gaussian_blur(im,1.42);
% subplot(3,1,3)
% imagesc(log(1+abs(fftshift(fft(blur_im)))));
% colormap gray
% correct_illi_im = im-blur_im;
% 
% % Using Freq Domain
% 
% 
% % 
% ift = abs(ifft2(orig_fft));
% figure
% mesh = uint8(ift);
% imagesc(mesh);
% colormap gray
% axis image;

%% Functions
function correct_illi_im = correctIlli(im)
    [filt_g, blur_im] = gaussian_blur(im,1.4);
    correct_illi_im = im-blur_im;
end

function [filt_g,blur_im] = gaussian_blur(im, std)
    guass_filt =@(x,y,std)exp(-(x.^2+y^2)/(2*std^2))*1/(2*pi*std^2);
    x = -4:1:4;
    y= x;
    len = length(x);
    for i=1:len
        for j=1:len
            filt_g(i,j) = guass_filt(x(i),y(j),4);
        end
    end
%     [X,Y] = meshgrid(linspace(-4,4,len),linspace(-4,4,len));
%     surf(X,Y,filt_g);
    blur_im = uint8(conv2(im,filt_g,'same'));
end


function [orig_fft, fft_im, filt_fft_im, mesh] = extract_mesh(im, threshold)

    %% Take 2 Dimensional Fourier Transform
    orig_fft = fft2(im);
    shifted_and_scaled_fft = log(1+fftshift(orig_fft));
    X = abs(shifted_and_scaled_fft);
    fft_im = X;
    

    %% Filter for Testing
    tempX = X;
    tempX(tempX<threshold) = 0;
    filt_fft_im = tempX;

    %% Apply to orignal Image
    orig_fft(orig_fft<exp(threshold)) = 0.0;

    %% Inverse Fourier Transform
    ift = abs(ifft2(orig_fft));
    mesh = uint8(ift);
    
  
end

function blur_f = create_blur_filter(n)
    blur_f = (1/n^2)*ones(n,n);
end

function lpf_filt = lpf_filter(im, cut_off_radius)
    lpf_filt = ones(size(im));
    [u,v] = size(im);
    center_u = u-u/2;
    center_v = v-v/2;
    
    
    [X,Y]=meshgrid(1:size(im,2),1:size(im,1));
    disk_locations=sqrt((X-center_v).^2+(Y-center_u).^2) <= cut_off_radius;
    lpf_filt = lpf_filt.*disk_locations;
    
    figure
    imagesc(lpf_filt);
    colormap gray
    title('Low Pass Filter');
end

function hpf_filt = hpf_filter(im, cut_off_radius, value)
    hpf_filt = ones(size(im));
    [u,v] = size(im);
    center_u = u-u/2;
    center_v = v-v/2;
    
    
    [X,Y]=meshgrid(1:size(im,2),1:size(im,1));
    disk_locations=sqrt((X-center_v).^2+(Y-center_u).^2) <= cut_off_radius;
    disk_locations = ~disk_locations;
    hpf_filt = hpf_filt.*disk_locations;
    hpf_filt(hpf_filt == 0) = value;
    
    figure
    imagesc(hpf_filt);
    colormap gray
    title('High Pass Filter');
end
