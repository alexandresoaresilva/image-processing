clc; clear;
im = imread('Proj2.tif');
figure
imagesc(im);
colormap gray
title('Orignal Image');


%% Matlab Correct Illumination
% background = imopen(im,strel('disk',15));
% mat_corrected_illi = im-background;
% figure
% imagesc(mat_corrected_illi);
% colormap gray
% title('Matlab Image');

%% Extract Pattern
[pattern, correct_ilu] = (extractPattern(im));
pattern = uint8(pattern);
figure
imagesc(pattern);
colormap gray
title('Periodic Pattern in Image');

%% Correct Illumination
% im2 = im+50;
% for j=1:size(im2,2)
%     im2(:,j) = (im2(:,j)-mean(im2(:,j)))/(std(correct_ilu(:,j)));
% end
% for j=1:size(im2,1)
%     im2(j,:) = (im2(j,:)-mean(im2(j,:)));
% end
% figure
% imagesc(im2);
% colormap gray
% title('Correct Illumination Profile Depth');

% correct_illi_im = correctIlli(im);
% figure
% imagesc(correct_illi_im);
% colormap gray
% title('Correct Illumination of Pattern');

% figure
% imagesc(correct_ilu);
% colormap gray
% title('Correct Illumination of Image');


%% Functions
function [filt_img, filt_img2] = extractPattern(im)
    

    orig_fft = fftshift(fft2(im));
    figure
    subplot(2,2,1)
    imagesc(abs(orig_fft));
    title('Original FFT');
    colormap gray
    axis image;

    shifted_and_scaled_fft = log(1+(orig_fft));
    X = abs(shifted_and_scaled_fft);
    subplot(2,2,2)
    imagesc(X);
    title('Shifted and Scaled FFT');
    colormap gray
    axis image;

    amplitudeThreshold = 11.9;
    nobrightSpikes = X < amplitudeThreshold; 
    subplot(2,2,3)
    imagesc(nobrightSpikes);
    title('High Pass Filter');
    colormap gray
    axis image;

    orig_fft_high_filter = orig_fft;
    orig_fft_high_filter(nobrightSpikes) = 0;
    scaled_orig_high_fft = log(abs(orig_fft_high_filter));
    minValue = min(min(scaled_orig_high_fft));
    maxValue = max(max(scaled_orig_high_fft));
    subplot(2,2,4)
    imagesc(scaled_orig_high_fft, [minValue maxValue]);
    title('High Pass Filter Masked Image');
    colormap gray
    axis image;
    
    amplitudeThreshold2 = 11.8;
    brightSpikes = (X > amplitudeThreshold2) ; 
    figure
    subplot(1,2,1)
    imagesc(brightSpikes);
    title('Low Pass Filter');
    colormap gray
    axis image;
    
    orig_fft_low_filter = orig_fft;
    orig_fft_low_filter(brightSpikes) = 0;
    orig_fft_low_filter(1:175,273) = orig_fft_low_filter(1:175,272);
    orig_fft_low_filter(235:size(orig_fft,1),273) = orig_fft_low_filter(235:size(orig_fft,1),272);
    scaled_orig_low_fft = log(abs(orig_fft_low_filter));
    minValue = min(min(scaled_orig_low_fft));
    maxValue = max(max(scaled_orig_low_fft));
%     figure
    subplot(1,2,2)
    imagesc(scaled_orig_low_fft, [minValue maxValue]);
    title('High Pass Filter Masked Image');
    colormap gray
    axis image;

    filteredImage = ifft2(fftshift(orig_fft_high_filter));
    filt_img = abs(filteredImage);
    
    filteredImage2 = ifft2(fftshift(orig_fft_low_filter));
    filt_img2 = abs(filteredImage2);
end

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


