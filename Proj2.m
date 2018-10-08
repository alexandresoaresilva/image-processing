%% Clear Workspace
clc; clear;

%% Read Image
im = imread('Proj2.tif');
% imagesc(im);
% title('Orignal Image');
% colormap gray
% axis image;
%% Correct Illumination
correct_illi_im = correctIlli(im);

%% Extract Periodic Pattern
periodic_pattern_normal = extract_mesh(correct_illi_im, 11.3);

%% Correct Illumination 2
% fft_image = fft2(im);
% X = log(1+abs(fftshift(fft_image)));
% figure;
% imagesc(X);
% title('FFT');
% colormap gray
% axis image;
% 
% tempX = X(180:230,265:280);
% tempX(tempX>12) = 7.5917;
% tempX = X(180:230,265:280);
% tempX(:,:) = 0;
% X(180:230,265:280) = 0;
% tempX(tempX>15) = 0;
% figure
% imagesc(X);
% title('High Pass Filtered FFT')
% colormap gray
% axis image;
% 
% figure
% imagesc(tempX);
% title('High Pass Filtered FFT')
% colormap gray
% axis image;


% fft_image(fft_image>exp(12)) = 0.0;
% 
% ift = abs(ifft2(exp(ifftshift(tempX))-1));
% 
% figure
% mesh = uint8(ift);
% imagesc(mesh);
% colormap gray
% axis image;

%% Functions
function correct_illi_im = correctIlli(im)
    [filt_g, blur_im] = gaussian_blur(im,1.42);
    correct_illi_im = im-blur_im;
    figure
    imagesc(correct_illi_im);
    title('Illumination Corrected Image');
    colormap gray
    axis image;
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


function mesh = extract_mesh(im, threshold)

    %% Take 2 Dimensional Fourier Transform
    fft_image = fft2(im);
    X = log(1+abs(fftshift(fft_image)));
    figure;
    imagesc(X);
    title('FFT');
    colormap gray
    axis image;

    %% Filter for Testing
    tempX = X;
    tempX(tempX<threshold) = 0;
    figure
    imagesc(tempX);
    title('Low Pass Filtered FFT')
    colormap gray
    axis image;

    %% Apply to orignal Image
    fft_image(fft_image<exp(threshold)) = 0.0;

    %% Inverse Fourier Transform
    ift = abs(ifft2(fft_image));
    
    figure
    mesh = uint8(ift);
    imagesc(mesh);
    colormap gray
    axis image;
  
end

function blur_f = create_blur_filter(n)
    blur_f = (1/n^2)*ones(n,n);
end
