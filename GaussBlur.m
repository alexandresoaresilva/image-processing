clc; clear;
guass_filt =@(x,y,std)exp(-(x.^2+y^2)/(2*std^2))*1/(2*pi*std^2);

x = -4:1:4;
y= x;
len = length(x);

for i=1:len
    for j=1:len
        filt_g(i,j) = guass_filt(x(i),y(j),4);
    end
end

[X,Y] = meshgrid(linspace(-4,4,len),linspace(-4,4,len));
surf(X,Y,filt_g);
%Z = sin(X) + cos(Y);

a = imread('TestImage1.tif');
imshow(a);
figure
blur_a = imshow(uint8(conv2(a,filt_g,'same')));

