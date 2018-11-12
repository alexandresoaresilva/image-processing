clc, clear, close all

I = rgb2gray(imread('6.jpg'));

x = fspecial('motion', 25, 0) ;
y = fspecial('motion', 25, 90) ;

for i = 1:1
    I2 = conv2(I,x,'same');
%     I2 = conv2(I,y,'same');
end

figure(1)
imshow(I2,[])
