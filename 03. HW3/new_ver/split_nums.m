% split
clear, clc, close all;

SIDE = 240;
addpath('numbers');
% addpath('VLFEATROOT');
I = imread('0-5.JPG');
I2 = imread('6-9.JPG');
[d, d_bin] = extract_digits(I, 2, SIDE);
[d2, d2_bin] = extract_digits(I2, 2, SIDE);

L = length(d_bin);
L2 = length(d2);
digits_bin = cell(1,L +L2);
digits = cell(1,L +L2);

j = 1;
for i=1:L
     digits_bin{i} = d_bin{i};
     digits{i} = d{i};
%     d{i} = d2{j};
    j = j + 1;
end

j = 1;
L = L + 1;
L2 = length(d2) + L-1;

for i=L:L2
    digits_bin{i} = d2_bin{j};
    digits{i} = d2{j};
    j = j + 1;
end



index = [1 6 9 3 4 11 7 10 2 5 12 8 13 15 17 19 20 16 18 21];
final_digits = digits(index);
final_digits_bin = digits_bin(index);

sub_plots(final_digits);
sub_plots(final_digits_bin);

%11, 12 
keySet = [0 11 12 21 22 23 31 32 41 42 43 51 52 61 62 71 72 81 82 9];
numbers = ["zero", "one", "one", "two", "two", "two", "three", "three",...
    "four", "four", "four", "five", "five", "six", "six", "seven",...
     "seven", "eight", "eight", "nine"];

digits = containers.Map(keySet,final_digits);
number_imgs{1} = digits;
number_imgs{2} = keySet;
number_imgs{3} = numbers;

digits_bin = containers.Map(keySet,final_digits_bin);
number_imgs_bin{1} = digits_bin;
number_imgs_bin{2} = keySet;
number_imgs_bin{3} = numbers;

save('number_imgs.mat','number_imgs');
save('number_imgs_bin.mat','number_imgs_bin');