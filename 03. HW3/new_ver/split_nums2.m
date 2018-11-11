clear, clc, close all;
addpath('numbers');

% for 0s
a = load('image_Templates.mat');
b = load('number_imgs.mat');
c = load('image_Templates2');

d0 = cell.empty;
d1 = cell.empty;
d2 = cell.empty;
d3 = cell.empty;
d4 = cell.empty;
d5 = cell.empty;
d6 = cell.empty;
d7 = cell.empty;
d8 = cell.empty;
d9 = cell.empty;
%% ZERO
d0{1} = a.image_Templates{1};
d0{2} = b.number_imgs{1}(0);
d0{3} = c.image_Templates{1};

%% ONE
d1{1} = a.image_Templates{2};
d1{2} = b.number_imgs{1}(11);
d1{3} = b.number_imgs{1}(12);
d1{4} = c.image_Templates{2};
%% TWO
d2{1} = a.image_Templates{3};
d2{2} = c.image_Templates{3};
d2{3} = b.number_imgs{1}(21);
d2{4} = b.number_imgs{1}(22);
%% THREE
d3{1} = a.image_Templates{4};
d3{2} = c.image_Templates{4};
d3{3} = b.number_imgs{1}(31);
d3{4} = b.number_imgs{1}(32);
%% FOUR
d4{1} = a.image_Templates{5};
d4{2} = c.image_Templates{5};
d4{3} = b.number_imgs{1}(41);
d4{4} = b.number_imgs{1}(42);
%% FIVE
d5{1} = a.image_Templates{6};
d5{2} = c.image_Templates{6};
d5{3} = b.number_imgs{1}(51);
d5{4} = b.number_imgs{1}(52);
%% SIX
d6{1} = a.image_Templates{7};
d6{2} = c.image_Templates{7};
d6{3} = b.number_imgs{1}(61);
d6{4} = b.number_imgs{1}(62);
%% SEVEN
d7{1} = a.image_Templates{8};
d7{2} = c.image_Templates{8};
d7{3} = b.number_imgs{1}(71);
d7{4} = b.number_imgs{1}(72);
%% EIGHT
d8{1} = a.image_Templates{9};
d8{2} = c.image_Templates{9};
d8{3} = b.number_imgs{1}(81);
d8{4} = b.number_imgs{1}(82);
%% NINE
d9{1} = a.image_Templates{10};
d9{2} = c.image_Templates{10};
d9{3} = b.number_imgs{1}(9);

%% 
% digits = cell.empty;
digit_imgs = {d0, d1, d2, d3, d4, d5, d6, d7, d8, d9};
digits_strg = ["zero", "one", "two", "three", "four", "four", "five", "six",...
    "seven", "eight", "nine"];
keySet = [0 1 2 3 4 5 6 7 8 9];

digit = containers.Map(keySet,digit_imgs);
plot_i = 251;
for i=1:length(c.image_Templates)
    subplot(plot_i);
    imshow(c.image_Templates{i});
    plot_i = plot_i + 1;
end