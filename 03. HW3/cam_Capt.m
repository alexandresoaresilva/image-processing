clc, clear, close all
% webcam setup
% webcamlist_conn = webcamlist;
% webcam_connected = webcam(webcamlist_conn{1});
% webcam_connected.Resolution ='640x480';
% 
% webcam_connected.Resolution = '320x240';
% prv = preview(webcam_connected);
% prv.Parent.Visible = 'off';
% prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
% prv.Parent.Visible = 'on';
% prompt = 'Press enter to capture a frame';
% x = input(prompt);
% I = snapshot(webcam_connected);
% delete(prv);

%% load file
% I = imread('fuzzy1_2.jpg');
I = imread('small_Digits.jpg');

if length(size(I)) > 2
    I = rgb2gray(I);
end

%% procesing
% I = imresize(I,[240 320]);
[M, N] = size(I);
LPF = avg_filt(2);
% I1 = I;
I1 = conv2(I,LPF,'valid');

I_bin = ~imbinarize(uint8(I1),.4);
for i=1:2
    I_bin = medfilt2(I_bin);
end
%imshow()
%I1 = preprocess(I)

% I1 = double(I)-I1;
% Iw1 = imgaussfilt(I,3);
% Iw2 = wiener2(I1,[5 5]);
% Iw2 = single(edge(I2,'Canny'));
%% matlab nonuniform illumination
% [I2, I3, I4] = matlab_nonuni_illu_remov(I,7, 15);
se = strel('disk',6);
% I_bin = conv2(double(I_bin)*255,LPF,'valid');
% I_bin = uint8(I_bin);
% I_bin = imbinarize(I_bin,.5);
I_binclosed = imclose(I_bin,se);
I_cell = {I, I1, I_bin, I_binclosed};
I_props = regionprops(I_binclosed);

figure
imshow(I_binclosed);
hold on;
digits = cellmat(0);
j = 1;
for i=1:length(I_props)
    if I_props(i).Area > (M*N / 330)
        if I_props(i).BoundingBox(4) > N/8
            rect = rectangle('Position',I_props(i).BoundingBox,'EdgeColor','r','LineWidth',3);
            
            digit = imcrop(I_binclosed, I_props(i).BoundingBox);
%             digits{i} = imcrop(I_binclosed, I_props(i).BoundingBox);
            [m,n] = size(digit);
            [max_x,i] = max([m,n]);            
%           0  100
            
            digits{j} = pad_digit_img(digit ,ceil(max_x*1.2));
            j = j + 1;
        end
    end
end

figure
sub_plots(digits);


%% functions
function I = pad_digit_img(I,square_side)
    [x,y] = size(I);
    [x_pre, x_pos] = padding_pre_post_calc(x,square_side);
    [y_pre, y_pos] = padding_pre_post_calc(y,square_side);

    I = padarray(I,[x_pre y_pre],'pre');
    I = padarray(I,[x_pos y_pos],'pos');


    function [n_pre, n_pos] = padding_pre_post_calc(n,dim_size)
        
        n_is_not_even = rem(dim_size-n,2);
        n_pre = (dim_size-n)/2;
%         n_pre_is_not_even = rem(n_pre,2);
        
        if n_is_not_even && n_pre %n_pre not zero
            n_pre = n_pre - 0.5;
            n_pos = n_pre + 1;
        else
            n_pos = n_pre;
        end
        n_pos = double(n_pos);
        n_pre = double(n_pre);
    end
end

function sub_plots(I_cell)
    n = length(I_cell);
    no_of_plots_per_row = ceil(n/2);
    for i=1:n
        subplot(2,no_of_plots_per_row,i)
        imshow(I_cell{i},[]);
    end
end

function [I2, I3, I4] = matlab_nonuni_illu_remov(I,disk_size, bin_thresh)
    se = strel('disk',disk_size);
    background = imopen(I,se);
%     LPF = avg_filt(50);
%     background = conv2(I,LPF,'same');
    %imshow(background);
    I2 = I - background;
    %imshow(I2,[]);
    I3 = imadjust(I2);
    %imshow(I3,[]);
    bw = imbinarize(I3);
    I4 = bwareaopen(bw,bin_thresh);
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end