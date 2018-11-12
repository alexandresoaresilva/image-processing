% split
clear, clc, close all;


% old = cd('E:\Google Drive\02. Image Processing\hw3_database');
% [file_names2, file_names2_char] = save_file_names_in_folder(pwd,'jpg');
% load('E:\_local Home\Documents\Git\image-processing\03. HW3\new_ver\numbers\digit.mat')
% number = cell(10,1);
% numberbin = cell(10,1);


keySet = [];
for i=1:size(file_names2,1)
    SIDE = 240;
    I = imread(file_names2(i));
    [d, dbin] = extract_digits(I, 4, SIDE);
    keySet(i) = i - 1;
    number{i} = d;
    numberbin{i} = dbin;
    L(i) = length(d);
    Lbin(i) = length(dbin);
end

min_L = min(L);
% dig = cell(1,min_L);
for i=1:length(number)
    d = digit(i-1);
    
    number{i} = {number{i}{1:min_L}, d{1:end}};
    L(i) = length(number{i});
end

min_L = min(L);
% dig = cell(1,min_L);
for i=1:length(number)
    number{i} = number{i}(1:min_L);
end



% min_LBin = min(Lbin);
% % digbin = cell(1,min_Lbin);
% for i=1:length(number)
%     numberbin{i} ={numberbin{i}(1:min_LBin), ;
% end


digits_map = containers.Map(keySet,number);
% digits_bin_map = containers.Map(keySet,numberbin);

cd(old);


for i=0:9
    di = digits_map(i);
    figure
    for j=1:length(di)
        subplot(2,5,j);
        imshow(di{j});
    end
end
% 
% for i=0:9
%     di = digits_bin_map(i);
%     figure
%     for j=1:length(di)
%         subplot(2,5,j);
%         imshow(di{j});
%     end
% end

function [file_names2, file_names2_char] = save_file_names_in_folder(img_folder,extension)
    %gets file names with the selected extension
    current_folder = pwd; %saving so the program can return to the original  folder

    cd(img_folder);
    if extension(1) ~= '*'
        if extension(1) ~= '.'
            extension = char(strcat('*.',extension));
        else
            extension = char(strcat('*',extension));
        end
    end

    file_names = struct2cell(dir(extension));
    file_names2 = string.empty(0, length(file_names(1,:)) );

    for i=1:size(file_names,2)%no. of columns
        %file_name_dummy = cell2mat(file_names(1,i));
        file_name_dummy = file_names{1,i}(1,:);
        file_name_dummy = string(file_name_dummy);
        if i == 1
            file_names2 = file_name_dummy;
        else
            file_names2 = [file_names2; file_name_dummy];
        end
    end
    file_names2_char = char(file_names2);
    cd(current_folder);
end