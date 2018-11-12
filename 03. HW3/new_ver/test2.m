% addpath('numbers');
% load('digits_map.mat');
% 
% I = cell.empty;
% for i=1:3
%    [img,~] = extract_digits(cam_capt, 3, 240);
%    I{i} = img{1};
% %    pause
% end
% number = cell(1,11);
% 
% for i=0:9
%     number{i+1} = digits_map(i);
% end
% number{11} = I;
% 
% keySet = 0:10;
% digits_map = containers.Map(keySet,number);

load('digits_map.mat');
for i=0:10
    d = digits_map(i);
    figure;
    for j=1:length(d)
        subplot(2,4,j)
        imshow(d{j});
    end
   
end