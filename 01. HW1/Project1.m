clc; clear;

%% Read the image
imageName = input('Enter the image name: ', 's');
image = imread(imageName);
imageOrg = image; % Copy of orignal Image
imro = image;
subplot(121)
imshow(image);

%% Binarize the Image
binaryThreshold = 150;
binImage =  binarizeImage(binaryThreshold, image);
imshow(binImage)
% figure
% imshow(imageOrg);

%% Extract the Corners
[row,col] = find(binImage == 255);
indices = [row,col];

mins = min(indices);
maxs = max(indices);

minY = mins(1);
minX = mins(2);
maxY = maxs(1);
maxX = maxs(2);

[cornerX, cornerY] = find(binImage(minY,:) == 255);
corner1 = [cornerY(1),minY];

[cornerX, cornerY] = find(binImage(maxY,:) == 255);
corner2 = [cornerY(1),maxY];

[cornerX, cornerY] = find(binImage(:,minX) == 255);
corner3 = [minX,cornerX(1)];

[cornerX, cornerY] = find(binImage(:,maxX) == 255);
corner4 = [maxX,cornerX(1)];

%hold on
% text(corner1(1),corner1(2), '1', 'FontSize', 30, 'Color', 'red');
% text(corner2(1),corner2(2), '2', 'FontSize', 30, 'Color', 'red');
% text(corner3(1),corner3(2), '3', 'FontSize', 30, 'Color', 'red');
% text(corner4(1),corner4(2), '4', 'FontSize', 30, 'Color', 'red');
% plot(corner1(1),corner1(2),'rx', 'MarkerSize', 20);
% plot(corner2(1),corner2(2),'gx', 'MarkerSize', 20);
% plot(corner3(1),corner3(2),'go', 'MarkerSize', 20);
% plot(corner4(1),corner4(2),'ro', 'MarkerSize', 20);

%% Check if image is straight
thresholdDistanceX = 100;
if(abs(corner4(1)-corner2(1))<thresholdDistanceX)
    
else
    %% Rotate Image
    % angle = atand(abs(corner4(2)-corner2(2))/(corner4(1)-corner2(1)));
    distanceThreshold = 220; % Threshold distance between two adjacent sides
    if(abs(corner4(1)-corner2(1)) < distanceThreshold)
        angle2 = rad2deg(atan2((corner3(2)-corner2(2)),(corner3(1)-corner2(1))));
%         imro = imrotate(imageOrg,90+angle2,'crop');
    else
        angle2 = rad2deg(atan2((corner4(2)-corner2(2)),(corner4(1)-corner2(1))));
        
    end
    imro = imrotate(imageOrg,90+angle2,'crop');
end


%% Binarize rotated image
binImage2 =  binarizeImage(binaryThreshold, imro);
% figure
% imshow(binImage2);

%% Crop Image
blackRowsXtop = 0;
blackRowsXbottom = 0;
blackRowsYleft = 0;
blackRowsYright = 0;

for i=1:size(binImage2,1)
    if sum(binImage2(i,:) == 255) > 0
        break
    else
        blackRowsXtop = blackRowsXtop+1;
    end
end
for i=blackRowsXtop:size(binImage2,1)
    if sum(binImage2(i,:) == 255) == 0
        blackRowsXbottom = blackRowsXbottom+1;
    end
end

for i=1:size(binImage2,2)
    if sum(binImage2(:,i) == 255) > 0
        break
    else
        blackRowsYleft = blackRowsYleft+1;
    end
end
for i=blackRowsYleft:size(binImage2,2)
    if sum(binImage2(:,i) == 255) == 0
        blackRowsYright = blackRowsYright+1;
    end
end

imageCropped = imro(blackRowsXtop:size(binImage2,1)-blackRowsXbottom,...
    blackRowsYleft:size(binImage2,2)-blackRowsYright);
subplot(122);



imshow(imageCropped);

%% Check if the picture is upside down
black_total_top_half =length(find(imageCropped(1:round(end/2),:) < 100));
black_total_bottom_half =length(find(imageCropped(ceil(end/2):end,:) < 100));

if black_total_bottom_half > black_total_top_half
    imageCropped = imrotate(imageCropped,180,'crop');
end

% filter_int=@(n)1/(n^2)*ones(n);
% imageCropped2 = imageCropped;
% imageCropped2 = binarizeImage(100, imageCropped2);
% imageCropped2(1:end,1:end)=255;
% binImageCropped =  binarizeImage(100, imageCropped);
% 
% [y,x] = find(binImageCropped == 0);
% p = [y x];
% 
% a = size(binImageCropped);


%% Create new Images
% imwrite(imro,'Testimage4.tif');

%% Functions
function binImage =  binarizeImage(threshold, image)
    binImage = image;
    binImage(image<(threshold)) = 0;
    binImage(image>(threshold)) = 255;
end

