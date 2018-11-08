clc, clear, close all
% addpath('VLFEATROOT')
run('vlfeat/toolbox/vl_setup');
vl_version verbose
%for the closing operation
SELECT_I = 8;

AVG_FILTER_SIZE = 3;
[file_names, file_names_char] = save_file_names_in_folder(pwd,'jpg');

%webcam setup
webcamlist_conn = webcamlist;
webcam_connected = webcam(webcamlist_conn{1});
%webcam_connected.Resolution ='640x480';

% webcam_connected.Resolution = '320x240';
prv = preview(webcam_connected);
prv.Parent.Visible = 'off';
prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
prv.Parent.Visible = 'on';
prompt = 'Press enter to capture a frame';
x = input(prompt);
I = snapshot(webcam_connected);
delete(prv);



%% Pre-Process Image
Ia_orig = I;
Ib_orig = imread('7_3.jpg');

Ia = ~imbinarize(rgb2gray(Ia_orig));
Ib = ~imbinarize(rgb2gray(Ib_orig));

Ia = medfilt2(Ia);
Ib = medfilt2(Ib);

Ia = (edge(Ia,'Canny'));
Ib = (edge(Ib,'Canny'));

se = strel('disk',2);
Ia = imdilate(Ia,se);
Ib = imdilate(Ib,se);

Ia = imbinarize(imgaussfilt(double(Ia),10));
Ib = imbinarize(imgaussfilt(double(Ib),10));


% Clip the Images
[row,col] = find(Ib == 1);
indices = [row,col];

mins = min(indices);
maxs = max(indices);

minY = mins(1);
minX = mins(2);
maxY = maxs(1);
maxX = maxs(2);

Ib = Ib(minY:maxY,minX:maxX);

[row,col] = find(Ia == 1);
indices = [row,col];

mins = min(indices);
maxs = max(indices);

minY = mins(1);
minX = mins(2);
maxY = maxs(1);
maxX = maxs(2);

Ia = Ia(minY:maxY,minX:maxX);


if size(Ib,1)+size(Ib,2) > size(Ia,1)+size(Ia,2)
    Ia = imresize(Ia,[size(Ib,1),size(Ib,2)]);
else
    Ib = imresize(Ib,[size(Ia,1),size(Ia,2)]);
end

figure(1)
subplot(121)
imagesc(Ia)
colormap gray
subplot(122)
imagesc(Ib)
colormap gray

%% SIFT 


[fa,da] = vl_sift(im2single((Ia))) ;
[fb,db] = vl_sift(im2single((Ib))) ;
% 
% interest_points = floor(length(fa)/2);
% [~,fa_sorted_indi] = sort(fa(3,:),'descend');
% fa_sorted_indi = fa_sorted_indi(1:interest_points);
% fa = fa(:,fa_sorted_indi);
% da = da(:,fa_sorted_indi);
% 
% [~,fb_sorted_indi] = sort(fb(3,:),'descend');
% fb_sorted_indi = fb_sorted_indi(1:interest_points);
% fb = fb(:,fb_sorted_indi);
% db = db(:,fb_sorted_indi);

figure(2)
subplot(121)
imshow(Ia)
perm = randperm(size(fa,2)) ;
sel = perm ;
h1 = vl_plotframe(fa(:,sel)) ;
h2 = vl_plotframe(fa(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;
h3 = vl_plotsiftdescriptor(da(:,sel),fa(:,sel)) ;
set(h3,'color','g') ;

figure(2)
subplot(122)
imshow(Ib);
perm = randperm(size(fb,2)) ;
sel = perm;
h1 = vl_plotframe(fb(:,sel)) ;
h2 = vl_plotframe(fb(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;
h3 = vl_plotsiftdescriptor(db(:,sel),fb(:,sel)) ;
set(h3,'color','g') ;

[matches, scores] = vl_ubcmatch(da,db) ;
[drop, perm] = sort(scores, 'ascend') ;
perm = perm(1:5);
matches = matches(:, perm) ;

scores  = scores(perm) ;


figure(3) ; clf ;
imagesc(cat(2, Ia, Ib)) ;
colormap gray
axis image off ;



xa = fa(1,matches(1,:)) ;
xb = fb(1,matches(2,:)) + size(Ia,2) ;
ya = fa(2,matches(1,:)) ;
yb = fb(2,matches(2,:)) ;

hold on ;
h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;

vl_plotframe(fa(:,matches(1,:))) ;
fb(1,:) = fb(1,:) + size(Ia,2) ;
vl_plotframe(fb(:,matches(2,:))) ;
colormap gray
axis image off ;

sum(scores)