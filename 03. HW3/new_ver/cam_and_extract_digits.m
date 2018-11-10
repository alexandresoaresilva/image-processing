clc, clear, close all
addpath('VLFEATROOT')

run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose
%for the closing operation
AVG_FILTER_SIZE = 3; %for extracting digits
DIGIT_SIZE = 120; %image size
addpath('numbers'); %where images/mat files are stored are stored
%% Loading and pre-processing images
%image A (reference)
% a = load('image_Templates.mat');
% Ia_orig = a.image_Templates{8};
a = load('number_imgs.mat');

Ia_orig = a.number_imgs{1}(72);

[Ia, Ia_bin] = extract_digits(Ia_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
Ia = Ia{1};
Ia_bin = process_bin_num(Ia_bin{1});
% to compare (camera can be used, as well)
Ib_orig = imread('7_2.jpg');
[Ib, Ib_bin] = extract_digits(Ib_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
Ib = Ib{1};
Ib_bin = process_bin_num(Ib_bin{1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(221);
imshow([Ia,Ib],[]);
title('Ia X  Ib');
hold on
%% SIFT
[fa,da] = run_SIFT(histeq(Ia),0);
[fb,db] = run_SIFT(Ib, size(Ia,2));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(222);
imshow([Ia_bin,Ib_bin]);
title('Ia bin X Ib bin');
[fa_bin,da_bin] = run_SIFT(Ia_bin,0);
x_offset = size(Ia_bin,2);
[fb_bin,db_bin] = run_SIFT(Ib_bin,x_offset);
%%  SIFT MATCHES ///////////////////////////////////////////////////////////
[matches_normal, scores_normal] = vl_ubcmatch(da, db) ;
m  = median(scores_normal);
std_n = std(scores_normal);
s = 0 % .5*std_n;
index_m = find(scores_normal >= (m + s) );
selec_scores_norm = scores_normal(index_m);
selec_matches_norm = matches_normal(:,index_m);

[matches_bin, scores_bin] = vl_ubcmatch(da_bin, db_bin) ;
m_bin  = median(scores_bin);
std_bin = std(scores_normal);
index_m_bin = find(scores_bin >= (m_bin + 2.5*std_bin));
selec_scores = scores_bin(index_m_bin);
selec_matches = matches_bin(:,index_m_bin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(223);
imshow([Ia,Ib]);
title('Ia X  Ib');
draw_lines(da,fa,db,fb,selec_matches_norm,selec_scores_norm,x_offset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(224);
imshow([Ia_bin,Ib_bin]);
title('Ia bin X Ib bin');

draw_lines(da_bin, fa_bin, db_bin, fb_bin,...
    selec_matches, selec_scores, x_offset);

%% 
function [f,d] = run_SIFT(I, x_offset)
    % A frame is a disk of center f(1:2), scale f(3) and orientation f(4)
    [f,d] = vl_sift(single(I)) ;
    perm = randperm(size(f,2)) ;
    sel = perm;
    %two circles because one is larger thant the other and thus, becomes
    %   a little contour (plotting circles with angles
    % the offset is present for the second picture
    f_for_plotting = f;
    f_for_plotting(1,:) = f_for_plotting(1,:) + x_offset;
    
    h1 = vl_plotframe(f_for_plotting(:,sel)) ;
    h2 = vl_plotframe(f_for_plotting(:,sel)) ;
    set(h1,'color','k','linewidth',3) ;
    set(h2,'color','y','linewidth',2) ;
    
    %plots square frame 4x4 with gradients
    h3 = vl_plotsiftdescriptor(d(:,sel),f_for_plotting(:,sel)) ;
    set(h3,'color','g') ;
end

function draw_lines(da,fa,db,fb,matches,scores,x_offset)
    xa = fa(1,matches(1,:)) ;
    xb = fb(1,matches(2,:)) + x_offset;
    ya = fa(2,matches(1,:)) ;
    yb = fb(2,matches(2,:)) ;
     
    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'b') ;

    vl_plotframe(fa(:,matches(1,:)));
    fb(1,:) = fb(1,:) + x_offset;
    vl_plotframe(fb(:,matches(2,:)));
    colormap gray
    axis image off ;
    hold off

end

% function 
% 
% end

function I = process_bin_num(I_bin)
    LPF = avg_filt(3);
    I_bin = uint8(I_bin*255);
    % I1 = I;
    I = conv2(I_bin,LPF,'valid');
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end