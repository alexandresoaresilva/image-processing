%imread('paris_879.jpg');
n = 400;
n_half = n/4;
n_quart = n_half/4;
n_eight = n_quart/2;

n_quart = round(n_quart);
n_eight = round(n_eight);
n_half = round(n_half)

A = zeros(n);
% square
A(n_quart:(n_quart+n_half),n_quart:(n_quart+n_half)) = 1;
A(n_quart:(n_quart+n_half),n_quart:(n_quart+n_half)) = 1;

% two triangles
A(n_quart:(n_quart+n_half),n_quart:(n_quart+n_half)) = 1;
A(n_quart:(n_quart+n_half),n_quart:(n_quart+n_half)) = 1;

imshow(A);
r = round(n_quart/2);
SE = strel('square',r)
J = imerode(A,SE);
imshow([A,J]) 