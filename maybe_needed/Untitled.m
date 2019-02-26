rgb = imread('peppers.jpg');
imshow(rgb);
I = rgb2gray(rgb);
hold on
h = imshow(I); % Save the handle; we'll need it later
hold off

set(h, 'AlphaData', alpha_data);

