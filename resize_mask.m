%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 24.11.2017                                                              %
%                                                                         %
% Resizes a mask to match the image size of the image to be masked        %
%                                                                         %
% Input:  mask_size_x                                                     %
%         mask_size_y                                                     %
%         size_x                                                          %
%         size_y                                                          %
%         mask_sum_roi_before_resize                                      %
%         mask_roi_before_resize                                          %
%         number_of_roi                                                   %
% Output: mask                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[mask, flag_resize] = ... 
                        resize_mask(mask_size_x, mask_size_y, ...
                                    image_size_x, image_size_y, ...
                                    mask_before_resize)

%check if the mask image and the image to be masked have the same size
if mask_size_x-image_size_x>0 | mask_size_x-image_size_x<0 | ...
   mask_size_y-image_size_y>0 | mask_size_y-image_size_y<0

   flag_resize = 1;
   % Bilinear interpolation: the output pixel value is a weighted average 
   % of pixels in the nearest 2-by-2 neighborhood
   mask = imresize(mask_before_resize,[image_size_x, image_size_y],'bilinear');   
   % After resize, not all mask values are 0 or 1, 
   % we have values = 0.5, 0.6 ...
   % Solution: values < 0.5 become zero, values > 0.5 become 1.
   
   for i=1:image_size_x
       for j=1:image_size_y
           if mask(i,j) > 0.5   
              mask(i,j) = 1;
           else
              mask(i,j) = 0;
           end
       end
   end
   end
end