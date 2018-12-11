%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 25.11.2017                                                              %
%                                                                         %
% Segmentation of an image using intensity threshold.                     %
% In the context of the CEST_analysis code, the segmentation is done      %
% on the anatomical image.                                                %
%                                                                         %
%  Input: size_x                                                          %
%         size_y                                                          %
%         image (image_2dseq)                                             %
% Output: mask.intensity                                                  %
%         par_images.intensity_masked_image                               %
%         par_images.intensity_masked_in_grey_scale                       %
% and conditionally:  mask.user_segmented                                 %
%                     par_image.user_masked_image                         %
%                     par_image.user_masked_in_grey_scale                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[mask,par_image] = mask_from_image(size_x,size_y,image)

% 'thresh_tool' opens a dialog window where the anatomical image gets
% segmented (masked) by changing the pixel intensity threshold.
% That way the noise from outside the brain is removed in this first step.
[level,segmentation(:,:)]=thresh_tool(mat2gray(image(:,:)),'gray');

segmented_mask = segmentation;     % segmented_mask: binary image!
segmented_image = immultiply(segmentation,image); % makes zero the non selected pixels
segmented_in_grey_scale = mat2gray(segmented_image);

mask.intensity = segmented_mask;
par_image.intensity_masked_image = segmented_image;
par_image.intensity_masked_in_grey_scale = segmented_in_grey_scale;

segmented_image_figure = figure;
imshow(segmented_image, [min(segmented_image(:)), max(segmented_image(:))]);
     if size_x ~= size_y
            axis square
     end
title('SEGMENTED IMAGE');

%% manual segmentation by the user

flag_seg = questdlg('Is the segmentation correct?','','yes','No','yes');
if strcmp(flag_seg,'No')   
    close(segmented_image_figure);
    clear image_seg;
    
    %opens a figure of the initial image to manually segment
    original_image_figure = figure; 
    imshow(image,[min(image(:)), max(image(:))]); 
    title('ORIGINAL IMAGE TO BE MANUAL SEGMENTED')
    
    % User interface manual segmentation; returns a binary image!
    segmented_mask = roipoly; 
    close(original_image_figure);

    segmented_image = immultiply(segmented_mask,image); %masks image 
    
    %opens a figure of the user segmented image
    segmented_image_figure = figure; 
    imshow(segmented_image, [min(segmented_image(:)), max(segmented_image(:))]);
    if size_x-size_y>0 | size_x-size_y<0
       axis square
    end
    title('MANUALLY SEGMENTED IMAGE')

    segmented_in_grey_scale = mat2gray(segmented_image);

    mask.user_segmented = segmented_mask;
    par_image.user_masked_image = segmented_image;
    par_image.user_masked_in_grey_scale = segmented_in_grey_scale;
end
clear flag_seg
end