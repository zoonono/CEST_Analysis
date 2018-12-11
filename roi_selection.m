%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 24.11.2017                                                              %
%                                                                         %
% ROI selection in an image (preferably the anatomical)                   %
%                                                                         %
%  Input: directory                                                       %
%         size_x                                                          %
%         size_y                                                          %
%         FOV                                                             %
%         image                                                           %
%         slice_thickness                                                 %
%         slice_index (for now 1, at the begining of the program)         %
%         directory_images                                                %
% Output: par_mask.rois                                                   %
%         par_mask.sum_roi                                                %
%         par_mask.number_of_rois                                         %
%         par_image.image_rois                                            %
%         par_image.image_sum_roi                                         %
%         par_image.image_rois_in_grey_scale                              %
%         par_image.image_sum_roi_in_grey_scale                           %
%         par_image.fig_roi_image                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[par_mask, par_image] = roi_selection(size_x, size_y, fov,image, ...
                                              slice_thickness, slice_index,...
                                              directory_images) 

flag_roi = questdlg('Do you have a mask with saved ROIs?','','Yes','No','No');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there IS a previous ROI file   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(flag_roi,'Yes')
  cd(directory_images)
  % User interface to find previously created ROI path and ROI:
  [saved_mask,path]=uigetfile('*.img','Saved ROI files'); 
  cd(path);
  fid=fopen(saved_mask,'r'); % opens roi file
  if fid <0
     msgbox('file not found!'); % error mesage if file not found
     return
  end
  
  %this 'name' will be used for saving of the ROI image outside the if statement 
  name = saved_mask(1:size(saved_mask,2)-4);

  % Read the matrix size of the mask and the number of ROIs
  [DIM, ~, ~, ~, ~, ~, ~] = spm_hread([name '.hdr']); 
  size_mx = DIM(2);
  size_my = DIM(1);
  number_of_rois = DIM(3);
  
  % Open and store the already created ROI-mask in 'mask_roi'.
  % 'mask_roi' is a stack of n mask images,
  % the n-th mask image is the mask of the n-th ROI.
  for k=1:number_of_rois   % number of rois
      for i=1:size_mx     % x axis
	      for j=1:size_my % y axis
	          mask_rois(i,j,k)=fread(fid,1,'uint8');
          end
      end
  end
  fclose(fid);
   
  % ROI-mask image with all ROIs together in one image (and not in stack)
  % if there is only one ROI, it will be stored in sum_roi too.
  mask_sum_roi = zeros(size_mx, size_my);
  for i=1:number_of_rois
      mask_sum_roi = mask_sum_roi + mask_rois(:,:,i);
  end
  
  % Display figure 
  fig_roi_image = figure;
  imshow(image, [min(image(:)), max(image(:))] );
  colormap gray
  freezeColors
  hold on
  imshow(mask_sum_roi/3, [min(mask_sum_roi(:)), max(mask_sum_roi(:))]);
  colormap hot
  alpha 0.5
  title('PREVIOUS SELECTED ROI') 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there is NO previous ROI file that can be used %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else    
   number_of_rois = str2double(char(inputdlg ...
                   ('How many ROIs do you want to select?','',1 )));
   %zeros
   mask_sum_roi = zeros(size_x, size_y);
   mask_rois    = zeros(size_x, size_y, number_of_rois);
   
   % display figure (MAYBE I SHOULD DO IT OUTSIDE THE IF statement)
   fig_roi_image = figure;
   imshow(image,[min(image(:)) max(image(:))]);
   if size_x - size_y > 0 | size_x - size_y < 0
      axis square
   end
   title('Manually selected ROIs')
   %select the ROIs manually
   for i=1:number_of_rois %loop over all ROIs!
       if i>1
          title(['Select', num2str(i), 'ROI manually'])
       end
       % 'roipoly' is the tool to draw the ROI,
       % x_p is the x coordinate of the point of the polygonal ROI,
       % y_p is the y coordinate of the point of the polygonal ROI.
       [roi,x_p,y_p]=roipoly; 
       mask_rois(:,:,i) = double(roi);
       % Connects the coordinates of the ROI with a line
       for n=1:(size(x_p,1)-1)
           h_line=line([x_p(n) x_p(n+1)],[y_p(n) y_p(n+1)]);
           set(h_line,'color','g');
           set(h_line,'LineWidth',2);
       end
       % Finds the 'center' of the ROI
       cx=mean(x_p);
       cy=mean(y_p);
       % Display the numbering of the ROI in the center of the ROI
       text(cx, cy, num2str(i), 'FontSize', 16, 'Color','red') 
       
       % ROI-mask image with all ROIs together in one image (and not in stack)
       mask_sum_roi = mask_sum_roi + mask_rois(:,:,i);
   end
   clear x_p y_p cx cy; 
   %this is to be used later to save the ROI image
   name = ['manually_selected_ROI_slice',num2str(slice_index),'_',...
            datestr(datetime,'ddmmyyyyHHMMSS')] ;   
end

cd(directory_images)
%Write the header file .hdr using the spm99 package.
%spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP)
% P       - filename 	     (e.g 'spm' or 'spm.img')
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [sec])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - [i j k] of origin  (default = [0 0 0])
% DESCRIP - description string (default = 'spm compatible')
% s       - number of elements successfully written (should be 348)
%That is the ROI mask, not the roi image.

%'name' was defined within elseif command where there was NO previous saved ROI 
% or was kept the same in the if command where there was a previous ROI 
[s] = spm_hwrite(name, [size_y size_x number_of_rois],...
                 [(fov(2)/size_y) (fov(1)/size_x) (slice_thickness*10)],...
                 1,2,0,[0 0 0]);
disp(['Saved ',name,'.hdr in',pwd]);
  
%save the image with roi mask in .img format in folder images_and_figures
fid=fopen([name, '.img'],'w');
for k=1:number_of_rois
    for i=1:size_x
        for j=1:size_y
            fwrite(fid,mask_rois(i,j,k),'uint8');
        end
    end
end
fclose(fid);
disp(['Saved ',name,'.img in ',pwd]);

%Save the figure of the image with the selected ROI in folder images_and_figures
title('ROI image');
saveas(fig_roi_image,[name,'.tif']);
disp(['Saved ',name,'.tif in ',pwd]);

for i = 1 : number_of_rois 
    image_roi(:,:,i) = immultiply(mask_rois(:,:,i),image); %masks image
end
image_roi_in_grey_scale = mat2gray(image_roi);
image_sum_roi = immultiply(mask_sum_roi,image);
image_sum_roi_in_grey_scale = mat2gray(image_sum_roi);

par_mask.number_of_rois = number_of_rois;
if number_of_rois > 1
  par_mask.rois = mask_rois;
end
par_mask.sum_roi = mask_sum_roi;
par_image.image_rois = image_roi;
par_image.image_sum_roi = image_sum_roi;
par_image.image_rois_in_grey_scale = image_roi_in_grey_scale;
par_image.image_sum_roi_in_grey_scale = image_sum_roi_in_grey_scale;
par_image.fig_roi_image = fig_roi_image;

clear DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP
end