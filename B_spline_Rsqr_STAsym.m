%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 05.12.2017                                                              %
%                                                                         %
% Cubic smoothing spline!                                                 %
% Fit the average values of the selected roi with the B-splines           %
% x is the ppm frequencies at which the images were taken                 %
% y are the normalized average values of the images in the selected range %
%                                                                         %
%  Input: ppm_values                                                      %
%         values_vector (raw z spectrum of ROI or pixel)                  % 
%         regularization_factor                                           %
%         auto_or_custom                                                  %
%         sign_st (in the ppm of the molecule is positive of negative     %
%                                                                         %
% Output: .Rsqr_of_z_spectrum                                               %
%         .cs                                                             %
%         .z_spectrum                                                     %
%         .ppm                                                            %
%         .st_positive                                                    %
%         .st_negative                                                    %
%         .non_interp_normalized_data                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = B_spline_Rsqr_STAsym(ppm_values, values_vector, regularization_factor, ...
                               auto_or_custom,sign_st)

normalized_vector = values_vector / max(values_vector);                           
                           
%automatic regularization factor (p)
if strcmp(auto_or_custom, 'auto') || strcmp(auto_or_custom, 'automatic')||...
   strcmp(auto_or_custom, 'Auto') || strcmp(auto_or_custom, 'Automatic')||...
   strcmp(auto_or_custom, 'AUTO') || strcmp(auto_or_custom, 'AUTOMATIC')||...
   strcmp(auto_or_custom, 'a') || strcmp(auto_or_custom, 'A')
   
  [~,p]=csaps(ppm_values,normalized_vector); % p is the calculated regularization factor
  %csaps(x,y,p,xx,w), uses the new calculated regularization factor
  cs = csaps(ppm_values,normalized_vector, p, [], []); 
  %finds the minimun of the interpolation fucntion
  [~,min_site]=fnmin(cs);
  
%custom regularization factor  
elseif strcmp(auto_or_custom, 'custom') || strcmp(auto_or_custom, 'Custom')||...
       strcmp(auto_or_custom, 'CUSTOM') || strcmp(auto_or_custom, 'C')||...
       strcmp(auto_or_custom, 'c')
   
  cs = csaps(ppm_values,normalized_vector, regularization_factor, [], []); 
  [~,min_site]=fnmin(cs);  
end

% uses the interpolated values of the B spline to get the y values at steps
% of 0.1 ppm by centering zero on the absolute minimum value of B spline.
step=0.1;
i_max=round(max(ppm_values)/step);
% first, predefine variables
positive_interpolated = zeros(1,i_max);
negative_interpolated = zeros(1,i_max);
ppm_positive_interpolated = zeros(1,i_max);
ppm_negative_interpolated = zeros(1,i_max);

%then, calculate z-spectrum values and ppms for positive side
for i = 1:i_max
    positive_interpolated(i)  = (fnval(cs, (i*step + min_site)));
    ppm_positive_interpolated(i) = step * i;
end
%also calculate z-spectrum value and ppm for zero
st_zero_interpolated = (fnval(cs, min_site));
ppm_zero_interpolated = 0.0 ;
%and last but not least, calculate z-spectrum values and ppms for negative side
for i = 1:i_max
    negative_interpolated(i)  = (fnval(cs,(-i*step + min_site)));
    ppm_negative_interpolated(i) = step * (-i);
end

%calculate ST Asymmetry curve
% calculates the value of ST% as 1-[(Is/Io)*100]^sign_st
for i = 1:i_max
    STAsym (i) = 1 - (positive_interpolated(i)/negative_interpolated(i))^sign_st;
end
%put ZERO too
STAsym = horzcat(0.0,STAsym);

%also put ZERO in zero_positive_interpolated
positive_interpolated = horzcat(st_zero_interpolated,positive_interpolated);
ppm_positive_interpolated = horzcat(ppm_zero_interpolated,ppm_positive_interpolated);

%put interpolated spectrum and ppm values in one variable each
z_spectrum_interpolated = horzcat(flip(negative_interpolated),...
                                  positive_interpolated);
ppm_interpolated = horzcat(flip(ppm_negative_interpolated),...
                           ppm_positive_interpolated);

                       
% calculation of R SQUARED, Dario Longo's (UNITO) method to check the goodness 
% of the interpolated curve to pass from i experimental points:  
% R^2= 1 - Sum(Yi-Yi,b)^2/Sum(Yi-Ymean)^2
% where:
%       Yi: the experimental points             (values_vector)
%     Yi,b: the interpolated points             (fnval(cs,ppm_values))
%    Ymean: the mean of the experimental points (mean(values_vector))
%
   numerator=(normalized_vector-fnval(cs,ppm_values)).^2;  %(arithmitis)
   numerator=sum(numerator(:));                            %(arithmitis)
   denominator=(normalized_vector-mean(normalized_vector)).^2; %(paronomastis)
   denominator=sum(denominator(:));                            %(paronomastis)
   Rsqr_of_z_spectrum = 1-(numerator/denominator);                         

%save variables in "parameters." struct   
parameters.Rsqr_of_z_sprectrum = Rsqr_of_z_spectrum;   
parameters.cs = cs;
parameters.z_spectrum = z_spectrum_interpolated;
parameters.ppm = ppm_interpolated;
parameters.ppm_positive = ppm_positive_interpolated;
parameters.st_positive = positive_interpolated;
parameters.st_negative = negative_interpolated;
parameters.non_interp_normalized_data = normalized_vector;
parameters.STAsym = STAsym;
end
