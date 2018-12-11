%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 14.01.2018                                                              %
%                                                                         %
% Area under the curve                                                    %
%                                                                         %
%  Input: x (ppm_values)                                                  %
%         y (values_vector (raw z spectrum of ROI or pixel))              % 
%         st_ppm (middle point of integral)                               %
%         width_integral                                                  %
%         endpoint_auc_from_zero                                          %
%                                                                         %
% Output: .auc_around_mid_point                                           %
%         .auc_from_zero                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = auc(x, y, mid_point_of_integral, ...
                           width_integral,end_point_auc_from_zero)

%calculate endpoints of integral
left_endpoint  = mid_point_of_integral - width_integral/2;
right_endpoint = mid_point_of_integral + width_integral/2;
%make sure the endpoints equal to values of the variable "ppm_values", 
%otherwise go to the closest values of the "ppm_values" vector
[~, index_left] = min(abs(x(:)-left_endpoint));
[~, index_right] = min(abs(x(:)-right_endpoint));
[~, index_end_point_auc_from_zero] = min(abs(x(:)- end_point_auc_from_zero));
[~, index_zero] = min(abs(x(:)-0.0));

%change endpoints and mid point to the closest x value
left_endpoint = x(index_left);
right_endpoint = x(index_right);
zero_endpoint = x(index_zero);
end_point_auc_from_zero = x(index_end_point_auc_from_zero);

%intergral within [midpoint-width/2, midpoint+width/2]
integral_x_points = x(index_left:index_right);
integral_y_points = y(index_left:index_right);                    
auc_around_mid_point = trapz(integral_x_points, integral_y_points);

%intergral within [0, end_point_auc_from_zero] 
integral2_x_points = x(index_zero:index_end_point_auc_from_zero);
integral2_y_points = y(index_zero:index_end_point_auc_from_zero);                    
auc_from_zero = trapz(integral2_x_points, integral2_y_points);

%put variables is "parameters." struct
parameters.left_endpoint = left_endpoint;
parameters.right_endpoint = right_endpoint;
parameters.zero_endpoint = zero_endpoint;
parameters.end_point_auc_from_zero = end_point_auc_from_zero;
parameters.auc_around_mid_point = auc_around_mid_point;
parameters.auc_from_zero = auc_from_zero;
end