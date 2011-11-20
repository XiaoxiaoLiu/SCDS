function [error] = distanceH(h,H,dims, spacing)
%PURPOSE: calcualte the distance  measure between two  displacement vector fields
%------------------------------------------------------------------------



if nargin <3
    s=size(h);
    dims=s(2:end);
end


if nargin <4
    spacing =[1 1 1];
end

h = reshape(h, [3 dims]);

H = reshape(H,[3 dims]);



%% mask the image for ROI region

% x_range = 60:460;
% y_range = 140:380;
% z_range = 1:dims(3);
% 
% h = h(1:3,x_range,y_range,z_range);
% H = H(1:3,x_range,y_range,z_range);


%% displacement vector field distance/distance
distance = abs(h - H);
distance_n = sqrt((spacing(1)*distance(1,:)).^2 + (spacing(2)*distance(2,:)).^2 +(spacing(3)* distance(3,:)).^2);

error.mean = mean(distance_n);
error.quan99 = quantile(distance_n,.99);
error.quan95 = quantile(distance_n,.95);
error.max = max(distance_n);



%
% %% check the motion ROIs
% norm_H = sqrt((spacing(1)*H(1,:)).^2 + (spacing(2)*H(2,:)).^2 +(spacing(3)* H(3,:)).^2);
%
% error_motion = distance_n(norm_H(:) >5); % >5mm, 2 pixel
%
% error.motion = error_motion; % 12MB, saved for later analysis
% error.motion_mean = mean(error_motion);
% error.motion_max = max(error_motion);

