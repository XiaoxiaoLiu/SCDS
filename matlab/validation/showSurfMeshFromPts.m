function[]=showSurfMeshFromPts(pts,displaySizes)
% PURPOSE:
%   INPUT: 
%  OUTPUT: 
%    Note: 
%  AUTHOR: Xiaoxiao liu, CS of UNC-CH
%    DATE: 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pts N*3

%[t]=MyRobustCrust(pts);
[t]= MyCrustOpen(pts);

%% plot the points cloud
axis equal
trisurf(t,pts(:,1),pts(:,2),pts(:,3),'facecolor','c','edgecolor','b')%plot della superficie

view(180,-20);
axis vis3d
axis(displaySizes);
axis equal
title('Output Triangulation','fontsize',14)




