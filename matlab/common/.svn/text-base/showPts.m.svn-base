function[]= showPts(p,displaySizes,color, orientation)
% PURPOSE:
%   INPUT: 
%  OUTPUT: 
%    Note: 
%  AUTHOR: Xiaoxiao liu, CS of UNC-CH
%    DATE: 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'position',[0,0,1280,800]);
%subplot(1,2,1)

if nargin<3
    color='b';
end

if nargin<4
    orientation='RAI';
end

switch  orientation 
case 'RAI'
        scatter3(p(:,1),p(:,2),p(:,3),25,color,'filled');
case 'RAS'
        maxZ = max(p(:,3));  
        minZ= min(p(:,3)); 
        scatter3(p(:,1),p(:,2),maxZ-p(:,3)+minZ,25,color,'filled');
end


 
axis(displaySizes);

view(-10,20);



