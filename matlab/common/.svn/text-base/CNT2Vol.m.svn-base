function binI= CNT2Vol(pts,dims)
%Read 3d points from the contour file (from  MSKCC), scan convert into
%the binary volume
%Input:  pts --- contour points , dimension :3*N
%        dims --- image dimensions
%Output: binI ---the binary image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pts = ceil(pts);% from image coordinates to matlab coordinates
binI = zeros(dims);

z_curr = pts(3,1);
slice_x =[];% pts(1,1);
slice_y = [];%pts(2,1);

[X,Y]= meshgrid(1:dims(1), 1:dims(2));

for i = 1: size(pts,2)
    
    z = pts(3,i);
    
    if z == z_curr
        slice_x = [slice_x,pts(1,i)];
        slice_y = [slice_y,pts(2,i)];
    else
        %close the polygon
        slice_x=[slice_x,slice_x(1)];
        slice_y=[slice_y,slice_y(1)];
        
        IN = inpolygon(X(:),Y(:),slice_x,slice_y);
        if (z_curr>0)
        binI(:,:,z_curr) = reshape(IN,[dims(1), dims(2)])';
        end
        %% crucial!!!  x-y swtich
        
        z_curr = z;
        slice_x = pts(1,i);
        slice_y = pts(2,i);
        
    end
    
end




%% load byu
% X = pts';
% [K,v]= convhulln(X);
% figure(2);trisurf(K,X(:,1),X(:,2),X(:,3));
%
% [x,y,z]= meshgrid(1:512,1:512,1:99);
% I = [x(:),y(:),z(:)];
% bI = tsearch(X,K,I);


% byu = BYU();
% byu = set(byu,'nVertices',size(pts,2));
% byu = set(byu,'vertices',{pts'});
% byu = set(byu,'adjacency',K);
% byu= set(byu,'nPolygons',size(K,1));
% byu= set(byu,'nConnectors', 3*size(K,1));
% byu = set(byu,'nParts',1);
%
% write(byu,byuFN);