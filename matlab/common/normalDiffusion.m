function normalDiffusion(normals, byu, img,outfile)
%By Rohit
%Fitting the m-rep model to the binary image
%---------------------------------------------


% pad image
padding = 1;

%%
figure

% some black magic
az = -121.5;
el = 6;

N = get(normals,'vertices',1);
pts = get(byu,'vertices',1);
neighbors = get(byu, 'neighbors');

% converting from world  to model coordinates
%R:\proj\lung3\mskcc_20070612_P41\Patient1\adapt\image\binary\right\right-lung-bin-cinePhase50.mhd
max_extent = max(size(img)-2*padding);

% spacing = [0.97656 0.97656 0.97656];
% origin  = [ -24.9512 41.1133 10.724];
%pts = (pts - repmat(origin, [size(pts,1) 1]) ) ./ repmat(spacing*max_extent, [size(pts,1) 1]);

byu = set(byu,'vertices',{pts});

clear path;
path(:,:,1) = N;

%%
for i = 1:100
    clf;
    draw(byu,[1 1 0], 1);
    camlight, lighting gouraud;
    view(3); axis equal; axis off; axis vis3d; hold on; view(az, el);
    quiver3(pts(:,1),pts(:,2),pts(:,3), N(:,1), N(:,2), N(:,3));
    pause(0.2);
    N = diffuseCoupledTrihedronNormals(N, neighbors, 10);
    path(:,:,i+1) = N;
end

%% Now fit

max_extent = max(size(img)-2*padding);
ends = (size(img)-2*padding)/max_extent;
step = [1 1 1]./(size(img)-2*padding - 1);

% +1 is because counting in matlab starts from +1
modelToImage = @(pt) pt*(max_extent) - 0.5  + 1 + padding;
lookupImage  = @(img, pt, ipt) interp3( ipt(1):ipt(1)+1, ipt(2):ipt(2)+1, ipt(3):ipt(3)+1, ...
    single(img( ipt(1):ipt(1)+1, ipt(2):ipt(2)+1, ipt(3):ipt(3)+1 ))/8.0, ...
    pt(1), pt(2), pt(3) );

for i = 1:size(pts,1)
    if( mod(i,10) == 0 )
        fprintf('%5.2f%%\n',i/size(pts,1)*100.0);
    end
    % Step in +ve direction:
    pt = pts(i,:);
    if( min(pt) < 0.0 || sum(pt > ends) > 0 )
        old_dist = +10000;
    else
        old_dist = lookupImage(img, modelToImage(pt), floor(modelToImage(pt)) );
    end
    % Choose a step size that takes at least 2 steps per voxel.
    if( old_dist > 0 )
        % We are starting from outside, use reverse direction
        nstep = -1/max_extent/5;
    else
        %> In path at 110
        nstep = 1/max_extent/5;
    end
    for t = 1:size(path,3);
        pt2 = pt + nstep * path(i,:,t);
        if( min(pt2) < 0.0 || sum(pt2 > ends) > 0 )
            dist = +10000;
        else
            dist = lookupImage(img, modelToImage(pt2), floor(modelToImage(pt2)) );
        end
        if( sign(dist) ~= sign(old_dist) )
            % Choose a weighted average between pt and pt2
            pts(i,:)    = (pt2*abs(old_dist)+pt*abs(dist))/(abs(old_dist)+abs(dist));
            break;
        end
        old_dist = dist;
        pt = pt2;
    end
end

% convert from model to world coords.

new_byu = set(byu,'vertices',{pts});


%new_byu = set(byu,'vertices',{pts .* repmat(spacing*max_extent, [size(pts,1) 1]) + repmat(origin, [size(pts,1) 1])});
%byu = set(byu,'vertices',{get(byu,'vertices',1) .* repmat(spacing*max_extent, [size(pts,1) 1]) + repmat(origin, [size(pts,1) 1])});

%%
figure
PN=draw(invertz(new_byu),[1 1 0], 1);
camlight, lighting gouraud;
view(3); axis equal; axis off; axis vis3d; hold on; view(az, el);
set(PN, 'FaceAlpha', 0.5 );
hold on
P=draw(invertz(byu),[1 0 0], 1);
set(P, 'FaceAlpha', 0.5 );

write(new_byu, outfile);