function [X]	= slicez(byu,z)
%
% function [X]	= slicez(byu,z)
%
% Slices the tileset in byu by the plane z and returns the <x,y,z>
% co-ordinates of the points at the slices.
% Right now, the returned points are not ordered in any way.
%

vertices= byu.vertices{1};
faces   = byu.adjacency;

selector = repmat(false,[size(faces,1) 1]);

for fi = 1:size(faces,1)
	plusc = sum(vertices(faces(fi,:),3) >= z);
	minusc = sum(vertices(faces(fi,:),3) <= z);
	selector(fi) = ( plusc ~= 0 && minusc ~= 0 );
end

faces = faces(selector,:);

pts = 0;
X	= zeros(size(faces,1)*2,3); % Upto 2 per face.

for fi = 1:size(faces,1)
	for si = 1:size(faces,2)
		v0 = vertices(faces(fi, si),:);
		v1 = vertices(faces(fi, mod(si,size(faces,2))+1 ),:);
		if( v0(3) == v1(3) && v0(3) == z )
			% This entire segment is in the plane ...
			% add both end pts, and go to the next triangle.
			% We don't deal with cases where the entire triangle might be
			% lying inside the plane ...
			X(pts+1,:) = v0;
			X(pts+2,:) = v1;
			pts = pts + 2;
			break;
		end
		t = (z - v0(3)) / (v1(3) - v0(3));
		if( t >= 0 && t <= 1 )
			% This segment intersects the plane.
			X(pts+1,:) = [ (1-t) * v0(1:2) + t * v1(1:2), z ];
			pts = pts + 1;
		end
	end
end

X = X(1:pts,:);

return;
