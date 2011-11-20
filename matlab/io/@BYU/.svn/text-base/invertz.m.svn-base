function [byu]	= invertz(byu, flipTiles)
%
% function [byu]	= invertz(byu, flipTiles)
%
% Just negates all the z co-ordinates, does not do any tile flips right
% now.
%

for i = 1:length(byu.vertices)
	vertices = byu.vertices{i};
	vertices(:,3) = -vertices(:,3);
	byu.vertices{i} = vertices;
end

return;
