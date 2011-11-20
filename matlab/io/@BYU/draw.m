function [P] = draw(byu, color, index)
%
% function [P] = draw(byu, color, index)
%
% Draws a tileset in the current figure as a patch object and returns the
% handle in P.
% It optionally also uses color for vertex coloring and index to select
% which of the tilesets to draw. If index is not specified, then 1 is
% assumed.
%

if (nargin <= 2)
	index = 1;
end

if (nargin <= 1)
	color = [1.0 1.0 1.0];
end

FV = struct('vertices', byu.vertices{index}, 'faces', abs(byu.adjacency) );
%'facevertexcdata', color );

P = patch(FV);
set(P, 'FaceLighting', 'phong');
set(P, 'BackFaceLighting', 'lit');
set(P, 'EdgeColor', 'none');
set(P, 'FaceColor', color(1:3) );
if( numel(color) == 4 )
	set(P, 'FaceAlpha', color(4) );
end
set(P, 'SpecularExponent', 5)
set(P, 'SpecularStrength', 0.15)

return;
