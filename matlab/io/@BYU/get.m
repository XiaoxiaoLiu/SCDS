function [val] = get(byu, property, index)
%
% function [val] = get(byu, property, index)
%
% property could be one of:
% nVertices
% nPolygons
% nConnectors
% adjacency
% polygons
% filenames
% vertices
% neighbors
% inPolygons
% numFiles
% in which case the corresponding property is returned.
%
% Certain properties can be indexed such as vertices. If index is specified
% then instead of returning byu.vertices, this method returns
% byu.vertices{index}
% This makes for simpler code outside.
%
switch (property)
	case 'nVertices'
		val = byu.nVertices;
	case 'nPolygons'
		val = byu.nPolygons;
	case 'nConnectors'
		val = byu.nConnectors;
	case 'adjacency'
		val = byu.adjacency;
	case 'polygons'
		val = byu.polygons;
	case 'filenames'
		val = byu.filenames;
	case 'vertices'
		if( nargin == 3)
			val = byu.vertices{index};
		else
			val = byu.vertices;
		end
	case 'neighbors'
		val = byu.neighbors;
	case 'inPolygons'
		val = byu.inPolygons;
	case 'numFiles'
		val = numel(byu.filenames);
	otherwise
		error(['BYU::get(): Unknown property : <' property '>']);
end

return;
