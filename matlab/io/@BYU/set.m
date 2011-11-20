function [byu] = set(byu, varargin)
%
% function [byu] = set(byu, [property, value]* )
%
% property could be one of:
% nVertices
% nPolygons
% nConnectors
% adjacency
% polygons
% filenames
% vertices
% in which case the corresponding property is set to the passed in value.
%

for i = 1:2:length(varargin)
	val	= varargin{i+1};
	switch (varargin{i})
		case 'nVertices'
			byu.nVertices = val;
		case 'nPolygons'
			byu.nPolygons = val;
		case 'nConnectors'
			byu.nConnectors = val;
		case 'adjacency'
			byu.adjacency = val;
		case 'polygons'
			byu.polygons = val;
        case 'polyList'
            byu.polyList = val;
		case 'filenames'
			byu.filenames = val;
		case 'vertices'
			byu.vertices = val;
        case 'neighbors'
            byu.neighbors = val;
        case 'inPolygons'
            byu.inPolygons = val;
        case 'nParts'
            byu.nParts = val;
	otherwise
		error(['BYU::set(): Unknown property : <' varargin{i} '>']);
	end
end

return;
