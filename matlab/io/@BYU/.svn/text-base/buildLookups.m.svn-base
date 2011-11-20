function [neighbors, polygons] = buildLookups(byu)
%
% function [neighbors, polygons] = buildLookups(byu)
%
% Build lookup tables for neighbor information and polygon containment.
%
neighbors = cell([byu.nVertices,1]);
polygons  = cell([byu.nVertices,1]);

for ipoly=1:byu.nPolygons
	poly = byu.adjacency(ipoly,:);
	pv   = abs(poly(end));
	for iv = 1:length(poly)
		v  = abs(poly(iv));
		neighbors{v}  = [ neighbors{v}, pv ];
		neighbors{pv} = [ neighbors{pv}, v ];
		pv = v;
		polygons{v}   = [ polygons{v}, ipoly ];
	end
end

for i = 1:byu.nVertices
	neighbors{i} = int32(unique(neighbors{i}));
	polygons{i}  = int32(unique(polygons{i}));
end