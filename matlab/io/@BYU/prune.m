function [byu]	= prune(byu, selector)
%
% function [byu]	= byuPrune(byu selector)
%
% byu - byu header and points list
% selector - a matrix representing which rows are to be selected.
%
% Returns a tileset [byu] which contains points and connectivity
% information only for the selected bits.
%
% This function does not know how to handle more than 1 part tilesets
%

%
% the new index for the vertex after the pruning.
%
renamed_vertex	= zeros([byu.nVertices,1]);
poly_selector	= repmat(true, [byu.nPolygons,1]);

nVertices       = 0;
%nConnectors     = byu.nConnectors;
for iv = 1:byu.nVertices
	if(selector(iv))
		nVertices = nVertices + 1;
		renamed_vertex(iv) = nVertices;
	else
		renamed_vertex(iv) = -1;
		polygons    = byu.inPolygons{iv};
%		nConnectors = nConnectors - length(polygons);
		for ipoly = 1:length(polygons)
			poly_selector(polygons(ipoly)) = false;
		end
	end
end

for i = 1:numel(byu.adjacency)
	old_iv  = byu.adjacency(i);
	if( old_iv < 0 )
		byu.adjacency(i) = -renamed_vertex(-old_iv);
	else
		byu.adjacency(i) = renamed_vertex(old_iv);
	end
end

byu.adjacency = byu.adjacency(poly_selector,:);

%byu.nConnectors = nConnectors;
% nConnectors is behaving weirdly.
byu.nConnectors = byu.nConnectors/byu.nPolygons * size(byu.adjacency, 1);

byu.nVertices = nVertices;
byu.nPolygons = size(byu.adjacency,1);
byu.polyList  = [ 1 byu.nPolygons ];
for i = 1:length(byu.vertices)
	byu.vertices{i} = byu.vertices{i}(selector,:);
end
% Re-build lookup tables.
[n,p] = buildLookups(byu);
byu.neighbors  = n;
byu.inPolygons = p;

return;
