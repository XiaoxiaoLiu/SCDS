function movie2subplot( M, i, j)

nPlots = i*j;
step = floor(length(M)/nPlots);
frames = 1:step:length(M);

clf;
for k=1:nPlots
    subplot( i,j,k);
    imagesc( M(frames(k)).cdata );
end