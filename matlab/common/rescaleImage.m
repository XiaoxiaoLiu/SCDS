function im = rescaleImage(fn,range)
%Rescale Image intensity range to [range(1), range(2)]
%----------------------------------------------------

im = loadMETA(fn);
minIM = min(im(:));
maxIM = max(im(:));
im = (im -minIM)/(maxIM-minIM);
im = im *( range(2)-range(1)) + range(1);