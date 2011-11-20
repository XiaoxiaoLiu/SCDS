function [cog]= cog(imFN)
%Calcualte the center of gravity of binary volume, in world coordinate
%Input:  binary image file name
%Output: the 1 by 3 coordinate of the COG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[im, orig, spacing]=loadMETA(imFN);
dims = size(im);
x=[];y=[];z=[];

for i=1:dims(1)
    for j= 1:dims(2)
        for k=1:dims(3)
            if im(i,j,k)>0
                x=[x,i];
                y=[y,j];
                z=[z,k];
            end
        end
    end
end

cog(1)=mean(x);
cog(2)=mean(y);
cog(3)=mean(z);


cog = cog.*spacing + orig;