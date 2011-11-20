

function[im, LM] =  circleGen(im, r, center)
dims= size(im);

dvec = repmat((center)',[1 dims]) - eyeHField(dims);
d = vFieldPointwiseNorm(dvec,2);

im( d <= r) = 1;


% lm :[rwo, col]
% 5 landmarks
a= center(1);
b= center(2);
%LM = [a,b];
LM=[];
for i = 0: pi/4: 2*pi-pi/4
    LM = [LM; a+r*cos(i), b+r*sin(i)];
end


% for i = 0: pi/4: 2*pi-pi/4
%     LM = [LM; a+r/2*cos(i), b+r/2*sin(i)];
% end

LM = LM';


