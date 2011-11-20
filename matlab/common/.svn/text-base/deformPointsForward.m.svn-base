function dpts = deformPointsForward(pts, h, extrapVal,method)
% returns dpts = h \circ pts

if nargin < 4
  method = 'linear';
end
if nargin < 3
  extrapVal = 0;
end

switch size(pts,1)
  case 3,
    dpts(1,:) = interp3(squeeze(h(1,:,:,:)),pts(2,:),pts(1,:),pts(3,:),...
      method);
    dpts(2,:) = interp3(squeeze(h(2,:,:,:)),pts(2,:),pts(1,:),pts(3,:),...
      method);
    dpts(3,:) = interp3(squeeze(h(3,:,:,:)),pts(2,:),pts(1,:),pts(3,:),...
      method);
  case 2,
    dpts(1,:) = interp2(squeeze(h(1,:,:)),pts(2,:),pts(1,:),method);
    dpts(2,:) = interp2(squeeze(h(2,:,:)),pts(2,:),pts(1,:),method);
  otherwise,
    error('fields must be 2 or 3 dimensional');
end
