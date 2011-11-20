function [H] = loadDisplacementFields(diffeoDir, outputDir, list, type,prefix)
%PURPOSE: Load displacement fields calculated from deformable registration.
%OUPUT: H -- displacement fields data matrix: N *(3*K), where N is the
%       number of training sample, K is the number of voxels in the image.
%-------------------------------------------------------------------------


list_phaseNum = [0,1,2,3,4,5,6,7,8,9];

if nargin<4
    type='fluid';
end
if nargin <5
    switch type
        case 'avants'
            prefix ='InverseWarp';
        case 'fluid'
            prefix='hfield-cinePhase';
    end
    
end

disp('load Hfields');
switch type
    
    case 'avants'
        
        for i = 1: length(list)
            
            idx = list(i);
            if (list_phaseNum(idx)==5)
                % zeros
                H(i,:) = zeros(size(H(1,:)));
            else
                
                pad = sprintf('%d0',list_phaseNum(idx));
                
                diffeoPath =[diffeoDir,'/phase',pad,'/',prefix];
                
                h = loadDiffeoFieldxyz(diffeoPath);
                
                
                H(i,:) = h(:)';
                
            end
            
        end
        
        %save ( [outputDir,'/hfields.mat'], 'H','-v7.3');
        
        
    case 'fluid'
        
        for i = 1: length(list)
            
            idx = list(i);
            
            if (list_phaseNum(idx)==5)
                % atlas Phase: zeros
                H(i,:) = zeros(size(H(1,:)));
            else
                
                pad = sprintf('%d0',list_phaseNum(idx));
                
                diffeoPath =[diffeoDir,'/',prefix,pad,'.mhd'];
                
                h = loadMETA(diffeoPath);
                s = size(h);
                h = h-eyeHField(s(2:4));
                
                H(i,:) = h(:)';
                
            end
        end
        
        
        
end

%H = H(list,:);