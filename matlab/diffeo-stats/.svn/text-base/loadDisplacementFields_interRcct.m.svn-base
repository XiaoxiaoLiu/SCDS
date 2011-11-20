function [H] = loadDisplacementFields_interRcct(diffeoDir, outputDir, list, type)
%modified from loadDisplacementFields.m
% to load secondary RCCT data, where hfield of phase50 is not zero 
%-------------------------------------------------------------------------


list_phaseNum = [0,1,2,3,4,5,6,7,8,9];


%dims = readMetaDimSize([diffeoDir,'/phase00/InverseWarpxvec.mhd']);

d = 3;


if nargin <4
    type ='avants';
end

disp('load Hfields');
switch type
    
    case 'avants'  %not completed yet
%         % s = [d,dims];
%         
%         for i = 1: length(list)
%             
%             idx = list(i);
%            
%                 % h = zeros(s);
%                 pad = sprintf('%d0',list_phaseNum(idx));
%                 
%                 diffeoPath =[diffeoDir,'/phase',pad,'/InverseWarp'];
%                 
%                 h = loadDiffeoFieldxyz(diffeoPath);
%                 
%                 
%                 H(i,:) = h(:)';
%          
%             
%         end
        
        %save ( [outputDir,'/hfields.mat'], 'H','-v7.3');
        
        
    case 'fluid'
        for i = 1: length(list)
            
            idx = list(i);
            
        
                % h = zeros(s);
                pad = sprintf('%d0',list_phaseNum(idx));
                
                diffeoPath =[diffeoDir,'/hfield-phase',pad,'.mhd'];
                
                h = loadMETA(diffeoPath);
                s = size(h);
                h = h-eyeHField(s(2:4));
                H(i,:) = h(:)';
                
         
        end
        
        
        
end

%H = H(list,:);