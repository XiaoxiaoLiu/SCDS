function [h_1d, h_scores] = reconHfromP(p,stats,methodType)
%PURPOSE: Reconstruct/predict the image deformation(Hfield) from the shape (P).
%INPUT: p --  points set of the model; stats -- correlation stats;
%       methodType -- MLR/CCA/PLS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



PC_num = stats.PC_num;

switch methodType
    
    
    case {'MLR_XY'} % MLR
        p_scores =  getProjScore( p, stats.P);
        
        p_scores = p_scores(1:PC_num);
        
        h_scores = mapMLR(p_scores, stats.corr);
        
        h_1d = HFieldRecon(h_scores,stats.H);
        
    case{'PLS_XY'}
        
        
        p_scores =  getProjScore( p, stats.P);
        
        p_scores = p_scores(1:PC_num);
        
        h_scores =   mapPLS(p_scores, stats.corr);
        
        
        h_1d = HFieldRecon(h_scores,stats.H);
        
        
    case 'CCA_XY'
        
        %  U = (X - repmat(mean(X),N,1))*A ;
        %  V = (Y - repmat(mean(Y),N,1))*B;
        p_scores =  getProjScore( p, stats.P,PC_num);
      
        
        U = (p_scores - stats.corr.meanX) * stats.corr.A;
      
        V= U.* stats.corr.r;
        
        h_scores = V * inv(stats.corr.B) + stats.corr.meanY;
        h_1d = HFieldRecon(h_scores,stats.H);
        
        
        
        
    case {'MLR_X'}% MLR
        p_scores =  getProjScore( p, stats.P);
        
        p_scores = p_scores(1:PC_num);
        
        h_1d = mapMLR(p_scores, stats.corr);
        
        
    case{'PLS_X'} %PLS
        p_scores =  getProjScore( p, stats.P);
        
        p_scores = p_scores(1:PC_num);
        
        h_1d = mapPLS(p_scores, stats.corr);
        
        
    case {'MLR_Y'}% MLR_Y no good
        
        h_scores = mapMLR(p, stats.corr);
        
        h_1d = HFieldRecon(h_scores,stats.H);
        
        
        
    case 'PLS_Y'%PLS
        
        h_scores =   mapPLS(p, stats.corr);
        
        
        h_1d = HFieldRecon(h_scores,stats.H);
        
        
        
        %     case 'CCA_HDLSS'
        %         U= (p-stats.corr.meanP)*stats.corr.A;
        %
        %         V = U * stats.corr.M ;
        %
        %         h_1d = V/stats.corr.B + stats.corr.meanH; %Y=B\V; V=Y*B
        %         %% error rank deficient
        %
        %         %for testing purposes
        %         h_scores =  (h_1d - stats.H.mean)*stats.H.PCs;
        
    case 'CCA_Y'
        
        %  U = (X - repmat(mean(X),N,1))*A ;
        %  V = (Y - repmat(mean(Y),N,1))*B;
        
        
        U = (p - stats.corr.meanX) * stats.corr.A;
        
        
        V= U.* stats.corr.r;
        
        h_scores = V * inv(stats.corr.B) + stats.corr.meanY;
        h_1d = HFieldRecon(h_scores,stats.H);
        
        
end


