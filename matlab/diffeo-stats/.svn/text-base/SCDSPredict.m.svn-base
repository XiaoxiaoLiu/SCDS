function [h, h_score] = SCDSPredict(p, stats,PC_num)
%Reconstruct the h from the prediction from P
%h:hfield in 1D
%stats is gnerated from CCA 

if nargin<3
    PC_num=1;
end


%CCA
p_scores =  getProjScore( p, stats.P,PC_num);

U = (p_scores - stats.corr.meanX) * stats.corr.A;

V= U.* stats.corr.r;

h_score = V /stats.corr.B + stats.corr.meanY;

h= HFieldRecon(h_score,stats.H);


% %% MLR
%   p_scores =  getProjScore( p, stats.P);
%         
%         p_scores = p_scores(1:PC_num);
%         
%         h_score = mapMLR(p_scores, stats.corr);
%         
%         h = HFieldRecon(h_score,stats.H);