patNo='105';

type='CS';%'GM'
varX=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/xVariance.txt']);
varY=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/yVariance.txt']);
varZ=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/zVariance.txt']);

varP=(varX+varY+varZ)/3;

% varGauss=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/gaussCurvatureVariance.win.txt']);
% 
% 
% varMean=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/meanCurvatureVariance.win.txt']);

varGauss=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/CVariance.txt']);
varMean=readVarianceFile(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/shape/',type,'/SVariance.txt']);


weight = [ 1/sqrt(varP),1/sqrt(varP),1/sqrt(varP), 1/sqrt(varGauss),1/sqrt(varMean)]*sqrt(varP)
