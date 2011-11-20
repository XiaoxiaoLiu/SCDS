function [A,B, r,U,V] = cca_HDLSS(X,Y)
% input : X: d1 *N, with N observations/samples, d1 dimenstional data
%         Y: d2 *N
% output: U: CC score    N*L
%         V: CC score    N*L
%         alpha: CC directions   d1 * L (first L eigenvectors)
%         beta: CC directions   d2 * L (L<N)
%         r: correlation
% Xiaoxiao Liu, 2009 Sept




d1 = size(X,1);
d2 = size(Y,1);



N = size(X,2);
X = X-repmat(mean(X')',1,N);
Y = Y-repmat(mean(Y')',1,N);


%% method 1


alpha = X'*X;
beta =double( Y'*Y);

% % remove linear dependent data
% % alpha = alpha- repmat(mean(alpha,1), N, 1);
% % beta = beta- repmat(mean(beta,1), N, 1);

% [U,S,V]=svd(alpha);
% idx=find(diag(S) > eps(abs(S(1)))*size(alpha,1));
% alpha1=U(idx,idx)*S(idx,idx)*V(idx,idx)';
% % 
% [U,S,V]=svd(beta);
% idx=find(diag(S) > eps(abs(S(1)))*size(beta,1));
% beta1=U(idx,idx)*S(idx,idx)*V(idx,idx)';

 addpath('/home/sharonxx/matlab/utilites/LSCCA/package/utilities/');
[a,b,r] = CCA(alpha,beta);
a = normc(a);
b = normc(b);
u = (alpha-repmat(mean(alpha),size(alpha,1),1)) *a;
v = (beta-repmat(mean(beta),size(beta,1),1)) *b;
% [a,b,r,u,v] = canoncorr(alpha,beta);


% idx=find(sum(a,2)==0)
% %delete zero row
% a=a([1:idx-1,idx+1:end],:);
% b=b([1:idx-1,idx+1:end],:);


A = X*a;
B = Y*b;
U = X'*A;
V = Y'*B;


% %% method 2
% 
% 
% % HDLSS
% alpha = X'*X+ 10^(-12)*eye(N);
% beta = Y'*Y+ 10^(-12)*eye(N);
% 
% 
% opts.issym = 1;
% [alpha,lamda] = eigs( double(inv(alpha*alpha)*(alpha*beta)*inv(beta*beta)* (beta'*alpha')),N,'lr',opts); % alp
% lamda=diag(lamda);
% 
% L=size(alpha,2);
% % --- select first severl real eigen values ---
% for i=1:L
%     if (~isreal(lamda(i)))
%         L=i-1;
%         break;
%     end
% end
% 
% if L<2
%     error('somthing wrong with the eigen decomposition in CCalpha!')
% end
% 
% 
% 
% b = inv(beta*beta)* (beta'*alpha')*a;     % beta: N*L
% b = normc(b);
% 
% %% calculate the correlations
% for i = 1:L
%     
%     A(:,i) = X*a; B(:,i) = Y*b;
%     r(i) = a'*alpha*beta*b/ (sqrt(a'*alpha*alpha*a)*sqrt(b'*beta*beta*b));
% end
% 
% 
% 
% 
