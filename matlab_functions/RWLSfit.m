function [beta,fobj,r,rN] = RWLSfit(u,e,w)
%UNTITLED Summary of this funtion goes here
%   Detailed explanation goes here

    if nargin<3,
        error('Error in the number of parameters in function RWLSfit');
    end
    %   Data size
    n=length(u);
    %   Design matrix
    X = [ones(n,1) u];
    Y = e;
    %   Initial weights
    W=spdiags(w,0,n,n);

    %   Optimal estimates
    beta = (X'*W*X)\(X'*W*Y);
      
    r = Y-X*beta;
    %   Once the optimal solution os obtained we can get many statistics
    fobj = r'*W*r;
    
    sigres = sqrt(fobj/(n-2));
    %   Residual variance-covariance matrix
    %   Hat or projection matrix
    P = X*inv(X'*W*X)*X';
    %   Sensitivity matrix
    S= eye(n)-P*W;

     %   Internally studentized residual
    rN = (sqrt(diag(W)).*r)./(sigres*sqrt(1- diag(W).*diag(P))); 
    
end

