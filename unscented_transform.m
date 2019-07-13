function [y_mean,Y,P,error]=unscented_transform(u,sigma_x,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(sigma_x,2); % all the sigma points
y_mean=zeros(n,1);
Y=zeros(n,L); 
for k=1:L                   
    Y(:,k)=process(sigma_x(:,k),u,0.1);       
    y_mean=y_mean+Wm(k)*Y(:,k);       
end
error=Y-y_mean(:,ones(1,L));
P=error*diag(Wc)*error'+R; 
