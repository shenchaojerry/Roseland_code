function projections=shepp_logan_proj_2(theta,t)
%
% function projections=shepp_logan_proj_2(theta,t)
%
% Compute projections of the Shepp Logan phantom that correspond to angle
% theta, with samples at distance t from the origin.
%
% Each row corresponds to a fixed t. Each column corresponds to a fixed
% theta. 
% theta should be in the range [0,2*pi].
%
% Example:
%
%   theta=0:0.01:2*pi;  
%   t=-1.5:0.01:1.5;
%   P=shepp_logan_proj_2(theta,t);
%   I=iradon(P,theta*180/pi);
%   imagesc(I);
%   colormap(gray);
%
% This example computes the projections of the Shepp-Logan phanom for the
% given theta and t and then reconstructs the phantom using iradon.
%
% Yoel Shkolnisky, January 2007

n=length(t);
m=length(theta);

projections=zeros(n,m);

ellipses=create_ellipses;
ne=length(ellipses);


for i=1:m    
    for k=1:ne
        p=project_ellipse(ellipses{k},theta(i),t);
        projections(:,i)=projections(:,i)+p(:);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliry functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=project_ellipse(ellipse_params,theta,t)
%
% Compute the projection of the ellipse given by ellipse_params for rays
% with direction theta and distance from the origiin t.
%
% Input parameters:
%
%    ellipse_params - Struct with the ellipse parameters (
%    theta - Projection angle. The projection is perpendicular to the
%            vector of angle theta with the x-axis
%    t - Vector of sampling points along the direction theta.
%
% This structure follows the description in Kak pp. 52-56.
%
% Yoel Shkolnisky, October 2006

c=ellipse_params.center;
cX=c(1); cY=c(2);

alpha=ellipse_params.rotation;

s=sqrt(cX*cX+cY*cY);
if s>0
    gamma=atan2(cY,cX);
else
    gamma=0;
end

theta1=theta-alpha;
t1=t-s*cos(gamma-theta);

A=ellipse_params.major_axis;
B=ellipse_params.minor_axis;
rho=ellipse_params.refractive_index;
aa=A*A*cos(theta1)*cos(theta1)+B*B*sin(theta1)*sin(theta1);

I=find(abs(t1)<=sqrt(aa));
p=zeros(size(t1));

if ~isempty(I)
    p(I)=2.*rho.*A.*B.*sqrt(aa-t1(I).*t1(I))/(aa);
end



function ellipses=create_ellipses

shep = soft_head_data;
n=size(shep,1);
ellipses=cell(n,1);

for k=1:n
    ellipses{k}.refractive_index=shep(k,1);
    ellipses{k}.major_axis=shep(k,2);
    ellipses{k}.minor_axis=shep(k,3);
    ellipses{k}.center=[shep(k,4) shep(k,5)];
    ellipses{k}.rotation=shep(k,6)*pi/180;
    
end
    
