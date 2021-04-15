function E=soft_head_data
%
% The phantom is the same as the standard Shepp-Logan phantom, but with
% softer external skull, so that the gradient of the projection will have
% less variations.
%
% Yoel Shkolnisky, January 2007.
% 

% Ellipses parameters for the soft head phantom:
%
%        A    a     b    x0    y0    phi
%      ---------------------------------
E  = [ 0.3   .69   .92    0     0     0   
       -.2  .6624 .8740   0  -.0184   0
     %  .2  10 0  -.22    10    -1		% random choice
      -.2  .1100 .3100  .22    0    -18	% orig pixel
       -.2  .1600 .4100 -.22    0     18
        .1  .2100 .2500   0    .35    0
        .1  .0460 .0460   0    .1     0
        .1  .0460 .0460   0   -.1     0
        .1  .0460 .0230 -.08  -.605   0 
        .1  .0230 .0230   0   -.606   0
        .1  .0230 .0460  .06  -.605   0   ];
