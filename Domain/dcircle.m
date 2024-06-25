function d=dcircle(p,varargin)
% Signed distance function for a disc
xc=varargin{1}(1); yc=varargin{2}(1); r=varargin{3}(1);
d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;
