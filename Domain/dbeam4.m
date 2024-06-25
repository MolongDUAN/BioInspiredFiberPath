function fd=dbeam4(p)
% Signed distance function for a square with hole
r1=70.067164;
xc1=7.5;
yc1=45.567164;
dc1=dcircle(p,xc1,yc1,r1);

r2=112.5;
xc2=7.5;
yc2=70;
dc2=dcircle(p,xc2,yc2,r2);

dr1=drectangle(p,7.5,82.5,-42.5,-24.5);

dr2=drectangle(p,7.5,122.5,-24.5,45.5);

dr3=drectangle(p,-82.5,127.5,42.5,187.5);

dr4=drectangle(p,-112.5,-82.5,2.5,142.5);

d1=min(dc2,dr1);

d2=min(dc1,min(dr2,min(dr3,dr4)));

fd=max(d1,-d2);