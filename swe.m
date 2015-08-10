clear all;
Mc =96*4;        % number of cells in x
Nc =32*4;        % number of cells in y
ntimes =1000000;       % number of time steps

xmin =-24.0; xmax=24.0;
ymin = -8.0; ymax= 8.0;

dx = (xmax-xmin)/Mc;
dy = (ymax-ymin)/Nc;

xh = xmin+dx/2:dx:xmax-dx/2;
yh = ymin+dy/2:dy:ymax-dy/2;
xe = xmin:dx:xmax;
ye = ymin:dy:ymax;

%%%%%%%%%Open the FORTRAN binary files%%%%%%%%%%%%
fu  = fopen('./output_3/fort.11','r','b');
fv  = fopen('./output_3/fort.12','r','b');
fzt = fopen('./output_3/fort.13','r','b');

%%%%%%%%%%%%Contour levels for pressure%%%%%%%%%%
ztlevs = [-0.035 -0.025 -0.015 -0.005
           0.005  0.015  0.025  0.035
           0.045  0.055  0.065  0.075
           0.085  0.095  0.105  0.115
           0.125  0.135  0.145  0.155
           0.165  0.175  0.185  0.195
           0.205  0.215  0.225  0.235];

nx = 4; ny=4; scale=2;  % draw every fourth vector and scale it
for it = 1:ntimes
  it
  u = fortread( fu, Mc+1,Nc);
  v = fortread( fv, Mc  ,Nc+1);
 zt = fortread(fzt, Mc  ,Nc);
% Plot pressure field as filled contours
  hold off;
  contour(xh,yh,zt',ztlevs);
% imagesc(xh,yh,zt');
  colorbar
% Plot velocity vectors as arrows but first interpolate them to cell centers
  uh = 0.5*(u(1:Mc,1:Nc)+u(2:Mc+1,1:Nc));
  vh = 0.5*(v(1:Mc,1:Nc)+v(1:Mc,2:Nc+1));
  hold on;
  quiver(xh(1:nx:Mc),yh(1:ny:Nc), uh(1:nx:Mc,1:ny:Nc)',...
                                  vh(1:nx:Mc,1:ny:Nc)',...
                                  scale,'Color','k');
% plot(xe,u(:,16),'r',ye,v(48,:),'b')
  axis([xmin xmax ymin ymax]);
  axis image
  pause
end
