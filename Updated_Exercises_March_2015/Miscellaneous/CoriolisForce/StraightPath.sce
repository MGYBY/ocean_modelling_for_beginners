//*******************************************
// Scilab script for visualisation of the
// Coriolis force.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (update)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [600,600];

T = 24.*3600.; // period
fre = -2*%pi/T; // rotation rate
radius = 20; 
time = 0;

dt = T/200; // time step

xp = 0;
yp = 5;
xf = xp;
yf = yp;
uf = 0.0; // speed in m/s
vf = 0.15;

// time loop starts here

fact = 0.9;
fac2 = fact*radius;

xpp(1:200)= %nan;
ypp(1:200)= %nan;

for n = 1:200 // run for 1 period

drawlater; clf;

// isoview scaling 
plot2d(0,0,-1,"030"," ",[-20,-20,20,20])

time = time+dt;

xx = radius*sin(fre*time);
yy = radius*cos(fre*time);

x1 = [fact*xx xx];
y1 = [fact*yy yy];
xset("color",5)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

xx = radius*cos(fre*time);
yy = -radius*sin(fre*time);

x2 = [fact*xx xx];
y2 = [fact*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x2,y2,"lines",1)

xx = -radius*sin(fre*time);
yy = -radius*cos(fre*time);

x3 = [fact*xx xx];
y3 = [fact*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x3,y3,"lines",1)

xx = -radius*cos(fre*time);
yy = radius*sin(fre*time);

x4 = [fact*xx xx];
y4 = [fact*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x4,y4,"lines",1)

xset("color",2)
xfarc(-fac2,fac2,2*fac2,2*fac2,0,360*64)

x1 = [0 0];
y1 = [fact*fac2 fac2];
xset("color",5)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

x1 = [0 0];
y1 = [-fac2 -fact*fac2];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

y1 = [0 0];
x1 = [fact*fac2 fac2];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

y1 = [0 0];
x1 = [-fac2 -fact*fac2];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

//new ball position

xf = xf + dt*uf/1000;
yf = yf + dt*vf/1000;

xp = xf*cos(fre*time)+yf*sin(fre*time);
yp = yf*cos(fre*time)-xf*sin(fre*time);

xpp(n) = xp;
ypp(n) = yp;

xset("color",1)
xfarc(xp-1,yp+1,2.0,2,0,360*64)
xset("color",7)
xfarc(xp-0.8,yp+0.8,1.6,1.6,0,360*64)

xset('thickness',2)
plot2d(xpp,ypp,8,'000','',[-radius,-radius,radius,radius])

title("Rotating Frame of Reference","fontsize",3,"font_style",4);

xpause(5d2);

// creation of GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end
drawnow;

end;
