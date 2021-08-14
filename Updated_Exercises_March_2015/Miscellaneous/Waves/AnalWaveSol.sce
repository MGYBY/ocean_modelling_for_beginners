//*******************************************
// Scilab script for visualisation of the 
// dynamics of long surface gravity waves.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (update)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [800,400];

len = 500.0; // wavelength of wave
eta0 = 1.0; // amplitude of wave
g = 9.81; // acceleration due to gravity
h = 20.0; // water depth
c = sqrt(g*h); // phase speed
per = len/c; // period of wave
u0 = eta0*sqrt(g/h); // u amplitude

xrange = 2*len; //x-range shown in animation
x=[0:xrange/20:xrange]'; // discrete grid points in x-direction
t = 0.; // start time
trange = 2*per; // simulate 2 wave periods
dt = trange/100.; // time step
ntot=trange/dt; // number of iteration steps

// initial locations of fluid parcels
xpos1 = x; zpos1(1:21) = 1.0; xpos2 = x; zpos2(1:21) = 6.0; 
xpos3 = x; zpos3(1:21) = 11.0; xpos4 = x; zpos4(1:21) = 16.0; 

for n = 1:ntot // start of iteration

drawlater; clf();

eta = eta0*sin(2*%pi*(x/len-t/per)); // solution for eta 
u = u0*sin(2*%pi*(x/len-t/per));   // solution for u
dwdz = -2*%pi*u0/len*cos(2*%pi*(x/len-t/per)); // vertical gradient of w

// new locations
xpos1 = xpos1+dt*u; w = dwdz.*zpos1; zpos1 = zpos1+dt*w;
xpos2 = xpos2+dt*u; w = dwdz.*zpos2; zpos2 = zpos2+dt*w;
xpos3 = xpos3+dt*u; w = dwdz.*zpos3; zpos3 = zpos3+dt*w;
xpos4 = xpos4+dt*u; w = dwdz.*zpos4; zpos4 = zpos4+dt*w;

// draw graphs
xset("thickness",2)
plot2d(xpos1,-h+zpos1,-9);
plot2d(xpos2,-h+zpos2,-9);
plot2d(xpos3,-h+zpos3,-9);
plot2d(xpos4,-h+zpos4,-9);
plot2d(x,eta,2,'000');
b = gca();  b.font_size = 3; b.data_bounds = [0,-20;1000,2];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [-20 -15 -10 -5 0], ["-20" "-15" "-10" "-5" "0"]);

xset("thickness",1); xset("font size",3);

drawnow; xpause(2d4);
t = t+dt; // time progresses forward

//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end; // reference point for iteration loop
