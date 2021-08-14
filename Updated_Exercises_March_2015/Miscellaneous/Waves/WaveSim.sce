//*****************************************
// Animation of a sinusoidal wave
// made up by vertically moving bars.
//
// Use the SciLab help for more information
// on individual functions used.
//
// Author: J Kaempf, 2015 (update)
//******************************************

clf(); // clears graphic window 
scf(0); // opens new graphic window with ID = 0
a = gcf(); // loads handles of figure
a.figure_size = [600,200];

len = 100.0; // wavelength
per = 60.0; // period
a = 0.4; // amplitude
xrange = 2*len; //range shown in animation
x=[0:xrange/30:xrange]'; // locations
[m1 m2] = size(x); // returns dimension

t = 0.; // start time
trange = 2*per; // simulate two full periods
dt = trange/30.; // time step
ntot = 3*trange/dt; // number of iteration steps

k1 = 1; k2 = m1;// defines the range of locations shown
//k1 = 10; k2 = 10; // this would only display the wave at one location

for n = 1:ntot // start of iteration

drawlater; clf;

plot2d(0,0,1,"010"," ",[0,0.25,1,0.75]);

f = a*sin(2*%pi*(x/len-t/per)); // wave equation

for k = k1:k2; // loop for locations
 if f(k) > 0 // positive values of f
  xset("color",2); // set color to blue
  xfrect(x(k)/220,0.5+0.5*f(k),0.02,0.5*f(k)) // fill rectangle
 else
  xset("color",5); // set color to red
  xfrect(x(k)/220,0.5,0.02,-0.5*f(k)) // fill rectangle
 end;
end;
  drawnow; // pixmap sent to the graphics window 
xpause(1d5); //pause of 100000 micro seconds = 0.1 s

// save current picture frame as GIF file (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

t = t+dt; // time progresses forward

end;
