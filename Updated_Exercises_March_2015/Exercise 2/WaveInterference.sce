//***************************************************************************
// This Scilab script simulated the interference of two sinusoidal waves.
//
// Consult the SciLab help for more information on individual functions used.
// Author: J. Kaempf, 2015 (updated version)
//****************************************************************************

clf; // clears graphic window 
scf(0); //initiate new graphic window

len1 = 100.0; // wavelength of wave 1
per1 = 60.0; // period of wave 1
len2 = 95.0; // wavelength of wave 2
per2 = -30.0; // period of wave 2
amp = 1.0; // amplitudes of both waves
xrange = 10*len1; //x-range shown in animation
x=[0:xrange/500:xrange]'; // discrete grid points in x-direction

t = 0.; // start time
trange = 10*per1; // simulate ten periods of wave 1
dt = trange/200.; // time step
ntot=trange/dt; // number of iteration steps

for n = 1:ntot // start of iteration

f1 = amp*sin(2*%pi*(x/len1-t/per1)); // wave 1 equation
f2 = amp*sin(2*%pi*(x/len2-t/per2)); // wave 2 equation

f3 = f1+f2; // superposition of waves

drawlater; clf;  // store graphics for later display

subplot(311)

plot2d(x,f1,2);
b = gce(); // handle to curve properties
b.children.thickness = 2; // increase line thickness
c = gca(); // handle to axis properties
c.data_bounds = [0,-2;1000,2];
c.font_size = 3; // increase font size of axis labels
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; // creates subticks
// the following command creates axis labels at specified positions
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-2 -1 0 1 2], ["-2" "-1" "0" "1" "2"]); 
title("Wave 1","fontsize",3) // add title

subplot(312)
plot2d(x,f2,5);

b = gce(); b.children.thickness = 2; 
c = gca(); c.font_size = 3; c.data_bounds = [0,-2;1000,2];  
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-2 -1 0 1 2], ["-2" "-1" "0" "1" "2"]); 
title("Wave 2","fontsize",3)

subplot(313)
plot2d(x,f3,6)
b = gce(); b.children.thickness = 2; 
c = gca(); c.font_size = 3; c.data_bounds = [0,-2;1000,2]; 
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-2 -1 0 1 2], ["-2" "-1" "0" "1" "2"]); 
title("Wave 1+2","fontsize",3)
drawnow(); // draw all graphs
xpause(10000); // delay next steep by number of microseconds

t = t+dt; // time progresses forward

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

end; // end of iteration loop
