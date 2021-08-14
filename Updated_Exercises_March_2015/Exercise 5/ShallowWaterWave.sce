//*******************************************
// This is the Scilab script for Exercise 5.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf();scf(0);
a=gcf(); a.figure_size= [1000,600];

// read input data
eta1=read("eta.dat",-1,101); 
u1=read("u.dat",-1,101); 
w1=read("w.dat",-1,101);

[ntot nx] = size(eta1); // extract dimension

x = (500:10:1000)'; // location vector 

for n=1:70 // animation loop

time = 2*n; // time in seconds

// grab data line
u = u1(n,51:101)'; // horizontal velocity
w = w1(n,51:101)'; // vertical velocity
eta = eta1(n,51:101)'; // sea-level elevation

drawlater; clf;

subplot(311); // draw graph of eta
plot2d(x,eta,2);
b = gce(); b.children.thickness = 2;
c = gca(); c.data_bounds = [500,-1;1000,1];
c.font_size = 3; 
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [500 600 700 800 900 1000], ["500" "600" "700" "800" "900" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-1 -0.5 0 0.5 1.0], ["-1" "-0.5" "0" "0.5" "1"]);  
xstring(820, 0.8,"Time = "+string(int(time))+" seconds"); //title
b = gce(); b.clip_state = "off";
xstring(510, 0.8,"eta (m)");  // y label
b = gce(); b.clip_state = "off";

subplot(312); // draw graph of u
plot2d(x,u,5); 
b = gce(); b.children.thickness = 2;
c = gca(); c.data_bounds = [500,-1;1000,1];
c.font_size = 3; 
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [500 600 700 800 900 1000], ["500" "600" "700" "800" "900" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-1 -0.5 0 0.5 1.0], ["-1" "-0.5" "0" "0.5" "1"]); 
xstring(510, 0.8,"u (m/s)");  // y label
b = gce(); b.clip_state = "off";

u(1:51) = 0.0;
subplot(313); // draw vertical velocity stickplot
champ(x,0,u,w,arfact = 3.0);
c = gca(); c.data_bounds = [500,-1;1000,1];
c.font_size = 3; 
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [500 600 700 800 900 1000], ["500" "600" "700" "800" "900" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-1 -0.5 0 0.5 1.0], ["-1" "-0.5" "0" "0.5" "1"]); 

xstring(720, -1,"Distance (m)"); b = gce(); b.clip_state = "off";
xstring(510,0.8,"Vertical velocity stickplot"); b = gce(); b.clip_state = "off";

drawnow(); xpause(2d4);

// create GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end

end // end of animation

