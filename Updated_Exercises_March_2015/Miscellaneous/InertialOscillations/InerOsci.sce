//*******************************************
// Scilab animation of inertial oscillations
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (update)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [600,600];

x=read("output2.txt",-1,3); // read input data

fre = x(1,1); dt = x(1,2); ntot = x(1,3);

for n=2:ntot

drawlater; clf;

time = x(n+1,3); xr = x(n,1); yr = x(n,2);

//draw trajectory
xset('thickness',2)
plot2d(x(2:n,1),x(2:n,2),5,'011');
b = gca(); b.font_size = 3; b.data_bounds = [0,0;30,30]; b.thickness = 1;
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30], ["0" "10" "20" "30"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 10 20 30], ["0" "10" "20" "30"]);

xstring(11, 30.5,"Time = "+string(0.1*int(10*time/(24*3600) ))+" days");
b = gce(); b.font_size =4; b.clip_state = "off";
xstring(13, 0.5,"x (km)");  // add x label
b = gce(); b.font_size =4;
xstring(0.5, 16,"y (km)");  // add z label
b = gce(); b.font_size =4;
xpause(5000.0);

drawnow;

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

end;
