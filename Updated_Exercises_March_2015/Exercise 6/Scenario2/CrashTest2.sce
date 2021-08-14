//*******************************************
// This is the Scilab script for Exercise 6
// (Scenario 2).
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************
clf; scf(0);
a=gcf(); a.figure_size= [700,350];

// read input data
h0=read("h0.dat",-1,101); 
eta1=read("eta.dat",-1,101);
u1=read("u.dat",-1,101); 

[ntot nx] = size(eta1); // extract dimension
x = (1:10:1010)'; // location vector 
zero = x;
zero(1:101) = 0.;

for n=1:ntot // start of animation loop

time = 2*n; // time in seconds

// grab data line
u = u1(n,1:101)'; // horizontal velocity
eta = eta1(n,1:101)'; // sea-level elevation

drawlater; clf;
  
plot2d(x,eta,5);  
b = gce(); b.children.thickness = 2;
c = gca(); c.data_bounds = [0,10;1000,20];
c.font_size = 3; 
c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
c.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
c.y_ticks = tlist(["ticks", "locations","labels"],..
 [10 12.5 15 17.5 20], ["10" "12.5" "15" "17.5" "20"]);  

plot2d(x(1:2:100),zero(1:2:100),0,'000');
b = gce(); b.children.thickness = 2;
plot2d(x,-h0,2,'000');
b = gce(); b.children.thickness = 2;

xstring(650, 20,"Time = "+string(int(time))+" seconds"); 
xstring(420, 10.4,"Distance (m)");  
xstring(20, 15,"eta (m)");  

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

drawnow; xpause(2d4);

end // end of animation loop

