//*******************************************
// This is the Scilab script for Exercise 10.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [700,700];

// read input data
eta=read("eta.dat",-1,51); 
h0=read("h0.dat",-1,51); h1=read("h.dat",-1,51);
u1=read("u.dat",-1,51); v1=read("v.dat",-1,51); 

x = (0:0.1:5)'; y = (0:0.1:5)'; // location vectors  
[ntot nx] =size(u1);
ntot = floor(ntot/nx); time = 0.0;

for n = ntot // animation loop (only shows last output)

time = n*2/24; 

// grab respective data blocks
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:51); 
v2 = v1(jtop:jbot,1:51); 
h2 = h1(jtop:jbot,1:51);

uh = u2;
vh = v2;

//calculate transports

for j = 1:51; for k = 1:50;
  uh(j,k) = u2(j,k)*0.5*(h2(j,k)+h2(j,k+1));
end; end;

uh(1:51,51) = 0.0;

for j = 1:50; for k = 1:51;
  vh(j,k) = v2(j,k)*0.5*(h2(j,k)+h2(j+1,k));
end; end;

vh(51,1:51) = 0.0;
  
u2 = uh; v2 = vh;

drawlater;clf;

col(1:11) = 1; 
xset("fpf"," "); xset("thickness",2); 
contour2d(x,y,h0',0:2:20,col); // bathymetry contours
b = gca(); b.font_size = 3; b.data_bounds = [0,0;5,5];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);

uu = u2(1:3:51,1:3:51);
vv = v2(1:3:51,1:3:51);
champ(x(1:3:51),y(1:3:51),uu',vv',0.7)// vector plot
xstring(3.0,5.05,"time = "+string(int(time))+" days"); //add time

xstring(2.2,-0.3,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.6,2.5,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off";b.font_size = 3;

drawnow();

// save frames as GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
// end

end; // end of animation loop
