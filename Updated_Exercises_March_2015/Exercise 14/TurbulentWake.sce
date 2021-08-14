//*******************************************
// This is the Scilab script for Exercise 14.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf(); scf(0); a=gcf(); a.figure_size= [1000,500];

// read input data
tin=read("t.dat",-1,101); h0=read("h0.dat",-1,101);
u1=read("u.dat",-1,101); v1=read("v.dat",-1,101); 

x = (0:0.1:10)'; y = (0:0.1:5)'; // location vectors  
[ntot nx] =size(u1);
ntot = floor(ntot/51); time = 0.0;

for n = 1:ntot // animation loop

time = real(n)*900/3600.;

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:101); v2 = v1(jtop:jbot,1:101); 
t = tin(jtop:jbot,1:101); 

drawlater;clf();

xset("fpf"," "); //suppress label output 
col1(1:10) = 1; contour2d(x,y,h0',10,col1); // bathymetry contours
b = gce(); b.thickness = 2;
xset("fpf"," "); 
col2 = 1:25; contour2d(x,y,t',25,col2);// tracer contours
//xset("fpf"," "); //suppress label output 

b = gca(); b.font_size = 3; b.data_bounds = [0,0;10,5];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 2 4 6 8 10], ["0" "2" "4" "6" "8" "10"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);

uu = u2(1:3:51,1:3:101);
vv = v2(1:3:51,1:3:101);
xset("thickness",1); 
champ(x(1:3:101),y(1:3:51),uu',vv',1.0,rect=[0,0,10,5]); // vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;10,5];
xstring(4,5.02,"time = "+string(int(time))+" hrs"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(4.8,-0.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.8,2.3,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
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
