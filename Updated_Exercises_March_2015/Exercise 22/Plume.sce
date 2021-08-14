//*******************************************
// This is the Scilab script for Exercise 22.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [1000,500];

 // read input data
h0=read("h0.dat",-1,101); e2=read("h2.dat",-1,101);
u1=read("u2.dat",-1,101); v1=read("v2.dat",-1,101);

x1 = (0:2:200)'; y1 = (0:2:100)'; // location vectors  
[ntot nx] =size(u1); ntot = floor(ntot/51); 

for n = 1:ntot // animation loop

time = real(1*n); 

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:101); v2 = v1(jtop:jbot,1:101);
eta2 = e2(jtop:jbot,1:101);

// interpolation of flow field onto coarser grid
u(1:11,1:21) = 0.06; v(1:11,1:21) = 0.06; // scaling 

for j = 1:10; for k = 1:20;
  j1 = 5*(j-1)+1; j2 = j1+5; k1 = 5*(k-1)+1; k2 = k1+5;
  uu = 0.; vv = 0.;
  for jstar = j1:j2; for kstar = k1:k2;
    uu = uu + u2(jstar,kstar); vv = vv + v2(jstar,kstar); 
  end; end;
  u(j,k) = uu/25.; v(j,k) = vv/25.;
end; end;

x = (6:10:206)'; y = (6:10:106)'; // location vectors  

// modify color map
mymap = jetcolormap(64); map = jetcolormap(64);
map = 1-graycolormap(64);
for k = 1:64
 mymap(k,:) =map(64-k+1,:);
end; 
a.color_map = mymap;

drawlater; clf;

Sgrayplot(x1,y1,eta2',zminmax=[0.,50.0]); // plot of relative plume elevation
xset("fpf"," ");
col(1:10) = 40; contour2d(x1,y1,(h0+eta2)',10,col);// // contours of absolute plume elevation

champ(x,y,u',v',1.0); // vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;200,100];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 50 100 150 200], ["0" "50" "100" "150" "200"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 50 100], ["0" "50" "100"]);

xstring(80,101,"time = "+string(int(100*time/24)/100)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(110,-8,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-14,55,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow;

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

end; // end of animation

// To restore the color map in subsequent applications either restart SciLab or use
// clf(gcf(), "reset");
