//*******************************************
// This is the Scilab script for Exercise 7.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf();scf(0);
a=gcf(); a.figure_size= [700,350];

// read input data
info = read("header.txt",1,5); // read header information
dtout = info(1); dx = info(2); nx = info(3); nz = info(4); hmax = info(5);
h0in=read("h0.dat",-1,nx); hin=read("h.dat",-1,nx); uin=read("u.dat",-1,nx); 

[ntot nx] = size(hin); ntot = floor(ntot/nz); // total number of frames
x = (0:dx:(nx-1)*dx); // location vector 
h0 = h0in; htot = h0; // bathymetry
time = 0.0; // time counter

for n=1:ntot // animation loop

time = time + dtout; // time in seconds

for i = 1:nz // convert data to individual matrices
  ii = nz*(n-1)+i;
  u(i,1:nx) = uin(ii,1:nx); // horizontal velocity
  h(i,1:nx) = hin(ii,1:nx); // layer thickness
end;

drawlater; clf; 

htot(1:nx) = 0.0; // calculate interface depths
for ii = 1:nz
  i = nz-ii+1; 
  htot(1:nx) = htot(1:nx)+h(i,1:nx);
  ic = i; if ic > 7; ic = ic-7; end; // line colours
  plot2d(x,htot-h0,ic) // draw interface depths
  b = gce(); b.children.thickness = 2;
  c = gca(); c.data_bounds = [0,-100;1000,20];
  c.font_size = 3; 
  c.auto_ticks = ["off","off","on"]; c.sub_ticks = [3,3]; 
  c.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
  c.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0 20], ["-100" "-80" "-60" "-40" "-20" "0" "20"]);  
end;

// draw bathymetry
plot2d(x,-h0,1);
 b = gce(); b.children.thickness = 2;

xstring(400, -95,"Distance (m)");  // x label
xstring(-100, 3,"z (m)"); b = gce(); b.clip_state = "off";
xstring(650, 7,"Time = "+string(int(time/3600))+" hours"); //title

// produce GIF files of each frame (optional)
// if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
// else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//   xs2gif(0,'ex1'+string(n)+'.gif')
//  end
// end
drawnow; xpause(2d4);

end // end of animation

