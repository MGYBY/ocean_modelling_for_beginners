//*******************************************
// This is the Scilab script for Exercise 20.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

// This script produces a graph showing float trajectories in the bottom layer
// together with the final flow field.

clf; scf(0); a=gcf(); a.figure_size= [600,600];

// read input data
hin=read("h1.dat",-1,51);
trx1=read("TRx2.dat",-1,600); try1=read("TRy2.dat",-1,600); 
uin=read("u2.dat",-1,51); vin=read("v2.dat",-1,51);

x = (0:1:50)'; y = (0:1:50)'; // location vectors  
[ntot nx] =size(hin); ntot = floor(ntot/nx); time = ntot*2/24; n = ntot;

// grab respective data block 
jtop = (n-1)*51+1; jbot = jtop+50; 
u = uin(jtop:jbot,1:51); v = vin(jtop:jbot,1:51);

// average velocity components at scalar grid point
uh = u; vh =v;
for j = 1:51; for k = 1:50; uh(j,k) = 0.5*(u(j,k)+u(j,k+1));end;end;
for j = 1:50; for k = 1:51; vh(j,k) = 0.5*(v(j,k)+v(j+1,k));end;end;

// remove velocities below a certain threshold value

speed = uh;
for j = 1:51; for k = 1:51;
speed(j,k) = sqrt(uh(j,k)*uh(j,k)+vh(j,k)*vh(j,k));
if speed(j,k) < 0.05 ; uh(j,k) = 0; vh(j,k) = 0; end;end;end;

drawlater; clf();

ii = 1; 
for i = 1:200 
ii = ii+1; // change of color
trax = trx1(:,i); // grab respective data line
tray = try1(:,i); // grab respective data line
if ii > 4; ii = ii-3; end;
plot2d(trax,tray,ii);//,rect=[0,0,50,50],nax=[4,6,4,6]);
end;
ii = 4; 
for i = 201:400
ii = ii+1;
trax = trx1(:,i); // grab respective data line
tray = try1(:,i); // grab respective data line

if ii > 7; ii = ii-3; end;
plot2d(trax,tray,ii);//,rect=[0,0,50,50],nax=[4,6,4,6]);
end;


uu = uh(1:2:51,1:2:51);
vv = uh(1:2:51,1:2:51);

// averaging
ua(1:26,1:26) = 0.0; va(1:26,1:26) = 0.0; 
for j = 1:25; for k = 1:25;
  j1 = 2*j-1; j2 = j1+1; k1 = 2*k-1; k2 = k1+1;
  uu = 0.; vv = 0.;
  for jstar = j1:j2; for kstar = k1:k2;
    uu = uu + uh(jstar,kstar); vv = vv + vh(jstar,kstar); 
  end; end;
  ua(j,k) = uu/4.; va(j,k) = vv/4.;
end; end;

xa = (1:2:51)'; ya = (1:2:51'); // location vectors  

champ(xa,ya,10*ua',10*va',1);// velocity vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;50,50];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);

xstring(30.0,50.5,"time = "+string(int(10*time)/10)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(22.5,-6.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-7.7,23.5,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow;

