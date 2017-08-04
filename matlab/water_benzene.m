% start fresh
clear
close all
clc

% inputs: dimensions of the system (approximate)
L = 22.0;

% nsolutex * nsolutey * nsoutez has to be even!
nsolutex = 2;
nsolutey = 4;
nsolutez = 4;
 
ntrials = 1e7; % for water packing

% WATER

% monomer TIP4-Ew
roh = 0.9572;
theta = 104.52*(pi/180);

% preliminary calculations for filling the box
doo = 2.6; % minimum oxygen-oxygen separation
dex = 2; % excluded distance between graphene and water
slack = 2;
nwatercent = round(1/1e21/1.66e-24/18*((L-slack)^3/1e3));
nwaterside = round(nwatercent/2);

coordsOxy = xyz_packing(L,L,L,nwatercent,doo,ntrials);

nwater = length(coordsOxy);

% coordinates = [index molid type charge x y z]
mon_xyz = [0  0  0; roh 0 0; roh*cos(theta) roh*sin(theta) 0];
xyz_water = zeros(3*nwater,3);
for i=1:nwater            
    rot_angle = pi*rand;
    
    Rx = [1 0 0; 0 cos(rot_angle) -sin(rot_angle); 0 sin(rot_angle) cos(rot_angle)];
    Ry = [cos(rot_angle) 0 sin(rot_angle); 0 1 0; -sin(rot_angle) 0 cos(rot_angle)];
    Rz = [cos(rot_angle) -sin(rot_angle) 0; sin(rot_angle) cos(rot_angle) 0; 0 0 1];
    
    Rtot = Rx*Ry*Rz;
    
    xyz_water(3*i-2:3*i,:) = mon_xyz*Rtot+repmat(coordsOxy(i,:),3,1);
end

% bonds = [index_b type ID1 ID2]
mon_bonds = [1 2; 1 3];
nbondsw = 2*nwater;
bondsw = zeros(nbondsw,2);
for i=1:nwater
    bondsw(2*i-1:2*i,:) = mon_bonds + 3*(i-1);
end

% angles = [index_a type ID1 ID2 ID3]
mon_angles = [2 1 3];
nanglesw = nwater;
anglesw = zeros(nanglesw,3);
for i=1:nwater
    anglesw(i,:) = mon_angles + 3*(i-1);
end

% molecule index
molid_water = zeros(3*nwater,1);
for i=1:nwater
    molid_water(3*i-2:3*i) = zeros(3,1)+i;
end

% atom types
type_water = repmat([1;2;2],nwater,1);

%%%%% benzene %%%%%%%%%%%%%%%%

xyz_benz_single = [0.000  1.396  0.000
           1.209  0.698  0.000
           1.209 -0.698  0.000
           0.000 -1.396  0.000
          -1.209 -0.698  0.000
          -1.209  0.698  0.000
           0.000  2.479  0.000
           2.147  1.240  0.000
           2.147 -1.240  0.000
           0.000 -2.479  0.000
          -2.147 -1.240  0.000
          -2.147  1.240  0.000];

rot_angle = pi/2;
Ry = [cos(rot_angle) 0 sin(rot_angle); 0 1 0; -sin(rot_angle) 0 cos(rot_angle)];
xyz_benz_single = xyz_benz_single*Ry;      
natbenz = length(xyz_benz_single);
      
bondsb_single = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1; 1 7; 2 8; 3 9; 4 10; 5 11; 6 12];
anglesb_single = [1 2 3; 2 3 4; 3 4 5; 4 5 6; 5 6 1; 6 1 2; 7 1 6; 7 1 2; 8 2 1; 8 2 3; 9 3 2; 9 3 4; 10 4 3; 10 4 5; 11 5 4; 11 5 6; 12 6 5; 12 6 1];
dihedralsb_single = [1 2 3 4; 2 3 4 5; 3 4 5 6; 4 5 6 1; 5 6 1 2; 6 1 2 3];
impropersb_single = [1 6 2 7; 2 1 3 8; 3 2 4 9; 4 3 5 10; 5 4 6 11; 6 5 1 12];

% many benzene
incrposx = L/nsolutex;
incrposy = L/nsolutey;
incrposz = L/nsolutez;

xyz_benz = [];

xyz_benz_single = xyz_benz_single + [zeros(natbenz,1) zeros(natbenz,1)+incrposy/2 zeros(natbenz,1)+incrposz/2];

count = 1;
for ii=1:nsolutex
    for jj=1:nsolutey
        for kk=1:nsolutez
            xyz_benz =[ xyz_benz; xyz_benz_single + [zeros(natbenz,1)+(ii-1)*incrposx zeros(natbenz,1)+(jj-1)*incrposy zeros(natbenz,1)+(kk-1)*incrposz]];
            count = count + 1;
        end
    end
end

xyz_benz = [xyz_benz(:,1)-mean(xyz_benz(:,1)) xyz_benz(:,2)-mean(xyz_benz(:,2)) xyz_benz(:,3)-mean(xyz_benz(:,3))]+L/2;


nsolutetotal = nsolutex*nsolutey*nsolutez;

bondsb = zeros(length(bondsb_single)*nsolutetotal,2);
anglesb = zeros(length(anglesb_single)*nsolutetotal,3);
dihedralsb = zeros(length(dihedralsb_single)*nsolutetotal,4);
impropersb = zeros(length(impropersb_single)*nsolutetotal,4);
for ii=1:nsolutetotal
    bondsb((ii-1)*length(bondsb_single)+1:ii*length(bondsb_single),:) = bondsb_single + (ii-1)*natbenz;
    anglesb((ii-1)*length(anglesb_single)+1:ii*length(anglesb_single),:) = anglesb_single + (ii-1)*natbenz;
    dihedralsb((ii-1)*length(dihedralsb_single)+1:ii*length(dihedralsb_single),:) = dihedralsb_single + (ii-1)*natbenz;
    impropersb((ii-1)*length(impropersb_single)+1:ii*length(impropersb_single),:) = impropersb_single + (ii-1)*natbenz;
end

nbondsb = length(bondsb);
nanglesb = length(anglesb);
ndihedralsb = length(dihedralsb);
nimpropersb = length(impropersb);

prov = repmat((1:nsolutetotal),12,1);
molid_benz = prov(:) + molid_water(end);
type_benz = repmat([zeros(natbenz/2,1)+3;zeros(natbenz/2,1)+4],nsolutetotal,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put together coordinates
xyz_total = [xyz_water; xyz_benz];
type_total = [type_water; type_benz];
molid_total = [molid_water; molid_benz];
index_total = (1:length(xyz_total))';
charge_total = zeros(length(xyz_total),1);
coordinates = [index_total molid_total type_total charge_total xyz_total];


% put together bonds, angles, and dihedrals
nbonds = nbondsw + nbondsb;
index_b = (1:nbonds)';
type_b = [zeros(nbondsw,1)+1; repmat([zeros(6,1)+2; zeros(6,1)+3],nsolutetotal,1)];
bondsAt = [bondsw; bondsb + length(xyz_water)];
bonds = [index_b type_b bondsAt];

nangles = nanglesw + nanglesb;
index_a = (1:nangles)';
type_a = [zeros(nanglesw,1)+1; repmat([zeros(6,1)+2; zeros(12,1)+3],nsolutetotal,1)];
anglesAt = [anglesw; anglesb + length(xyz_water)];
angles = [index_a type_a anglesAt];

type_d = repmat(zeros(6,1)+1,nsolutetotal,1);
dihedrals = [(1:ndihedralsb)' type_d dihedralsb+length(xyz_water)];

type_i = repmat(zeros(6,1)+1,nsolutetotal,1);
impropers = [(1:nimpropersb)' type_i impropersb+length(xyz_water)];

% inputs are: coordinates, bonds, angles, dihedrals, impropers
% Write LAMMPS data file open the file with write permission
% coordinates = [index molid type charge x y z]
% bonds = [index_b type ID1 ID2]
% angles = [index_a type ID1 ID2 ID3]

name = strcat('L',num2str(round(L)),'nwat',num2str(nwater),'nbenz',num2str(nsolutetotal),'.data');
fid = fopen(name, 'w');

fprintf(fid,'%s\n\n','water box');

natoms = length(coordinates);
nbonds = length(bonds);
nangles = length(angles);
ndihedrals = length(dihedrals);
nimpropers = length(impropers);
fprintf(fid,'%d atoms\n',natoms);
fprintf(fid,'%d bonds\n',nbonds);
fprintf(fid,'%d angles\n',nangles);
fprintf(fid,'%d dihedrals\n',ndihedrals);
fprintf(fid,'%d impropers\n\n',nimpropers);

natomtypes = 4;
nbondtypes = 3;
nangletypes = 3;
ndihedraltypes = 1;
nimpropertypes = 1;
fprintf(fid,'%d atom types\n',natomtypes);
fprintf(fid,'%d bond types\n',nbondtypes);
fprintf(fid,'%d angle types\n',nangletypes);
fprintf(fid,'%d dihedral types\n',ndihedraltypes);
fprintf(fid,'%d improper types\n\n',nimpropertypes);

% empty box padding in the x-dimension
minx = 0; 
miny = 0;
minz = 0;
maxx = L;
maxy = L;
maxz = L;

fprintf(fid,'%.2f %.2f xlo xhi\n',minx,maxx);
fprintf(fid,'%.2f %.2f ylo yhi\n',miny,maxy);
fprintf(fid,'%.2f %.2f zlo zhi\n\n',minz,maxz);

% O=1, H=2, C=3
fprintf(fid,'%s\n\n','Masses');
mass = [15.9994; 1.008; 12; 1.008];
id = [1; 2; 3; 4];
masses = [id mass];
for i=1:length(id)
    prov = masses(i,:);
    fprintf(fid,'%d %.2f\n',prov);
end

%print atoms
fprintf(fid,'\n%s\n\n','Atoms');
for i=1:length(coordinates)
    prov = coordinates(i,:);
    fprintf(fid,'%d %d %d %.2f %.2f %.2f %.2f\n',prov);
end

%print bonds
fprintf(fid,'\n%s\n\n','Bonds');
for i=1:length(bonds)
    prov = bonds(i,:);
    fprintf(fid,'%d %d %d %d\n',prov);
end

%print angles
fprintf(fid,'\n%s\n\n','Angles');
for i=1:length(angles)
    prov = angles(i,:);
    fprintf(fid,'%d %d %d %d %d\n',prov);
end

%print dihedrals
fprintf(fid,'\n%s\n\n','Dihedrals');
for i=1:length(dihedrals)
    prov = dihedrals(i,:);
    fprintf(fid,'%d %d %d %d %d %d\n',prov);
end

%print impropers
fprintf(fid,'\n%s\n\n','Impropers');
for i=1:length(impropers)
    prov = impropers(i,:);
    fprintf(fid,'%d %d %d %d %d %d\n',prov);
end

fclose(fid);
