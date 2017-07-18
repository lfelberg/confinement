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

% ions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
incrposx = L/nsolutex;
incrposy = L/nsolutey;
incrposz = L/nsolutez;

ntotalion = nsolutex*nsolutey*nsolutez;

count = 1;
xyz_ion = zeros(ntotalion,3);
for ii=1:nsolutex
    for jj=1:nsolutey
        for kk=1:nsolutez
            xyz_ion(count,:) = [(ii-1)*incrposx (jj-1)*incrposy (kk-1)*incrposz];
            count = count + 1;
        end
    end
end

xyz_ion = [xyz_ion(:,1)-mean(xyz_ion(:,1)) xyz_ion(:,2)-mean(xyz_ion(:,2)) xyz_ion(:,3)-mean(xyz_ion(:,3))]+L/2;

type_ion = repmat([3; 4],ntotalion/2,1);

molid_ion = (1:ntotalion)' + molid_water(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put together coordinates
xyz_total = [xyz_water; xyz_ion];
type_total = [type_water;  type_ion];
molid_total = [molid_water; molid_ion];
index_total = (1:length(xyz_total))';
charge_total = zeros(length(xyz_total),1);

coordinates = [index_total molid_total type_total charge_total xyz_total];

% put together bonds, angles, and dihedrals
bonds = [(1:nbondsw)' ones(nbondsw,1) bondsw];
angles = [(1:nanglesw)' ones(nanglesw,1) anglesw];

% inputs are: coordinates, bonds, angles, dihedrals, impropers
% Write LAMMPS data file open the file with write permission
% coordinates = [index molid type charge x y z]
% bonds = [index_b type ID1 ID2]
% angles = [index_a type ID1 ID2 ID3]

name = strcat('L',num2str(round(L)),'nwat',num2str(nwater),'nion',num2str(ntotalion),'.data');
fid = fopen(name, 'w');

fprintf(fid,'%s\n\n','water box');

natoms = length(coordinates);
nbonds = length(bonds);
nangles = length(angles);
ndihedrals = 0;
nimpropers = 0;
fprintf(fid,'%d atoms\n',natoms);
fprintf(fid,'%d bonds\n',nbonds);
fprintf(fid,'%d angles\n',nangles);
fprintf(fid,'%d dihedrals\n',ndihedrals);
fprintf(fid,'%d impropers\n\n',nimpropers);

natomtypes = 4;
nbondtypes = 1;
nangletypes = 1;

fprintf(fid,'%d atom types\n',natomtypes);
fprintf(fid,'%d bond types\n',nbondtypes);
fprintf(fid,'%d angle types\n\n',nangletypes);

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
mass = [15.999; 1.008; 22.989; 35.453];
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

fclose(fid);
