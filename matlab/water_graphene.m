% start fresh
clear
close all
clc

% inputs: dimensions of the system (approximate)
Lx = 12.75;
Ly = 46;
Lz = Ly;
ntrials = 1e7; % for water packing

% monomer TIP4-Ew
roh = 0.9572;
theta = 104.52*(pi/180);
% graphene
b = 1.4;

% correct box dimensions to fit graphene with PBCs; armchair in y, zigzag in z
prov = Ly/(b*(1+cos(pi/3)));
ny = 2*round(prov/2);
Ly = ny*(b*(1+cos(pi/3)));

prov = Lz/(b*sin(pi/3));
nz = 2*round(prov/2);
Lz = nz*(b*sin(pi/3));

% WATER

% preliminary calculations for filling the box
doo = 2.6; % minimum oxygen-oxygen separation
dex = 2; % excluded distance between graphene and water
slack = 2;
nwatercent = round(1/1e21/1.66e-24/18*((Lx-2*dex)*(Ly-slack)*(Lz-slack)/1e3));
nwaterside = round(nwatercent/2);

coordsLeft = xyz_packing(Lx/2-dex-slack,Ly-slack,Lz-slack,nwaterside,doo,ntrials);
coordsLeft(:,1) = coordsLeft(:,1) - Lx/2 + slack;
coordsCent = xyz_packing(Lx-2*dex,Ly-slack,Lz-slack,nwatercent,doo,ntrials);
coordsCent(:,1) = coordsCent(:,1) + dex;
coordsRight = xyz_packing(Lx/2-dex-slack,Ly-slack,Lz-slack,nwaterside,doo,ntrials);
coordsRight(:,1) = coordsRight(:,1) + Lx + dex;

coordsOxy = [coordsLeft; coordsCent; coordsRight];
nwater = length(coordsOxy);

% atom types
type_water = repmat([1;2;2],nwater,1);

% coordinates = [index molid type charge x y z]
% monomer (OHH) 
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

% molecule index
molid_water = zeros(3*nwater,1);
for i=1:nwater
    molid_water(3*i-2:3*i) = zeros(3,1)+i;
end

% bonds and angles

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

% graphene %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ngraph = nz*ny; % number of graphene atoms
incrz = b*sin(pi/3);
incry = b*(1+cos(pi/3));

y = zeros(ngraph,1);
z = zeros(ngraph,1);
x = zeros(ngraph,1);
count = 1;
for j=1:ny
    for i=1:nz      
        z(count) = (i-1)*incrz;
        if mod(j,2)==1 % j is odd
            if mod(i,2)==1 % i is odd
                incrincry = b*cos(pi/3);
            else % i is even
                incrincry =0;
            end
            y(count) = incrincry + (j-1)*incry;
        else % j is even
            if mod(i,2)==1 % odd
                incrincry = 0;
            else
                incrincry = b*cos(pi/3);
            end
            y(count) = incrincry + (j-1)*incry;
        end       
        count = count + 1;
    end
end

xyz_graph_l = [x y z];
xyz_graph_r = [x+Lx y z];

type_graph = zeros(2*ngraph,1)+3; % O-type1, H-type2
molid_graph = [ones(ngraph,1)+nwater; ones(ngraph,1)+nwater+1];

[bondsg_l,anglesg_l,dihedralsg_l,impropersg_l,dihtypel] = graph_topo(xyz_graph_l,b,Ly,Lz);
[bondsg_r,anglesg_r,dihedralsg_r,impropersg_r,dihtyper] = graph_topo(xyz_graph_r,b,Ly,Lz);

bondsg = [bondsg_l; bondsg_r+ngraph];
anglesg = [anglesg_l; anglesg_r+ngraph];
dihedralsg = [dihedralsg_l; dihedralsg_r+ngraph];
impropersg = [impropersg_l; impropersg_r+ngraph];

nbondsg = length(bondsg);
nanglesg = length(anglesg);
ndihedralsg = length(dihedralsg);
nimpropersg = length(impropersg);

% put together coordinates
xyz_total = [xyz_water; xyz_graph_l; xyz_graph_r];
index_total = (1:length(xyz_total))';
charge_total = zeros(length(xyz_total),1);
type_total = [type_water; type_graph];
molid_total = [molid_water; molid_graph];
coordinates = [index_total molid_total type_total charge_total xyz_total];

% put together bonds, angles, and dihedrals
nbonds = nbondsw + nbondsg;
index_b = (1:nbonds)';
type_b = [ones(nbondsw,1); ones(nbondsg,1)+1];
bonds = [index_b type_b [bondsw; bondsg+nwater*3]];

nangles = nanglesw + nanglesg;
index_a = (1:nangles)';
type_a = [ones(nanglesw,1); ones(nanglesg,1)+1];
angles = [index_a type_a [anglesw; anglesg+nwater*3]];

ndihedrals = ndihedralsg;
index_d = (1:ndihedrals)';
type_d = [dihtypel; dihtyper];
dihedrals = [index_d type_d dihedralsg+nwater*3];

nimpropers = nimpropersg;
index_i = (1:nimpropers)';
type_i = ones(nimpropers,1);
impropers = [index_i type_i impropersg+nwater*3];

% inputs are: coordinates, bonds, angles, dihedrals, impropers
% Write LAMMPS data file open the file with write permission
% coordinates = [index molid type charge x y z]
% bonds = [index_b type ID1 ID2]
% angles = [index_a type ID1 ID2 ID3]

name = strcat('d',num2str(Lx),'L',num2str(Ly),'_nobenzene.data');
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

natomtypes = 3;
nbondtypes = 2;
nangletypes = 2;
ndihedraltypes = 2;
nimpropertypes = 1;
fprintf(fid,'%d atom types\n',natomtypes);
fprintf(fid,'%d bond types\n',nbondtypes);
fprintf(fid,'%d angle types\n',nangletypes);
fprintf(fid,'%d dihedral types\n',ndihedraltypes);
fprintf(fid,'%d improper types\n\n',nimpropertypes);

% empty box padding in the x-dimension
minx = -Lx/2; 
miny = 0;
minz = 0;
maxx = Lx+Lx/2;
maxy = Ly;
maxz = Lz;

fprintf(fid,'%.2f %.2f xlo xhi\n',minx,maxx);
fprintf(fid,'%.2f %.2f ylo yhi\n',miny,maxy);
fprintf(fid,'%.2f %.2f zlo zhi\n\n',minz,maxz);

% O=1, H=2, C=3
fprintf(fid,'%s\n\n','Masses');
mass = [15.9994; 1.008; 12];
id = [1; 2; 3];
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
