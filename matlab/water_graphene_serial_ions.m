% start fresh
clear
close all
clc

% inputs: dimensions of the system (approximate)
d = 9;
Ly = 48;
Lz = Ly;
ntimes = 2;

% nsolutey*nsoutez has to be even!
nsolutey = 4;
nsolutez = 8;
 
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
nwatercent = round(1/1e21/1.66e-24/18*((d-2*dex)*(Ly-slack)*(Lz-slack)/1e3));
nwaterside = round(nwatercent/2);

coordsOxy = xyz_packing(d-2*dex,Ly-slack,Lz-slack,nwatercent,doo,ntrials);
coordsOxy(:,1) = coordsOxy(:,1) + dex;

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

xyz_graph = [x y z];

type_graph = zeros(ngraph,1)+3; % O-type1, H-type2
molid_graph = ones(ngraph,1)+nwater;

[bondsg,anglesg,dihedralsg,impropersg,dihtype] = graph_topo(xyz_graph,b,Ly,Lz);

nbondsg = length(bondsg);
nanglesg = length(anglesg);
ndihedralsg = length(dihedralsg);
nimpropersg = length(impropersg);

% ions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
incrposy = Ly/nsolutey;
incrposz = Lz/nsolutez;
ntotalion = nsolutey*nsolutez;

count = 1;
xyz_ion = zeros(ntotalion,3);
for ii=1:nsolutey
    for jj=1:nsolutez
        if mod(count,2)==0
            incrz = 0.70*d;
        else
            incrz = 0.30*d;
        end
        xyz_ion(count,:) = [incrz (ii-1)*incrposy (jj-1)*incrposz];
        count = count + 1;
    end
end

type_ion = repmat([4; 5],ntotalion/2,1);

molid_ion = (1:ntotalion)' + molid_graph(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put together coordinates
xyz_comb = [xyz_water; xyz_graph; xyz_ion];
nwatgraphbenz = size(xyz_comb,1);
type_comb = [type_water; type_graph; type_ion];
molid_comb = [molid_water; molid_graph; molid_ion];

xyz_total = [];
type_total =[];
molid_total = [];
for ii=1:ntimes
    xyz_total = [xyz_total; xyz_comb(:,1)+(ii-1)*d xyz_comb(:,2:3)];
    type_total = [type_total; type_comb];
    molid_total = [molid_total; molid_comb+(ii-1)*max(molid_comb)];
end
index_total = (1:length(xyz_total))';
charge_total = zeros(length(xyz_total),1);
coordinates = [index_total molid_total type_total charge_total xyz_total];

% put together bonds, angles, and dihedrals
nbonds = (nbondsw + nbondsg)*ntimes;
index_b = (1:nbonds)';
type_b = repmat([zeros(nbondsw,1)+1; zeros(nbondsg,1)+2],ntimes,1);
bondsAt = [];
for ii=1:ntimes
    pluswat = (ii-1)*nwatgraphbenz;
    plusgraph = pluswat + length(xyz_water);
    bondsAt = [bondsAt; bondsw+pluswat; bondsg+plusgraph];
end
bonds = [index_b type_b bondsAt];

nangles = (nanglesw + nanglesg)*ntimes;
index_a = (1:nangles)';
type_a = repmat([zeros(nanglesw,1)+1; zeros(nanglesg,1)+2],ntimes,1);
anglesAt = [];
for ii=1:ntimes
    pluswat = (ii-1)*nwatgraphbenz;
    plusgraph = pluswat + length(xyz_water);
    anglesAt = [anglesAt; anglesw+pluswat; anglesg+plusgraph];
end
angles = [index_a type_a anglesAt];

ndihedrals = (ndihedralsg)*ntimes;
index_d = (1:ndihedrals)';
type_d = repmat(dihtype,ntimes,1);
dihedralsAt = [];
for ii=1:ntimes
    pluswat = (ii-1)*nwatgraphbenz;
    plusgraph = pluswat + length(xyz_water);
    dihedralsAt = [dihedralsAt; dihedralsg+plusgraph];
end
dihedrals = [index_d type_d dihedralsAt];

nimpropers = (nimpropersg)*ntimes;
index_i = (1:nimpropers)';
type_i = repmat([zeros(nimpropersg,1)+1],ntimes,1);
impropersAt = [];
for ii=1:ntimes
    pluswat = (ii-1)*nwatgraphbenz;
    plusgraph = pluswat + length(xyz_water);
    plusbenz = plusgraph + length(xyz_graph);
    impropersAt = [impropersAt; impropersg+plusgraph];
end
impropers = [index_i type_i impropersAt];

% inputs are: coordinates, bonds, angles, dihedrals, impropers
% Write LAMMPS data file open the file with write permission
% coordinates = [index molid type charge x y z]
% bonds = [index_b type ID1 ID2]
% angles = [index_a type ID1 ID2 ID3]

name = strcat('d',num2str(d),'L',num2str(round(Ly)),'nwat',num2str(nwater),'nrep',num2str(ntimes),'nion',num2str(ntotalion),'.data');
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

natomtypes = 5;
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
minx = -d/2; 
miny = 0;
minz = 0;
maxx = (ntimes-1)*d+d/2;
maxy = Ly;
maxz = Lz;

fprintf(fid,'%.2f %.2f xlo xhi\n',minx,maxx);
fprintf(fid,'%.2f %.2f ylo yhi\n',miny,maxy);
fprintf(fid,'%.2f %.2f zlo zhi\n\n',minz,maxz);

% O=1, H=2, C=3
fprintf(fid,'%s\n\n','Masses');
mass = [15.9994; 1.008; 12; 39.0983; 35.453];
id = [1; 2; 3; 4; 5];
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
