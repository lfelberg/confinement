function [bonds,angles,dihedrals,impropers,dtype] = graph_topo(coords,b,Ly,Lz)
% Generate topology of a graphene sheet . The
% convention is zigzag in z and armchair in y. The sheet is in the yz
% plane.

ngraph = length(coords);
% Reduce to 2D
coords = coords(:,2:3);

% Calculate distances
dist = zeros(ngraph,ngraph);
for ii=1:ngraph
    ref = coords(ii,:);
    for j=ii+1:ngraph
        tar = coords(j,:);
        incry = abs(tar(1)-ref(1));
        incrz = abs(tar(2)-ref(2));
        if incry > Ly/2
            incry = Ly - incry;
        elseif incrz > Lz/2
            incrz = Lz - incrz;
        end
        dist(ii,j) = sqrt(incry^2+incrz^2);
    end
end

% Bonds: every atom at distance b from each other is bonded. 3 bonds/atom. 
ibond = find(dist<b+0.1&dist>0.1);
[row,col] = ind2sub(size(dist),ibond);
bonds = [row,col];

% Angles
angles = zeros(3*ngraph,3);
for ii=1:ngraph
    % find the 3 bonds associated to atom ii
    iprov = find(bonds==ii);
    [row,~] = ind2sub(size(bonds),iprov);
    % find the 3 indexes of the atoms bonded to ii 
    bb = bonds(row,:);
    cc = bb(bb~=ii);
    % define the 3 angles related to ii
    angles(3*ii-2:3*ii,:) = [ cc(1) ii cc(2); cc(1) ii cc(3); cc(2) ii cc(3)];
end

% Dihedrals: each bond has 4 dihedrals associated
nbonds = length(bonds);
dihedrals = zeros(4*nbonds,4);
dtype = zeros(4*nbonds,1);

for ii=1:nbonds
    % atoms in the bond
    at2 = bonds(ii,1);
    at3 = bonds(ii,2);

    % find the neighbors of at2 
    [aj,bj] = find(bonds==at2);
    bj(bj==1)=0; bj(bj==2) = 1; bj(bj==0) = 2;
    otherj = bonds(sub2ind(size(bonds),aj,bj));    
    cc2 = otherj(otherj~=at3);
    
    % find the neighbors of at3
    [aj,bj] = find(bonds==at3);
    bj(bj==1)=0; bj(bj==2) = 1; bj(bj==0) = 2;
    otherj = bonds(sub2ind(size(bonds),aj,bj));    
    cc3 = otherj(otherj~=at2);
    
    % define the 4 dihedrals per bond
    dihedrals(4*ii-3:4*ii,:) = [cc2(1) at2 at3 cc3(1);
                                cc2(1) at2 at3 cc3(2);
                                cc2(2) at2 at3 cc3(1);
                                cc2(2) at2 at3 cc3(2)];
                            
    % dihedral types: type 1 \_/ ; type 2 \_\
    incry = coords(cc2(1),1)-coords(cc3(1),1);
    incrz = coords(cc2(1),2)-coords(cc3(1),2);
    if incry > Ly/2
        incry = Ly - incry;
    elseif incrz > Lz/2
        incrz = Lz - incrz;
    end
    ddih1 = sqrt(incry^2+incrz^2);
    if ddih1 > 2*b+0.1
        tdih1 = 2;
        tdih2 = 1;
    else
        tdih1 = 1;
        tdih2 = 2;
    end
    
    incry = coords(cc2(2),1)-coords(cc3(1),1);
    incrz = coords(cc2(2),2)-coords(cc3(1),2);
    if incry > Ly/2
        incry = Ly - incry;
    elseif incrz > Lz/2
        incrz = Lz - incrz;
    end
    ddih3 = sqrt(incry^2+incrz^2);
    if ddih3 > 2*b+0.1
        tdih3 = 2;
        tdih4 = 1;
    else
        tdih3 = 1;
        tdih4 = 2;
    end    
    dtype(4*ii-3:4*ii,1) = [tdih1 tdih2 tdih3 tdih4];
    
    
end

% Impropers: each atoms has 1 improper
impropers = zeros(ngraph,4);

for ii=1:ngraph

    % atoms in the bond
    at_i = ii;
    
    % find the atoms bonded to ati
    [aj,bj] = find(bonds==at_i);
    bj(bj==1)=0; bj(bj==2) = 1; bj(bj==0) = 2;
    other = bonds(sub2ind(size(bonds),aj,bj));
    
    impropers(ii,:) = [at_i other(1) other(2) other(3)];
    
    
end
