function coords = xyz_packing(dx,dy,dz,numberOfPoints,minAllowableDistance,ntrials)

    % generate coordinates betweeen 0 and dx,dy, and dz respectively
    x = dx*rand(1, ntrials);
    y = dy*rand(1, ntrials);
    z = dz*rand(1, ntrials);
    
    % Initialize first point.
    coords = zeros(numberOfPoints,3);

    % Try dropping down more points.
    count = 2;
    k = 1;
    while count <= numberOfPoints
        % Get a trial point.
        thisX = x(k);
        thisY = y(k);
        thisZ = z(k);
        k = k+1;
        % See how far is is away from existing keeper points.
        distances = sqrt((thisX-coords(:,1)).^2 + (thisY - coords(:,2)).^2 + (thisZ - coords(:,3)).^2);
        minDistance = min(distances);
        if minDistance >= minAllowableDistance
            coords(count,1) = thisX;
            coords(count,2) = thisY;
            coords(count,3) = thisZ;
            count = count + 1;
        end
    end
    
end

