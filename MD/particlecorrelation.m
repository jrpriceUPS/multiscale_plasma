function out = radialdistribution(pos,dr)
%
%out = particlecorrelation(pos)
%
%A function to compute the particle correlation relative to an ideal gas
%for a distribution of ions in a 1x1x1 cubic box
%

%compute the number of particles
[Np,Dim] = size(pos);

distances = zeros(Np*(Np-1)/2,1);

index=1;
for i=1:Np
    for j=i+1:Np
        distances(index) = sum((pos(i,:)-pos(j,:)).^2);
        index = index+1;
    end
end

maxr = max(distances);
minr = min(distances);

rRange = minr:dr:maxr;

rCount = zeros(length(rRange)-1,1);

for k=1:length(rCount)
   rCount(k) =  sum(distances>=rRange(k) & distances<rRange(k+1));
end

