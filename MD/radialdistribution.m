function [g,rRange] = radialdistribution(pos,dr)
%
%out = particlecorrelation(pos)
%
%A function to compute the particle correlation relative to an ideal gas
%for a distribution of ions in a 1x1x1 cubic box
%

%compute the number of particles
[Np,Dim] = size(pos);

rRange=0:dr:sqrt(3);
rCount = zeros(1,length(rRange));

distance=zeros(Np,Np);

for i=1:Np
    i
    for j=i+1:Np
        distance(i,j) = norm(pos(i,:)-pos(j,:));
    end
end


g = 3*rCount/(Np^2)./((rRange+dr).^3-rRange.^3);