function out=Halton(index,base)
out=0;
f=1;
i=index;
while i>0
    f=f/base;
    out=out+f*mod(i,base)
    i=floor(i/base);
end
end