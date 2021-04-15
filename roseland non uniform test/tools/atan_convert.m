function v=atan_convert(y,x)
v=nan;
if x>0
    v=atan(y/x);
end
if y>=0 && x<0
    v=pi+atan(y/x);
end
if y<0 && x<0
    v=-pi+atan(y/x);
end
if y>0 && x==0
    v=pi/2;
end
if y<0 && x==0
    v=-pi/2;
end
if v<0
    v=v+2*pi;
end
end