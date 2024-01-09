% fuction for element stiffness
function K=ElementStiffness( x1,y1 ,x2 ,y2 ,k,AE)
l=sqrt((x2-x1)^2+(y2-y1)^2);
if x2==x1
    theta=pi/2;

else
theta=atan((y2-y1)/(x2-x1));
end

%c->cos(theta)    s->sin(theta)
c=cos(theta);
s=sin(theta);

T=[c s 0 0 0 0;
   -s c 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 c s 0;
   0 0 0 -s c 0;
   0 0 0 0 0 1];

A = k*[AE/(k*l) 0 0 -AE/(k*l) 0 0;
    0 12/(l*l*l) 6/(l*l) 0 -12/(l*l*l) 6/(l*l);
    0 6/(l*l) 4/l 0 -6/(l*l) 2/l;
    -AE/(k*l) 0 0 AE/(k*l) 0 0;
    0 -12/(l*l*l) -6/(l*l) 0 12/(l*l*l) -6/(l*l);
    0 6/(l*l) 2/l 0 -6/(l*l) 4/l;];

K=T'*A* T;
end

