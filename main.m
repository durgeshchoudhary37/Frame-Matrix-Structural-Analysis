% main script 
coordinates=load('coordi.txt');


% 3rd and 4th column contains values of EI and AE of corresponding element
% respectively
connectivity=load('connectivity.txt');


% known forces
kf=load('kf.txt');


% known displacements
kd=load('ku.txt');


% fixed end moments
fem=load('fem.txt');
a=size(coordinates);


% dof
m=3*a(1,1);


%size of global stiffness matrix denoted by n
n=m;
K_g=zeros(m,m);
a=size(connectivity);
m=a(1,1);

for i=1:m
  node_1=connectivity(i,1);
  node_2=connectivity(i,2);
  x1=coordinates(node_1,1);
  y1=coordinates(node_1,2);
  x2=coordinates(node_2,1);
  y2=coordinates(node_2,2);
  k=connectivity(i,3);
  AE=connectivity(i,4);

  k_e = ElementStiffness(x1,y1,x2,y2,k,AE);
  key=[1 2 3 4 5 6];
  value=[3*node_1-2 3*node_1-1 3*node_1 3*node_2-2 3*node_2-1 3*node_2];
  M=containers.Map(key,value);

  for j=1:6
      for k=1:6
      K_g(M(j),M(k))=K_g(M(j),M(k))+k_e(j,k);
      end
  end 
end


% construct matrixex for solution 
% defining two matrix c,d for lower division of K_g
a=size(kf);
b=size(kd);
p=a(1,1);
q=b(1,1);
C=zeros(p,n-p);
D=zeros(p,p);


% constructing C for the relation between known forces and known displacements
for i=1:p
    for j=1:n-p
        C(i,j)=K_g(kf(i,1),kd(j,1));
    end
end


% constructing D as a submatrix of the global stiffness matrix
for i=1:p
    for j=1:p
        D(i,j)=K_g(kf(i,1),kf(j,1));
    end
end

f=zeros(p,1);

for i=1:p
    f(i,1)=kf(i,2)-fem(kf(i,1));
    fem(kf(i,1)) = 0;
end

ku=zeros(q,1);

for i=1:q
    ku(i,1)=kd(i,2);
end

f=f-C*(ku);
D=inv(D);


% unknown displacements
d=D*f;

% displacement vector
disp=zeros(n,1);
i=1;
j=1;

for k=1:n
    if i<=q && kd(i,1)==k
        disp(k,1)=kd(i,2);
        i=i+1;
    else 
        disp(k,1)=d(j,1);
        j=j+1;
    end 
end


% force matrix
force=K_g*disp+fem;


% K_g
fprintf('The outputs for displacement and force are printed according to their serial degree of freedoms\n');


% displacement at node
disp


% forces
force


%Plot of undeformed structure 
s=size(connectivity);
e=s(1,1);

figure
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1),coordinates(node_2,1)];
    y=[coordinates(node_1,2),coordinates(node_2,2)];
    plot(x,y,'-ob','LineWidth',3);
    hold on 
end

title('ORIGINAL SHAPE (BEFORE DEFORMATION)');
hold off


%Plot of deformed structure 
figure
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1)+disp((node_1-1)*3+1,1),coordinates(node_2,1)+disp((node_2-1)*3+1,1)];
    y=[coordinates(node_1,2)+disp((node_1-1)*3+2,1),coordinates(node_2,2)+disp((node_2-1)*3+2,1)];
    plot(x,y,'-ob','LineWidth',3); 
    hold on
end

title('DEFORMED SHAPE','deformations in the scale of pow(10,-1)');
hold off


% bending moment and shear force diagrams
% representation of the member forces
M=zeros(e,6);
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
      x1=coordinates(node_1,1);
      y1=coordinates(node_1,2);
      x2=coordinates(node_2,1);
      y2=coordinates(node_2,2);
      k=connectivity(i,3);
      AE=connectivity(i,4);

      % gf->global forces
      gd=[disp(node_1*3-2,1);
          disp(node_1*3-1,1);
          disp(node_1*3,1);
          disp(node_2*3-2,1);
          disp(node_2*3-1,1);
          disp(node_2*3,1)];
      a=MF(x1,y1,x2,y2,k,AE,gd);
      for j=1:6
          M(i,j)=a(j,1);
      end
end


% plotting sfd & bmd
% (1). BMD
figure

for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1),coordinates(node_2,1)];
    y=[coordinates(node_1,2),coordinates(node_2,2)];
    plot(x,y,'-ob','LineWidth',3);
    hold on 
end

for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
   x1=coordinates(node_1,1);
   y1=coordinates(node_1,2);
   x2=coordinates(node_2,1);
   y2=coordinates(node_2,2);
   if x2==x1
    theta=pi/2;
    else
    theta=atan((y2-y1)/(x2-x1));
   end

    %c->costheta    s->sinetheta
    c=cos(theta);
    s=sin(theta);
    y=[coordinates(node_1,2)+M(i,3)*0.0001*c,coordinates(node_2,2)+M(i,6)*0.0001*c];
    x=[coordinates(node_1,1)+M(i,3)*0.0001*s,coordinates(node_2,1)+M(i,6)*0.0001*s];
    plot(x,y,'-og','LineWidth',2);
end

title('BENDING MOMENT DIAGRAM','BM in the scale of pow(10,-4)');


%(2). SFD
figure
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1),coordinates(node_2,1)];
    y=[coordinates(node_1,2),coordinates(node_2,2)];
    plot(x,y,'-ob','LineWidth',3);
    hold on
end

for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
   x1=coordinates(node_1,1);
   y1=coordinates(node_1,2);
   x2=coordinates(node_2,1);
   y2=coordinates(node_2,2);
   if x2==x1
    theta=pi/2;
    else
    theta=atan((y2-y1)/(x2-x1));
    end
    %c->costheta    s->sinetheta
    c=cos(theta);
    s=sin(theta);
    y=[coordinates(node_1,2)+M(i,2)*0.0001*c,coordinates(node_2,2)+M(i,5)*0.0001*c];
    x=[coordinates(node_1,1)+M(i,2)*0.0001*s,coordinates(node_2,1)+M(i,5)*0.0001*s];
    plot(x,y,'-og','LineWidth',2);
end
title('SHEAR FORCE DIAGRAM','SF in the scale of pow(10,-4)');
hold off




