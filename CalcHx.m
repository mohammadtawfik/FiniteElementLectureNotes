function H=CalcHx(x,y)

H=zeros(1,4);

f1=0;
f2=1;
g1=1;
g2=y;

H(1)=f1*g1;
H(2)=f2*g1;
H(3)=f1*g2;
H(4)=f2*g2;
