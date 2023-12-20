function[xc,d,d2]=Chevigood(nr,Rout,rin)
%use c,x,dx,dx2,dp,u,v,du
nr1=nr+1;
pi=acos(-1.0);

c=zeros(nr1);
for jm=2:nr
	c(jm)=1.0;
end
c(1)=2.0;
c(nr1)=2.0;

c1=4.0;
c2=1.0+2.0*(c1/Rout);

x=zeros(1,nr1);
xc=zeros(1,nr1);
dx=zeros(1,nr1);
dx2=zeros(1,nr1);
d=zeros(nr1,nr1);
d2=zeros(nr1,nr1);

x(1)=1.0;
xc(1)=rin;
dx(1)=-2.d0/Rout ;
	dx2(1)=0.d0;
for j=1:nr
    x(j+1)=cos((j*pi)/nr);
    xc(j+1)=rin+0.5d0*(1.d0-x(j+1))*Rout; 
    dx(j+1)=-2.d0/Rout;
    dx2(j+1)=0;
end

for l=0:nr
    for j=0:nr
        jm=j+1;
        lm=l+1;
        if (j==l)
            if( jm>1 && jm<nr1)
                d(lm,jm)=-x(lm)/(2.0*(1-(x(jm)^2.0)));
            end
        else
            d(lm,jm)=c(lm)*(-1.0)^(l+j)/(c(jm)*(x(lm)-x(jm)));
        end
    end
end
d(1,1)=(2.0*(nr^2.0)+1.0)/6.0;
d(nr1,nr1)=-d(1,1);



%physical dlj2 and dlj
dp=zeros(nr1,nr1);
for l=1:nr1
    for j=1:nr1
       
        d(l,j)=d(l,j)*dx(l);
    end
end
d2=d*d;
