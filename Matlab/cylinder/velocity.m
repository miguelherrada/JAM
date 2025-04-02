function velocity(G1,G2,F1,F2,w1,w2)
c1=max(max(w1));
c2=max(max(w2));
c01=max(c1,c2);

c1=min(min(w1));
c2=min(min(w2));
c02=min(c1,c2);
v=  linspace(c02,c01,10);

contourf(G1,F1,w1,v)
 hold on
contourf(G2,F2,w2,v)

xlabel('x')
ylabel('y')

