% 
% 

subplot(1,3,1) ;hold on ;
j1=1 %theta=0
 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(wA(:,:,j1));
hold on
 contourf(z,x,w); xlabel('z'), ylabel('x'),title('v_z')

j1=4

 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(wA(:,:,j1));
contourf(z,x,w) ;
axis('equal')



subplot(1,3,2) ;hold on ;
j1=1 %theta=0
 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(uA(:,:,j1));
hold on
 contourf(z,x,w); xlabel('z'), ylabel('x'), title('v_x')

% 
j1=4 %theta=0

 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(uA(:,:,j1));
contourf(z,x,w) ;
axis('equal')


subplot(1,3,3) ;hold on ;
j1=1
 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(pA(:,:,j1));
hold on
 contourf(z,x,w);xlabel('z'), ylabel('x'),title('p')


j1=4

 x=squeeze(X(:,:,j1));
 y=squeeze(Y(:,:,j1));
 z=squeeze(Z(:,:,j1));
 w=squeeze(pA(:,:,j1));
contourf(z,x,w) ;
axis('equal')










