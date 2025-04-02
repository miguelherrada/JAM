figure;

% Primer subplot
subplot(2,1,1)
hold on  
for i=1:3:nr(1)
    plot(G1(i,:),F1(i,:),'g-');
end
for j=1:3:nz(1)
    plot(G1(:,j),F1(:,j),'g-');
end

for i=1:3:nr(2)
    plot(G2(i,:),F2(i,:),'r-');
end
for j=1:3:nz(2)
    plot(G2(:,j),F2(:,j),'r-');
end  
xlabel('x')
ylabel('y')
axis equal % Mantiene proporciones correctas
xlim([min(G1(:)), max(G1(:))])
ylim([min(F1(:)), max(F2(:))])

% Segundo subplot
subplot(2,1,2)
hold on  
for i=1:3:nr(1)
    plot(z1,r1(i,:),'r-');
end
for j=1:3:nz(1)
    plot(z1(:,j),r1(:,j),'r-');
end

for i=1:3:nr(2)
    plot(z1(i,:),-r1(i,:),'g-');
end
for j=1:3:nz(2)
    plot(z1(:,j),-r1(:,j),'g-');
end  

% Ajustar tamaño para que ambos sean iguales en la figura
subplot(2,1,1)
pbaspect([1 1 1]) % Ajusta la relación de aspecto (1:1)
subplot(2,1,2)
pbaspect([1.8 1 1])

% Ajustar márgenes y tamaños
set(gcf, 'Position', [100, 100, 600, 800]); % Ajusta tamaño de la figura

xlabel('z_0')
ylabel('r_0')