

ns = nzA;
ntheta =nx;
neta = nrA;

% Discretize the parameters
s = z0A;
theta = xx;
eta = r0A;


% Define the radius of the center curve
Rc = 2 * L / pi;

% Initialize arrays for coordinates
[X, Y, Z] = deal(zeros(ns, neta, ntheta));

% Calculate the coordinates
for i = 1:ns
for j = 1:neta
for k = 1:ntheta
% Parametric equations for the volume of the tube
X(i, j, k) =(cos(s(i))*(2*L - pi*eta(j)*cos(theta(k))))/pi;
Y(i, j, k)= -eta(j)*sin(theta(k));
Z(i, j, k)=(sin(s(i))*(2*L - pi*eta(j)*cos(theta(k))))/pi;
end
end
end

% % Plot the mesh
% %figure;
% hold on;
% for i = 1:ns
% for j = 1:neta
% plot3(squeeze(X(i, j, :)), squeeze(Y(i, j, :)), squeeze(Z(i, j, :)), 'b');
% end
% end
% for j = 1:neta
% for k = 1:ntheta
% plot3(squeeze(X(:, j, k)), squeeze(Y(:, j, k)), squeeze(Z(:, j, k)), 'r');
% end
% end
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Mesh of the Cylindrical Tube Bent into a Semicircle');
% 
% %Midle plane \thetha=0
% surface(squeeze(X(:, :, 1)),squeeze(Y(:, :, 1)),squeeze(Z(:, :, 1)));
% %entrace s=0
% surface(squeeze(X(1, :, :)),squeeze(Y(1, :, :)),squeeze(Z(1, :, :)));
% %wall eta=1;
% surface(squeeze(X(:, neta, :)),squeeze(Y(:, neta, :)),squeeze(Z(:, neta, :)));
% grid on;
% axis equal;
% hold off;

X = permute(X, [2 1 3]);
Y = permute(Y, [2 1 3]);
Z = permute(Z, [2 1 3]);