clear
close all
sides = 8;
deg_inc = 360/sides;
beta = 57;
deg0 = deg_inc:deg_inc:360;
x = sind(deg0);
y = cosd(deg0);
z = ones([1 sides])*cosd(beta);
figure
hold on
grid on
% scatter3(x,y,z)
xlabel('X')
ylabel('Y')
zlabel('Z')
for i = 1:sides
    vec1 = [-x(i) -y(i) z(i)];
    vec2 = [0 0 -1];
    vec3 = cross(vec2,vec1);
    vec4 = cross(vec3,vec1);
    quiver3(x(i),y(i),z(i),vec3(1), vec3(2), vec3(3))
    quiver3(x(i),y(i),z(i),vec1(1),vec1(2),vec1(3))
    quiver3(x(i),y(i),z(i),vec4(1), vec4(2), vec4(3))
end
viscircles([0 0],1)
%%
align = [0.65	0.537	0.537
0.65	-0.537	0.537
-0.65	-0.537	0.537
-0.65	0.537	0.537]
x = align(:,1)
y = align(:,2)
z = align(:,3)
figure
hold on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
for i = 1:4
    quiver3(x(i),y(i),z(i),x(i),y(i),z(i))
end