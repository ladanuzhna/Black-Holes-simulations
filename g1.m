%Quick converters
parsec_to_AU=inline('parsec*206264.80624538','parsec')
AU_to_parsec=inline('AU*0.0000048481367817234','AU')
km2_to_AU2 = inline('km2*4.46837*pow(10,-17)','km2')

% %function definitions
syms r(T) p(T) T Y t

G= 4.30091252525*(10^(-3)); %grav constant, in (parsec*km^2)/(Ms*sec^2) 
M = 10^6; %mass of the Black Hole, in solar masses
m = 10^5; %mass of the object, in solar masses
reduced = m*M / (m+M);
c = 0.0020053761; % AU/sec, speed of 
rs = 2*G*M / c^2; %Schwarzchild radius
PI = 3.14;
h=m;

file1 = fopen('C:\Users\ladan\Documents\Research_Project\energy.csv','r');
file2 = fopen('C:\Users\ladan\Documents\Research_Project\angular_momentum.csv','r');
E = fscanf(file1, '%f', [1 Inf]);
L = fscanf(file2, '%f', [1 Inf]);
fclose(file1);
fclose(file2);

ode1 = diff(r,T) == sqrt((E^2/c^2 - c^2) + 2*parsec_to_AU(G)*M/r*(1+L^2/((c^2)*(r^2)))-(L^2)/(r^2));
ode2 = diff(p,T) == (1/r^2) * L;


% %function setup
[f1, Subs]= odeToVectorField([ode1 ode2]);
F1= matlabFunction(f1, 'Vars',{T,Y})

%
r_0= 1.496*10^8; %In KM
p_0 = PI;
t_0 = 0;
span = [0 0.008];
my_val = sqrt((E^2/c^2 - c^2) + 2*parsec_to_AU(G)*M/r_0*(1+L^2/((c^2)*(r_0^2)))-(L^2)/(r_0^2));
conditions = [r_0;p_0]; 
options = odeset('MaxStep', 9.5*10^(-7),'AbsTol',10^(-7));
[rsol,psol] = ode45(@(T,Y)F1(T,Y),span,conditions,options);
xfile = fopen('C:\Users\ladan\Documents\Research_Project\xgeodesic.csv','w');
yfile = fopen('C:\Users\ladan\Documents\Research_Project\ygeodesic.csv','w');
for i=1:length(rsol)
    PHI = psol(i,2);
%     if sin(PHI) > 0 
%         ysign = 1
%     else
%         ysign = -1
%     end
%     if cos(PHI) > 0 
%         xsign = 1
%     else
%         xsign = -1
%     end
    x(i) = psol(i,1)*cos(PHI);
    y(i) = psol(i,1)*sin(PHI);
    %fprintf(xfile,'%d\n',x(i));
    %fprintf(yfile,'%d\n',y(i));
end


x = smoothdata(x,'gaussian',30);
y = smoothdata(y,'gaussian',30);

for i=1:length(x)
    fprintf(xfile,'%d\n',x(i));
    fprintf(yfile,'%d\n',y(i));
end

fclose(xfile)
fclose(yfile)
% time_step= span(2)/length(rsol);
% time_vector = 0:time_step:span(2);
% figure
% plot(psol(:,2),rsol)
% title('Position and angle plotted against each other')
% 
% figure
% plot(time_vector(1:(length(time_vector)-1)),rsol)
% title('X axis - time, Y axis - position')
% 
% figure
% plot(time_vector(1:(length(time_vector)-1)),psol(:,2))
% title('X axis - time, Y axis - angle')

%figure
%animatedLine(x,y)
%comet(x,y,0.01)

