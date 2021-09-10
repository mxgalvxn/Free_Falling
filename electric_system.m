clear vars;
clf;
set(gcf, 'Position', get(0,'Screensize')); warning('off','all');
% Valores iniciales


m0 = 4*pi*(10^-7); % permeabilidad del vacío
miu3 = m0^3;
M = 0.01;                % masa (kg)
r = 2;                   % radio (m)
r2 = r^2;
r6 = r^6;
S = pi*r^2;              % area (m^2) 
%R = 0.00009              % (acercada)
R = 0.0009;             % resistencia (amplia)
g = 9.81;                % aceleración (m/s^2)
U = 1000000;             % miu
U2 = U^2;
x0 = 60;
dt = 0.01; 
h = dt;
v0 = 0;
t0 = 0;
tf = 6;
fi0 = 0;



xplot=t0:h:tf;
yplot = 0.*xplot;

% Aceleración inicial con U_B
a = @(x) -g + (27*(x^3)*S*(miu3)*(U2)*(r6)*(v0^2))/(8*M*(R^2)*(x^2+r2)^(15/2));

a0 = -g + (27*(x0^3)*S*(miu3)*(U2)*(r6)*(v0^2))/(8*M*(R^2)*(x0^2+r2)^(15/2));


%Fem inducida
f = @(x) (3*m0*S*U*(r2)*(x)*(v0))/(2*(x^2 + r2)^(5/2));

f0 = (3*m0*S*U*(r2)*(x0)*(v0))/(2*(x0^2 + r2)^(5/2));


% Creación de vectores de posición, velocidad, aceleración y tiempo
t = t0:h:tf;
X = zeros(1,length(t));
V = X;
A = X;
T = X;
Fem = X;

A(1) = a0;
V(1) = v0;
X(1) = x0;
T(1) = t0;
Fem(1) = fi0;
    
%empieza en y(2)
for i = 1:(length(t)-1)
    
    k1 = -g + (27*(X(i)^3)*S*miu3*U2*(r6)*(V(i)^2))/(8*M*(R^2)*(X(i)^2 + r2)^(15/2));
    k2 = -g + (27*((X(i)+0.5*h)^3)*S*miu3*U2*(r6)*(V(i)^2))/(8*M*(R^2)*((X(i)+0.5*h)^2 + r2)^(15/2));
    k3 = -g + (27*((X(i)+0.5*h)^3)*S*miu3*U2*(r6)*(V(i)^2))/(8*M*(R^2)*((X(i)+0.5*h)^2 + r2)^(15/2));
    k4 = -g + (27*((X(i)+h)^3)*S*miu3*U2*(r6)*(V(i)^2))/(8*M*(R^2)*((X(i)+h)^2 + r2)^(15/2));

    V(i+1) = V(i) + (1/6).*(k1+(2.*k2)+(2.*k3)+k4).*h;

    A(i) = k1;

    l1 = V(i) + A(i)*dt;
    l2 = V(i) + A(i)*(dt+0.5*h);
    l3 = V(i) + A(i)*(dt+0.5*h);
    l4 = V(i) + A(i)*(dt+h);

    X(i+1) = X(i) + (1/6) * (l1 + 2*l2 + 2*l3 + l4) * h;
    
    T(i+1) = T(i) + dt;
    
    
    b1 = (3*m0*S*U*(r2+50)*(X(i))*(V(i)))/(2*(X(i)^2 + r2+50)^(5/2));
    b2 = (3*m0*S*U*(r2+50)*(X(i)+0.5*h)*(V(i)))/(2*((X(i)+0.5*h)^2 + r2+50)^(5/2));
    b3 = (3*m0*S*U*(r2+50)*(X(i)+0.5*h)*(V(i)))/(2*((X(i)+0.5*h)^2 + r2+50)^(5/2));
    b4 = (3*m0*S*U*(r2+50)*(X(i)+h)*(V(i)))/(2*((X(i)+h)^2 + r2+50)^(5/2));
    
    Fem(i+1) = Fem(i) + (1/6).*(b1 + 2*b2 + 2*b3 + b4) * h;
    
end



%Se grafica el resultado
t = tiledlayout(1,4);
nexttile
plot(T,X,'b')
hold on;
plot(xplot, yplot);
grid on;
axis([t0 tf -10 x0+5]);
title('Posición')
xlabel('Tiempo (s)')
ylabel('Posicion (m)')

nexttile
plot(T,V,'b')
hold on;
plot(xplot, yplot);
grid on;
axis([t0 tf-1 -75 0]);
title('Velocidad')
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')

nexttile
plot(T,A,'b')
hold on;
plot(xplot, yplot);
grid on;
axis([t0 tf -100 600]);
title('Aceleración')
xlabel('Tiempo (s)')
ylabel('Aceleración (m/s2)')

nexttile
plot(T,Fem, 'b')
hold on;
plot(xplot, yplot);
grid on;
axis([t0 tf -1.5 1.5]);
title('Fem inducida')
xlabel('Tiempo (s)')
ylabel('Corriente (A)')
