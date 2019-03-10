clear all




nx = 150;
nz = nx;
c0 = 580; % c * dt <= dx
isx = 25;
isz = 25;
isx2 = 25;
isz2 = 75;
irx = 50;
irz = irx;
nt = 502;
dt = .0005;
dx = 1;

P = zeros(nx/dx, nz/dx);
Pold = zeros(nx/dx, nz/dx);
Pnew = zeros(nx/dx, nz/dx);

Pnew(isx,isz) = 3;
Pnew(isx2, isz2) = 3; 


j = (2:(nx-1));
k = (2:(nz-1));

for clock = 2:nt
t = clock;
Pold = P;
P = Pnew;
Pnew(j,k) = (dt^2/dx^2)*c0^2*(P((j+1),k) - 2*(P(j,k))+ ...
            P((j-1),k) + P(j,(k+1)) - 2*(P(j,k)) + P(j,(k-1)))+ ...
            2*P(j,k) - Pold(j,k); 
        
contour(Pnew)
axis equal
F(t-1) = getframe;
end

