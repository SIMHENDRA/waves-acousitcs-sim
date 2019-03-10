clear all
L = 5; 
J = 101; 
dx = L/(J-1); 
M = .1;
f = 220;
T = M*(2*L*f)^2;
nskip = ceil(f*2*(J-1)/8192);
dt = 1/(8192*nskip);
r = 0.01; 
xp = 0.2*L;
hp = 1.5;
V = zeros(1,J);
H = zeros(1,J);

for j = 1:J
    x = (j-1)*dx;
    if x < xp
        H(j) = hp*(x/xp);
    else
        H(j) = hp*(L-x)/(L-xp);
    end
end

tmax = 10;
clockmax = ceil(tmax/dt);
count = 0;
nt = clockmax;
if (clockmax > nt)    
    S = zeros(1,clockmax);
else
    S = zeros(1,clockmax + (nt-clockmax));
end
S2 = zeros(1,ceil(clockmax/nskip));
j = 2:(J-1);

for clock = 1:clockmax
    t = clock*dt;
    V(j) = V(j) + (dt/dx^2)*(T/M)*(H(j+1) - 2*H(j) + H(j-1)) + (dt/dx^2)*(r/M)*(V(j+1) - 2*V(j) + V(j-1));
    H = H + (dt*V);
    S(clock) = H(5);
    if (mod(clock, nskip) == 0)
        count = count + 1;
        S2(count) = H(5);
%         plot(H);
%         ylim([-hp - .5, hp + .5]);
%         F(count) = getframe;
    end
end

%can cause square
dt = dt;
c0 = 580;
dx = c0*dt +.5*c0*dt;



nx = 30*dx;
nz = nx; % c * dt <= dx
isx = 25;
isz = 25;
isx2 = 25;
isz2 = 75;
irx = 28;
irz = irx;
nt = clockmax;



P = zeros(ceil(nx/dx), ceil(nz/dx));
Pold = zeros(ceil(nx/dx), ceil(nz/dx));
Pnew = zeros(ceil(nx/dx), ceil(nz/dx));
Pnew(isx,isz) = S(1);
%can cause square

 
AcousSave = zeros(1,ceil(nt/nskip));

j1 = 2:ceil(nx/dx - 1);
k1 = 2:ceil(nz/dx - 1);

count = 2;

tracker = floor((.05*nt));
trackcount = 0;

for clock = 1:nt
t = clock;
Pold = P;
P = Pnew;

Pnew(j1,k1) = (dt^2/dx^2)*c0^2*(P((j1+1),k1) - 2*(P(j1,k1))+ ...
            P((j1-1),k1) + P(j1,(k1+1)) - 2*(P(j1,k1)) + P(j1,(k1-1)))+ ...
            2*P(j1,k1) - Pold(j1,k1);
Pnew(isx,isz) = S(count);
if (mod(clock, nskip) == 0)
    
    AcousSave(count) = Pnew(irx,irz);
    count = count + 1;
end
if (mod(clock, tracker) == 0)
    trackcount = trackcount + 1;
    disp([num2str(trackcount*5), '% complete'])
end
       
% contourf(Pnew)
% 
% 
% axis equal
% F(t) = getframe;
end

% xfft = abs(fftshift(fft(S)));
% df = 1:1:length(xfft);
% 
% plot(df,xfft)
% fig = figure;
% newfft = fft(AcousSave,2048);
% power = newfft.* conj(newfft) / 512;
% plot(power);
% testa = 1000*(0:1024)/2048;
% plot(testa,power(1:1025)) %frequency content





