clear all
L = 5; 
J = 101; 
dx = L/(J-1); 
M = 1;
f = 440;
T = M*(2*L*f)^2;
nskip = ceil(f*2*(J-1)/8192);
dt = 1/(8192*nskip);
r = 0.6; 
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


sizejawn = 20;
nx = sizejawn*dx;
nz = nx; % c * dt <= dx
isx = .3*sizejawn;
isx = floor(isx);
isz = .3*sizejawn;
isz = floor(isz);
isx2 = .7*sizejawn;
isx2 = floor(isx2);
isz2 = .3*sizejawn;
isz2 = floor(isz2);
irx = .5*sizejawn;
irx = floor(irx);
irz = irx;
irz = floor(irz);
nt = clockmax;
bdecay = 0.5;
airdecay = 0;


P = zeros(ceil(nx/dx), ceil(nz/dx));
Pold = zeros(ceil(nx/dx), ceil(nz/dx));
Pnew = zeros(ceil(nx/dx), ceil(nz/dx));
Pnew(isx,isz) = S(1);
Pnew(isx2,isz2) = S(1);
%can cause square

 
AcousSave = zeros(1,ceil(nt/nskip));
AcousSave(:,:,:) = zeros(ceil(nx/dx), ceil(nz/dx));



for ejawn = 1:10 %ceil(nt/nskip)
   for fjawn = 1:ceil(nx/dx)
       for gjawn = 1:ceil(nz/dx)
           AcousSave(ejawn,fjawn,gjawn) = 0;
       end
   end    
end







j1 = 2:ceil(nx/dx - 1);
k1 = 2:ceil(nz/dx - 1);

count = 2;

tracker = floor((.05*nt));
trackcount = 0;

for clock = 1:nt
t = clock;
Pold = P;
P = Pnew;
otherjawn = clock + 1;

Pnew(j1,k1) = ((dt^2/dx^2)*c0^2*(P((j1+1),k1) - 2*(P(j1,k1))+ ...
            P((j1-1),k1) + P(j1,(k1+1)) - 2*(P(j1,k1)) + P(j1,(k1-1)))+ ...
            2*P(j1,k1) - Pold(j1,k1));
           % ...
          %  - (((abs(isx - j1)/isx) + (abs(isz - k1)/isz))/2)*airdecay*P(j1,k1) ...
        
        
if (otherjawn <= nt)        
    Pnew(isx,isz) = S(otherjawn) + min(S2);
end
%Pnew(isx2,isz2) = S(count);

Pnew(min(j1),k1) = sqrt(bdecay)*Pnew(min(j1),k1);
Pnew(max(j1),k1) = sqrt(bdecay)*Pnew(max(j1),k1);
Pnew(j1,min(k1)) = sqrt(bdecay)*Pnew(j1,min(k1));
Pnew(j1,max(k1)) = sqrt(bdecay)*Pnew(j1,max(k1));

% plot(Pnew(isx,:))
% ylim([-2, 2]);
% G(t) = getframe;

if (mod(clock, nskip) == 0)
    
    AcousSave(count,:,:) = Pnew(:,:);
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

testerjawn = abs(fft(AcousSave));
tester2jawn = abs(fft(S2));

testerjawn(1) = 1;



