L = 1; 
J = 101; 
dx = L/(J-1); 
M = 1.8;
f = 220;
T = M*(2*L*f)^2;
nskip = ceil(f*2*(J-1)/8192);
dt = 1/(8192*nskip);
r = 0.01; 
xp = 0.1*L;
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
S = zeros(1,floor(clockmax/nskip));
j = 2:(J-1);

for clock = 1:clockmax
    t = clock*dt;
    V(j) = V(j) + (dt/dx^2)*(T/M)*(H(j+1) - 2*H(j) + H(j-1)) + (dt/dx^2)*(r/M)*(V(j+1) - 2*V(j) + V(j-1));
    H = H + (dt*V);
    if (mod(clock, nskip) == 0)
        count = count + 1;
        S(count) = H(5);
%         plot(H);
%         ylim([-hp - .5, hp + .5]);
%         F(count) = getframe;
    end
end

soundsc(S)
% writerObj = VideoWriter('testWave.avi');
% open(writerObj);
% writeVideo(writerObj, F);
% close(writerObj);



