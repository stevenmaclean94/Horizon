clear
clc

JD = ymdUT2JD(2004, 3, 3, 4.5);

p1 = lla2eci(0,0,0,JD)';

i = 1;
for lon = 0:360
    p2 = lla2eci(0,lon,0,JD)';
    l(i) = acos(p1'*p2)*180/pi;
    i = i+1;
end

plot(l)
hold on
plot(0:360,'r')