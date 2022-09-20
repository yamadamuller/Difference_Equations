function [y, ny] = delay(x,nx,n0)
if n0>=0
    ny = nx(1):nx(length(nx))+n0;
    y = [zeros(1,n0) x];
else
    ny = nx(1)+n0:nx(length(nx));
    y = [x zeros(1,abs(n0))];
end

