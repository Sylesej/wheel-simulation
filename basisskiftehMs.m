function hMs = basisskiftehMs(K)

hMs = zeros(size(K));

a = size(K)/6;
a = a(1);
for n = 1:a
    hMs(n*6-5:n*6,n*6-5:n*6) = delFunk(-2*pi/a*(n-1));
end
    