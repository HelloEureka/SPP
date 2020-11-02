function [B, L, H] = XYZ2BLH(ecef)
    % ECEF to BLH
    % WGS84椭球参数
    a = 6378137;
    e2 = 0.0066943799013;
    
    X=ecef(1);
    Y=ecef(2);
    Z=ecef(3);
    L = atan2(Y, X);
    B = 999;
    B_prime = atan2(Z, sqrt(X^2 + Y^2));
    while abs(B - B_prime) > 1e-8
        B = B_prime;
        W = sqrt(1 - e2 * sin(B)^2);
        N = a / W;
        H = sqrt(X^2 + Y^2) / cos(B) - N;
        B_prime = atan2(Z * (N + H), sqrt(X^2 + Y^2) * (N * (1 - e2) + H));
    end

end
