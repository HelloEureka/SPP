function enu = ECEF2ENU(d_ecef, B_0, L_0)
    % ECEF坐标系下与站心差量在ENU下投影
    % d_ecef 列向量
    A = [-sin(L_0), cos(L_0), 0;
        sin(B_0) * cos(L_0), -sin(B_0) * sin(L_0), cos(B_0);
        cos(B_0) * cos(L_0), cos(B_0) * sin(L_0), sin(B_0)];
    enu = A * d_ecef;
end
