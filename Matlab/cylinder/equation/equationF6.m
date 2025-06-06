function [FF6,dFF6] = equationF6(z0,r0,in3,in4)
%equationF6
%    [FF6,dFF6] = equationF6(Z0,R0,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    16-May-2025 19:37:25

pa_5 = in4(5,:);
yb1v1d1 = in3(1,:);
yb1v2d2 = in3(9,:);
yb1v2d3 = in3(10,:);
yb1v3d2 = in3(16,:);
yb1v3d3 = in3(17,:);
yb1v4d1 = in3(22,:);
yb1v4d2 = in3(23,:);
yb1v4d3 = in3(24,:);
yb1v5d2 = in3(30,:);
yb1v5d3 = in3(31,:);
yb1v5d5 = in3(33,:);
yb2v1d1 = in3(36,:);
yb2v2d1 = in3(43,:);
yb2v3d1 = in3(50,:);
yb2v4d1 = in3(57,:);
yb2v5d1 = in3(64,:);
t2 = yb1v4d2.*yb1v5d3;
t3 = yb1v4d3.*yb1v5d2;
t4 = yb1v5d2.^2;
t5 = yb1v5d3.^2;
t6 = -t3;
t7 = t2+t6;
t8 = 1.0./t7;
t9 = t8.^2;
t10 = t8.*yb1v5d2;
t11 = t8.*yb1v5d3;
t12 = -t10;
FF6 = [yb1v1d1;t11.*yb1v2d2+t12.*yb1v2d3;t11.*yb1v3d2+t12.*yb1v3d3;-pa_5+yb1v4d1;yb1v5d5;yb2v1d1;yb2v2d1;yb2v3d1;yb2v4d1;yb2v5d1];
if nargout > 1
    mt1 = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5.*t9.*yb1v2d2+t9.*yb1v2d3.*yb1v5d2.*yb1v5d3,-t5.*t9.*yb1v3d2+t9.*yb1v3d3.*yb1v5d2.*yb1v5d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t4.*t9.*yb1v2d3+t9.*yb1v2d2.*yb1v5d2.*yb1v5d3,-t4.*t9.*yb1v3d3+t9.*yb1v3d2.*yb1v5d2.*yb1v5d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t8.*yb1v2d3+t6.*t9.*yb1v2d3+t9.*yb1v2d2.*yb1v4d3.*yb1v5d3,-t8.*yb1v3d3+t6.*t9.*yb1v3d3+t9.*yb1v3d2.*yb1v4d3.*yb1v5d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t8.*yb1v2d2-t2.*t9.*yb1v2d2+t9.*yb1v2d3.*yb1v4d2.*yb1v5d2];
    mt3 = [t8.*yb1v3d2-t2.*t9.*yb1v3d2+t9.*yb1v3d3.*yb1v4d2.*yb1v5d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt4 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt5 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    dFF6 = reshape([mt1,mt2,mt3,mt4,mt5],10,70);
end
end
