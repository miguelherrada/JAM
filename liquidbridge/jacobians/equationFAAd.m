function [FAAd,dFAAd] = equationFAAd(s,y,in3,in4)
%equationFAAd
%    [FAAd,dFAAd] = equationFAAd(S,Y,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    20-Dec-2023 16:12:40

dfdy = in3(:,25);
dpdy = in3(:,18);
dwdy = in3(:,4);
f = in3(:,22);
u = in3(:,8);
t2 = 1.0./f;
FAAd = [dwdy.*t2,u,dpdy.*t2,dfdy];
if nargout > 1
    t3 = t2.^2;
    dFAAd = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-dwdy.*t3,0.0,-dpdy.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,28]);
end
