function [FAAt,dFAAt] = equationFAAt(s,y,in3,in4)
%equationFAAt
%    [FAAt,dFAAt] = equationFAAt(S,Y,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    20-Dec-2023 16:12:40

dfds = in3(:,23);
dfdt0 = in3(:,28);
dfdss = in3(:,24);
dpdy = in3(:,18);
duds = in3(:,9);
dudy = in3(:,11);
dwds = in3(:,2);
dwdy = in3(:,4);
f = in3(:,22);
p = in3(:,15);
pa_1 = in4(1,:);
pa_2 = in4(2,:);
u = in3(:,8);
w = in3(:,1);
t2 = dfdss.*f;
t3 = dfds.^2;
t4 = dfds.^3;
t5 = dwds.*2.0;
t6 = 1.0./f;
t7 = t6.^2;
t8 = t3+1.0;
t9 = -t2;
t10 = dwdy.*t6;
t12 = dfds.*dudy.*t6.*y;
t11 = dwdy.*t7;
t13 = dfds.*dudy.*t7.*y;
t14 = 1.0./t8;
t16 = dfds.*t10.*y.*2.0;
t17 = t8+t9;
t18 = -t12;
t21 = 1.0./t8.^(3.0./2.0);
t15 = t14.^2;
t19 = -t13;
t20 = -t16;
t22 = dfds.*pa_1.*t14.*2.0;
t24 = dfds.*t6.*t14.*2.0;
t25 = dudy.*t6.*t14.*2.0;
t27 = duds+t10+t18;
t28 = t3.*t6.*t14.*y.*2.0;
t23 = t5+t20;
t26 = t11+t19;
t29 = dfds.*t14.*t27.*2.0;
FAAt = [p-pa_1.*(t25-t29+t3.*t14.*t23)+pa_2.*s+t6.*t21.*(t2-t8),-pa_1.*(dudy.*t24+t14.*t27-dfds.*t14.*t23-t3.*t14.*t27),dpdy.*t6,dfdt0-u+dfds.*w];
if nargout > 1
    mt1 = [0.0,0.0,0.0,dfds,pa_1.*t3.*t14.*-2.0,t22,0.0,0.0,0.0,0.0,0.0,0.0,pa_1.*(t24+t4.*t6.*t14.*y.*2.0),-pa_1.*(t28+t6.*t14-t3.*t6.*t14),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t22,-pa_1.*(t14-t3.*t14),0.0,0.0,0.0,0.0,0.0,0.0,-pa_1.*(t28+t6.*t14.*2.0),-pa_1.*(t24-dfds.*t6.*t14.*y+t4.*t6.*t14.*y),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-pa_1.*(dfds.*t14.*t26.*2.0-dudy.*t7.*t14.*2.0+t4.*t11.*t14.*y.*2.0)+dfdss.*t6.*t21-t7.*t21.*(t2-t8)];
    mt2 = [pa_1.*(t14.*t26-t3.*t14.*t26+dfds.*dudy.*t7.*t14.*2.0+t3.*t11.*t14.*y.*2.0),-dpdy.*t7,0.0,pa_1.*(t12.*t14.*-2.0+t14.*t27.*2.0-dfds.*t14.*t23.*2.0+t4.*t15.*t23.*2.0-t3.*t15.*t27.*4.0+dfds.*dudy.*t6.*t15.*4.0+t3.*t10.*t14.*y.*2.0)-dfds.*t6.*t21.*2.0-dfds.*t6.*1.0./t8.^(5.0./2.0).*(t2-t8).*3.0,-pa_1.*(t25-t29-t14.*t23-dfds.*t15.*t27.*2.0+t3.*t15.*t23.*2.0+t4.*t15.*t27.*2.0-dudy.*t3.*t6.*t15.*4.0+dfds.*t10.*t14.*y.*2.0-dudy.*t6.*t14.*y+dudy.*t3.*t6.*t14.*y),0.0,w,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0];
    dFAAt = reshape([mt1,mt2],4,28);
end
