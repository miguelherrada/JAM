function [FAAb,dFAAb] = equationFAAb(s,in2,in3)
%equationFAAb
%    [FAAb,dFAAb] = equationFAAb(S,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    27-Nov-2023 13:27:28

dgds = in2(:,6);
dgdss = in2(:,7);
dgdtau = in2(:,8);
duds = in2(:,2);
dudss = in2(:,3);
dudtau = in2(:,4);
pa_1 = in3(1,:);
pa_2 = in3(2,:);
u = in2(:,1);
t2 = dgds.^2;
t3 = duds.^2;
t4 = duds.^3;
t5 = pa_2.^2;
t6 = 1.0./dgds;
t7 = 1.0./t2;
t8 = t6.^3;
t10 = t6.^5;
t11 = duds.*t6;
t12 = dudss.*t6;
t13 = t2.^(3.0./2.0);
t9 = t7.^2;
t14 = dgdss.*duds.*t7;
t15 = pa_2.*t3.*t7;
t16 = -t14;
t17 = t15+2.0./5.0;
t18 = 1.0./t17.^2;
t19 = 1.0./t17.^3;
t20 = t12+t16;
t21 = duds.*dudss.*pa_2.*t7.*t18.*2.0;
t22 = dgdss.*pa_2.*t3.*t8.*t18.*2.0;
t23 = -t22;
t24 = t21+t23;
FAAb = [dudtau-dgdtau.*t11+t11.*u-pa_1.*t6.*t20,dgds.*dgdss+t13.*t24];
if nargout > 1
    dFAAb = reshape([t11,0.0,-dgdtau.*t6+t6.*u+dgdss.*pa_1.*t8,t13.*(dudss.*pa_2.*t7.*t18.*2.0-dgdss.*duds.*pa_2.*t8.*t18.*4.0+dgdss.*t4.*t5.*t10.*t19.*8.0-dudss.*t3.*t5.*t9.*t19.*8.0),-pa_1.*t7,duds.*pa_2.*sqrt(t2).*t18.*2.0,1.0,0.0,0.0,0.0,pa_1.*t6.*(dudss.*t7-dgdss.*duds.*t8.*2.0)+dgdtau.*duds.*t7-duds.*t7.*u+pa_1.*t7.*t20,dgdss-t13.*(dgdss.*t3.^2.*t5.*t7.^3.*t19.*8.0+duds.*dudss.*pa_2.*t8.*t18.*4.0-dgdss.*pa_2.*t3.*t9.*t18.*6.0-dudss.*t4.*t5.*t10.*t19.*8.0)+dgds.*sqrt(t2).*t24.*3.0,duds.*pa_1.*t8,dgds-pa_2.*t3.*t8.*t13.*t18.*2.0,-t11,0.0],[2,8]);
end
end