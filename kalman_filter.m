clear;
close all;
clc;

x = csvread('gnss_data.csv');
Ts = 1;

sigma_w = zeros(1,length(x));
for i = 1:length(sigma_w)
    if(i < 40)
        sigma_w(i) = 10;
    elseif(i == 40)
        sigma_w(i) = 20;
    elseif(i == 41)
        sigma_w(i) = 30;
    elseif(i == 42)
        sigma_w(i) = 40;
    elseif(i == 43)
        sigma_w(i) = 50;
    elseif(i == 44)
        sigma_w(i) = 60;
    elseif(i <= 55)
        sigma_w(i) = inf;
    elseif(i == 56)
        sigma_w(i) = 60;
    elseif(i == 57)
        sigma_w(i) = 50;
    elseif(i == 58)
        sigma_w(i) = 40;
    elseif(i == 59)
        sigma_w(i) = 30;
    elseif(i == 60)
        sigma_w(i) = 20;
    else
        sigma_w(i) = 10;
    end
end

mi_a = zeros(1,length(x));
for i = 1:length(mi_a)
   if(i <= 30)
       mi_a(i) = 5;
   elseif(i > 70)
       mi_a(i) = -5;
   else
       mi_a(i) = 0;
   end
end

C = sigma_w.^2;

sigma_a = zeros(1,length(x));
for i = 1:length(sigma_a)
    if(i <=30)
        sigma_a(i) = 5/3;
    elseif(i > 70)
        sigma_a(i) = 5/3;
    else
        sigma_a(i) = 1/3;
    end
end

Q =@(i) [0 0 0;0 0 0;0 0 sigma_a(i)^2];

A = [1 Ts Ts^2/2;0 1 Ts;0 0 0];
B = [0;0;1];
H = [1 0 0];

s_est = zeros(3,1);
M_est = zeros(3,3);

s_estimirano = zeros(3,length(x)+1);
M_estimirano = zeros(3,length(x)+1);
K_zap = zeros(3,length(x));

%Inicijalizacija
s_estimirano(:,1) = s_est;
M_estimirano(:,1) = [M_est(1,1) M_est(2,2) M_est(3,3)];

%Iteracije Kalmanovog filtra
for i = 1:length(x)
    %predikcija
    s_pred = A*s_est + B*mi_a(i);
    M_pred = A*M_est*A' + B'*Q(i)*B;
    
    %estimacija
    if(C(i) ~= inf)
        K = M_pred*H'*(C(i) + H*M_pred*H')^-1;
        s_est = s_pred + K*(x(i) - H*s_pred);
        M_est = (eye(3) - K*H)*M_pred;
    else
        K = 0;
        s_est = s_pred;
        M_est = M_pred;
    end
    %pamcenje
    s_estimirano(:,i+1) = s_est;
    M_estimirano(:,i+1) = [M_est(1,1) M_est(2,2) M_est(3,3)];
    K_zap(:,i) = K;
end

t = 1:Ts:length(x)*Ts;

figure(1)
plot(t,x);
hold on;
plot(t,s_estimirano(1,2:end));
hold on;
plot(t,s_estimirano(1,2:end)+2*sqrt(M_estimirano(1,2:end)),'g--');
hold on;
plot(t,s_estimirano(1,2:end)-2*sqrt(M_estimirano(1,2:end)),'g--');
hold off;
legend('Mereno','Estimirano','2\sigma interval poverenja');
xlabel('t[s]');
ylabel('$$\hat{p}(t)$$','Interpreter','Latex')
title('Pozicija');

figure(2)
plot(t,s_estimirano(2,2:end));
hold on;
plot(t,s_estimirano(2,2:end)+2*sqrt(M_estimirano(2,2:end)),'g--');
hold on;
plot(t,s_estimirano(2,2:end)-2*sqrt(M_estimirano(2,2:end)),'g--');
hold off;
xlabel('t[s]');
ylabel('$$\hat{v}(t)$$','Interpreter','Latex');
title('Brzina');
legend('Estimacija','2\sigma interval poverenja');

figure(3)
plot(t,s_estimirano(3,2:end));
hold on;
plot(t,s_estimirano(3,2:end)+2*sqrt(M_estimirano(3,2:end)),'g--');
hold on;
plot(t,s_estimirano(3,2:end)-2*sqrt(M_estimirano(3,2:end)),'g--');
hold off;
xlabel('t[s]');
ylabel('$$\hat{a}(t)$$','Interpreter','Latex');
title('Ubrzanje');
legend('Estimacija','2\sigma interval poverenja');

figure(4)
plot(t,K_zap(1,:));
hold on;
plot(t,K_zap(2,:));
hold on;
plot(t,K_zap(3,:));
hold off;
xlabel('t[s]');
title('Kalmanovo pojacanje');
legend('Za poziciju','Za brzinu','Za ubrzanje');
