%Example 3.9 - Digital Integrator
%Based on Lathi_EDO

fs = 66.667*1e3; %frequÊncia de amostragem
Ts = 1/fs; %período de amostragem
t_final = 1e-3/6;
t = linspace(0,t_final,8000); %índices contínuos
x = t; %x(t)
desloc = 1; %deslocamento 

n_disc = -10:1:round(t_final*fs); %índices discretos
x_n = 0:1:round(t_final*fs); %x[n] p/ n>=0
impulso = (n_disc==0); %impulso unitário discreto no tempo
x_disc = conv(x_n,impulso); %x[n]
n_conv = (x_n(1)+n_disc(1):1:(x_n(length(x_n))+n_disc(length(n_disc))));

%entrada
y_diff = zeros(1,length(x_disc));
idx_max = find(x_disc==max(x_disc));
for i=0:1:n_conv(idx_max)
    idx = find(n_conv==i);
    y_diff(1,idx) = fs*(x_disc(1,idx)-x_disc(1,idx-desloc));
end

%integrador
y_int = zeros(1,length(y_diff));
idx_max2 = find(y_diff==max(y_diff));
for j=0:1:max(n_conv(idx_max2))
    idx2 = find(n_conv==j);
    y_int(1,idx2) = Ts*y_diff(1,idx2) + y_int(1,idx2-1);
end

subplot(2,1,1);
stem(n_conv,y_diff,'filled');
hold on
plot(n_conv,y_diff);
hold off
axis([0 max(n_conv(idx_max2)) 0 fs+1e4]);
xlabel('n');
ylabel('Amplitude')
title('Entrada')
grid on

subplot(2,1,2);
stem(n_conv,y_int,'filled');
hold on
plot(n_conv,y_int);
hold off
axis([0 n_conv(idx_max) 0 12]);
xlabel('n');
ylabel('Amplitude')
title('Saída')
grid on



