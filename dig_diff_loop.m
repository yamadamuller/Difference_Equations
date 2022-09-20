%Example 3.8 - Digital Differentiator
%Based on Lathi_EDO

fs = 66.667*1e3; %frequÊncia de amostragem
t_final = 1e-3/6;
t = linspace(0,t_final,8000); %índices contínuos
x = t; %x(t)
desloc = 1; %deslocamento 

n_disc = -10:1:round(t_final*fs); %índices discretos
x_n = 0:1:round(t_final*fs); %x[n] p/ n>=0
impulso = (n_disc==0); %impulso unitário discreto no tempo
x_disc = conv(x_n,impulso); %x[n]
n_conv = (x_n(1)+n_disc(1):1:(x_n(length(x_n))+n_disc(length(n_disc))));
%índices da convolução

y_diff = zeros(1,length(x_disc));
idx_max = find(x_disc==max(x_disc))
for i=0:1:n_conv(idx_max)
    idx = find(n_conv==i);
    y_diff(1,idx) = fs*(x_disc(1,idx)-x_disc(1,idx-desloc));
end

subplot(2,1,1);
stem(n_conv,x_disc,'filled');
hold on
plot(t.*fs,x.*fs);
hold off
axis([0 n_conv(idx_max) 0 x_disc(idx_max)]);
xlabel('n');
ylabel('Amplitude')
title('Entrada')
grid on

subplot(2,1,2);
stem(n_conv,y_diff,'filled');
hold on
plot(n_conv,y_diff);
hold off
axis([0 n_conv(idx_max) 0 fs+1e4]);
xlabel('n');
ylabel('Amplitude')
title('Saída')
grid on