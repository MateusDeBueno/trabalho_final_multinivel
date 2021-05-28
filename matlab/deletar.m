close all
clear all
clc
f_rede = 60; %freq da rede
Vref = 3.4; %a tensao de referencia maxima eh 3.4607
fs = 7200; %frequencia de amostragem (7 segmentos)
n_pontos_t = fs/f_rede;
n_pontos = n_pontos_t+1;
time = linspace(1/(2*fs),(1/f_rede+1/(2*fs)),n_pontos); %vetor de tempo
for k=1:length(time)
    VA(k) = Vref*sin(time(k)*2*pi*f_rede);
    VB(k) = Vref*sin(time(k)*2*pi*f_rede+2*pi/3);
    VC(k) = Vref*sin(time(k)*2*pi*f_rede+4*pi/3);
end
figure
hold on
plot(VA)
plot(VB)
plot(VC)
hold off