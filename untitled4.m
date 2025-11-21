clear all
clc
tau=0.02; %длительность импульса
n_p=10; %количество периодов несущей в импульсе
q=2; %скважность импульса (период/длит. импульса)
Tp=q*tau; %период импульса
f_nes=n_p/tau; %частота несущей
fs=500*f_nes; %частота дискретизации
dt=1/fs; %шаг дискретизации
n=3; %количество повторений импульса
T=0:dt:n*Tp; %длительность сигнала
T1=0:Tp:n*Tp; %повторения



mod=rectpuls(T,2*tau);%видеоимпульс ч/з ректпульс



nes=sin(2*pi*f_nes.*T);%несущая
videoimp=rectpuls(T,2*tau);%видеоимпульс
rad=videoimp.*nes; %радиоимпульс
pack=pulstran(T,T1,rad,fs); %пачка
pack_video = pulstran(T,T1,videoimp,fs);

%-----------------
%ЛЧМ
Fmin = 4000;
Fmax = 5000;
T_lfm = 0.02; %Длительность одиночного импульса
Fs_lfm = 400000;
dt_lfm = 1/Fs_lfm;
t_lfm = 0:dt_lfm:T_lfm;
deviation = Fmax - Fmin; 
u = deviation/T_lfm; %veloc изм частоты
lfm = cos(2*pi*(Fmin+u./2*t_lfm).*t_lfm); %Несущая

%ФКМ
tau_lfm=3.94*(10^-5);%0.001; %длительность символа
n_M=3; %количество периодов в символе
f_nes_M=n_M/tau_lfm; %частота несущей
Tp_M = 511*tau_lfm;
fd_M=300000; %частота дискретизации для М-посл
dt_M=1/fd_M; %шаг дискретизации для М-посл
t_M=0:dt_M:Tp_M; %временные отсчеты
imp_M=ones(1,length(t_M));

M = function_test;

for k=1:length(t_M)
imp_M(1,k)=M(1,ceil((k)/(length(t_M))*length(M))); %подгонка временных отсчетов под код Баркера
end

nes_M=sin(2*pi*f_nes_M.*t_M);%несущая
mposl=nes_M.*imp_M;

plot(T,nes);
xlabel('t,c');
ylabel('s(t)');
title("Гармонический");
grid on;
hold on;

plot(T,videoimp);
xlabel('t,c');
ylabel('s(t)');
title("Видеоимпульс");
grid on;


plot(T,rad);
xlabel('t,c');
ylabel('s(t)');
title("Радиоимпульс");
grid on;
legend('Гармонический', "Видео", "Радио");
hold off;



plot(T,pack);
xlabel('t,c');
ylabel('s(t)');
title("Пачка");
grid on;
hold on;
plot(T,pack_video);
legend('Радиоимпульсов', "Видеоимпульсов")
hold off;





plot(t_M,nes_M);
xlabel('t,c');
ylabel('s(t)');
title("Гармонический");
axis([0 0.02 -1.5 1.5]);
grid on;
hold off;

plot(t_M,imp_M);
xlabel('t,c');
ylabel('s(t)');
title("М-Видеоимпульс");
grid on;

plot(t_M,mposl);
xlabel('t,c');
ylabel('s(t)');
title("ФКМ импульс (М-послед.)");
grid on;

%legend("Гармонический", "М-видеоимпульс", "ФКМ")






plot(t_lfm, lfm);
xlabel('t,c');
ylabel('s(t)');
title("ЛЧМ");
grid on;




%Баркер
tau_bark=0.0015385;%=0.001; %длительность символа
n_B=15; 
f0_B=n_B/tau_bark; %частота несущей
Tp_B = 13*tau_bark; %период для пачки
fd_B=500000; %частота дискретизации
dt_B=1/fd_B; %шаг дискретизации
t_B=0:dt_B:Tp_B; %временные отсчеты
barker=[1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
imp=ones(1,length(t_B));

for i=1:length(t_B)
imp(1,i)=barker(1,ceil((i)/(length(t_B))*length(barker))); %подгонка временных отсчетов под код Баркера
end

nes_bark=sin(2*pi*f0_B.*t_B);%несущая
bark=nes_bark.*imp; 




plot(t_B,nes_bark);
xlabel('t,c');
ylabel('s(t)');
title("Гармонический");
axis([0 0.02 -1.5 1.5]);
grid on;


plot(t_B,imp);
xlabel('t,c');
ylabel('s(t)');
title("Код Баркера");

plot(t_B,bark);
xlabel('t,c');
ylabel('s(t)');
title("ФКМ импульс (Код Баркера)");
axis([0 0.02 -1.1 1.1])
grid on;

%legend("Гармонический", "Код баркера", "ФКМ")


hold off


%% Амплитудные спектры 
%Несущая
N=2^nextpow2(length(nes));
F_range=(0:N/2-1)*(fs/N);
spec_nes=fft(nes,N);
amp_nes=abs(spec_nes)./N;
%Видео
spec_mod=fft(mod,N);
mod0=fftshift(spec_mod);
F_range0=(-N/2:N/2-1)*(fs/N);
amp_mod=abs(mod0)./N;
%радио
spec_rad=fft(rad,N);
amp_rad=abs(spec_rad)./N;
%пачка
% радио
spec_pack=fft(pack,N);
amp_pack=abs(spec_pack)./N;
% видео
spec_pack_video=fft(pack_video,N);
amp_pack_video=abs(spec_pack_video)./N;


%лчм
N_lfm = 2^nextpow2(length(lfm));
F_range_lfm = (0:N_lfm/2-1)*(Fs_lfm/N_lfm);
spec_lfm = fft(lfm, N_lfm);
amp_lfm = abs(spec_lfm)./length(N_lfm);



% м-послед.
N_M = 2^nextpow2(length(mposl));
F_range_M = (0:N_M/2-1)*(fd_M/N_M);
spec_M = fft(mposl, N_M);
amp_M = abs(spec_M)./N_M;
%баркер
N_bark = 2^nextpow2(length(bark));
F_range_bark = (0:N_bark/2-1)*(fd_B/N_bark);
spec_bark = fft(bark, N_bark);
amp_bark = abs(spec_bark)/N_bark;


% отрисовка

plot(F_range,amp_nes(1:N/2),'LineWidth',2,'Color','blue')
title('Спектр гармонического сигнала');
xlabel('f, Гц');
ylabel('A');
axis([0 f_nes*2 -inf inf]);
grid on;

%disp("Эф. ширина спектра")



plot(F_range0,amp_mod,'LineWidth',2,'Color','blue')
title('спектр видеоимпульса');
xlabel('f, Гц');
ylabel('A');
axis([-5000 5000 0 inf]);
grid on;
disp("Эф. ширина спектра")
1/tau


plot(F_range,amp_rad(1:N/2),'LineWidth',2,'Color','blue')
title('Спектр радиоимпульса');
xlabel('f, Гц');
ylabel('A');
axis([0 f_nes*2 0 inf]);
grid on;
disp("Эф. ширина спектра")
2/tau


plot(F_range,amp_pack(1:N/2),'LineWidth',2,'Color','blue')
title('Спектр пачки радиоимпульсов');
xlabel('f, Гц');
ylabel('A');
axis([0 f_nes*2 0 inf]);
grid on;
disp("Эф. ширина спектра")
1/(tau)

plot(F_range,amp_pack_video(1:N/2),'LineWidth',2,'Color','blue')
title('Спектр пачки видеоимпульсов');
xlabel('f, Гц');
ylabel('A');
axis([0 inf 0 inf]);
grid on;
disp("Эф. ширина спектра")
2/(tau)

plot(F_range_lfm, amp_lfm(1:N_lfm/2), 'LineWidth', 1.5);
title('Спектр ЛЧМ'), grid on;
xlabel('Частота f, Гц'), ylabel('Амплитуда');
axis([3000 6000 0 inf]);
disp("Эф. ширина спектра")
Fmax-Fmin

plot(F_range_M,amp_M(1:N_M/2),'LineWidth',2,'Color','blue');
title ('Спектр М-последовательности');
axis([0 f_nes_M*2 0 inf]);
grid on;
disp("Эф. ширина спектра")
1/tau_lfm

plot(F_range_bark,amp_bark(1:N_bark/2),'LineWidth',2,'Color','blue');
title ('Спектр кода Баркера');
axis([0 f0_B*2 0 inf]);
grid on;
disp("Эф. ширина спектра")
2/tau_bark


%----------------------
% Фазовые спектры
% несущая
ph_nes=angle(spec_nes);
%видео
ph_mod=angle(mod0);
%радио
ph_rad=angle(spec_rad);
%пачка
% радио
ph_pack=angle(spec_pack);
% видео 
ph_pack_video=angle(spec_pack_video);
%лчм
ph_lfm=angle(spec_lfm);
%м
ph_M=angle(spec_M);
%баркер
ph_bark=angle(spec_bark);
%отрисовка
plot(F_range,ph_nes(1:N/2),'LineWidth',2,'Color','blue')
title('Фазовый спектр несущей');
xlabel('f, Гц');
ylabel('fi');
axis([0 f_nes -inf inf]);
grid on;


plot(F_range0,ph_mod,'LineWidth',2,'Color','blue')
title('Фазовый спектр видео');
xlabel('f, Гц');
ylabel('fi');
axis([-1000 1000 -inf inf]);
grid on;


plot(F_range,ph_rad(1:N/2),'LineWidth',2,'Color','blue')
title('Фазовый спектр радиоимпульсов');
xlabel('f, Гц');
ylabel('fi');
axis([0 f_nes -inf inf]);
grid on;

plot(F_range0,ph_pack_video,'LineWidth',2,'Color','blue')
title('Фазовый спектр пачки видеоимпульсов');
xlabel('f, Гц');
ylabel('A');
axis([0 f_nes -inf inf]);
grid on;

plot(F_range,ph_pack(1:N/2),'LineWidth',2,'Color','blue')
title('Фазовый спектр пачки радиоимпульсов');
xlabel('f, Гц');
ylabel('A');
axis([0 f_nes -inf inf]);
grid on;


stem(F_range_lfm, ph_lfm(1:N_lfm/2), 'LineWidth', 1.5);
title('Фазовый спектр ЛЧМ'), grid on;
xlabel('Частота f, Гц'), ylabel('Амплитуда');
axis([0 Fmax*2 -inf inf]);

plot(F_range_M,ph_M(1:N_M/2),'Color','blue');
title ('Фазовый спектр М-последовательности');
axis([0 f_nes_M*2 -inf inf]);
grid on;

stem(F_range_bark,ph_bark(1:N_bark/2),'LineWidth',1,'Color','blue');
title ('Фазовый спектр Баркера');
axis([0 f0_B*2 -inf inf]);
grid on;





function ans=get_ef(x,y)
    ind_zero = find(x==0);
    y_real = y(:,ind_zero:end);
    extremum = max(y_real);
    ans=sum(y_real)/extremum;
end

% --------------------
% video_xcorr = xcorr(videoimp);
% radio_xcorr = xcorr(rad);
% m_xcorr = xcorr(imp_M);
% lfm_xcorr = xcorr(lfm);
% bark_xcorr = xcorr(bark);
% 
% plot(video_xcorr)
% title("Видео XCORR");
% grid on;
% 
% plot(radio_xcorr)
% title("Радиоимпульс XCORR");
% grid on;
% 
% 
% plot(m_xcorr)
% title("М-Видеоимпульс Xcorr");
% grid on;
% 
% m_posl_xcorr =xcorr(mposl);
% plot(m_posl_xcorr)
% title("ФКМ импульс Xcorr");
% grid on;
% 
% 
% plot(lfm_xcorr);
% title("ЛЧМ XCORR");
% grid on;
% 
% plot(bark_xcorr);
% title("Cod Barker (видеоимп) Xcorr");
% grid on;
