clear; clc; close all;

%% ------------------ ПАРАМЕТРЫ ------------------
A = 1; % амплитуда

% Гармонический сигнал
f_car_harm = 500;
T_harm = 0.12;
fs = 2e8;
dt = 1/fs;

% Видеоимпульс
T_video = 0.02;

% Пачка импульсов
T_imp = 0.04;               
N_rep = 3;
T_pulse_pack = N_rep*T_imp;

% Радиоимпульс
f_car_radio = 500;

% M-последовательность
f_car_m = 100e3;
T_m = 0.02;
fs_m = 2e8; 
dt_m = 1/fs_m;
n_M = 511;
sym_T = T_m / n_M;

% Код Баркера 13
barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
f_car_b = 20000;
sym_T_b = T_video / 13;

% ЛЧМ-сигнал
f0 = 1000;
f1 = 2000;
T_lfm = 0.02;

%% ------------------ ВРЕМЕННЫЕ МАССИВЫ ------------------
t_harm = 0:dt:T_harm;
t_video = -T_video:dt:T_video;
t_pack = -T_pulse_pack:dt:T_pulse_pack;
t_m = 0:dt_m:T_m;
t_barker = 0:dt:T_video;
t_lfm = 0:dt:T_video;

%% ------------------ СИГНАЛЫ ------------------
% Гармонический
harm = A * sin(2*pi*f_car_harm*t_harm);

% Видеоимпульс
video = double(abs(t_video) <= T_video/2);

% Пачка видеоимпульсов
starts = (-floor(N_rep/2)*T_imp) : T_imp : (floor(N_rep/2)*T_imp);
video_pack = zeros(size(t_pack));
for k = 1:length(starts)
    video_pack = video_pack + double(abs(t_pack - starts(k)) <= T_video/2);
end

% Радиоимпульс
radio = video .* sin(2*pi*f_car_radio*t_video);

% Пачка радиоимпульсов
radio_pack = zeros(size(t_pack));
for k = 1:length(starts)
    radio_pack = radio_pack + ...
        double(abs(t_pack - starts(k)) <= T_video/2) .* sin(2*pi*f_car_radio*(t_pack - starts(k)));
end

% М-последовательность
m_seq = 2*(round(rand(1,n_M))-0.5);  
m_signal = zeros(size(t_m));
for i=1:n_M
    idx = t_m >= (i-1)*sym_T & t_m < i*sym_T;
    m_signal(idx) = m_seq(i);
end
m_signal = m_signal .* sin(2*pi*f_car_m*t_m);

% Код Баркера
barker_signal = zeros(size(t_barker));
for i=1:13
    idx = t_barker >= (i-1)*sym_T_b & t_barker < i*sym_T_b;
    barker_signal(idx) = barker(i);
end
barker_signal = barker_signal .* sin(2*pi*f_car_b*t_barker);

% ЛЧМ
k = (f1 - f0) / T_lfm;
lfm = A * sin(2*pi*(f0*t_lfm + 0.5*k*t_lfm.^2));

%% ------------------ АВТОКОРРЕЛЯЦИОННЫЕ ФУНКЦИИ ------------------
acf = @(x) xcorr(x, 'coeff');

[acf_harm, l_harm] = acf(harm);
[acf_video, l_video] = acf(video);
[acf_video_pack, l_vp] = acf(video_pack);
[acf_radio, l_radio] = acf(radio);
[acf_radio_pack, l_rp] = acf(radio_pack);
[acf_m, l_m] = acf(m_signal);
[acf_barker, l_b] = acf(barker_signal);
[acf_lfm, l_lfm] = acf(lfm);

%% ------------------ ОТОБРАЖЕНИЕ (НОРМИРОВАНИЕ ОСЕЙ) ------------------

% Гармонический сигнал
figure;
plot(l_harm/fs, acf_harm/max(acf_harm)); grid on; 
title('Нормированная АКФ гармонического сигнала');
xlabel('Нормированное время задержки'); ylabel('Нормированная АКФ');
xlim([-T_harm T_harm]/T_harm);

% Видеоимпульс + пачка видеоимпульсов
figure;
subplot(2,1,1)
plot(l_video/fs, acf_video/max(acf_video), 'LineWidth', 1.7);
grid on; title('Нормированная АКФ одиночного видеоимпульса');
xlabel('Нормированное время задержки'); ylabel('Нормированная АКФ');
xlim([-T_video T_video]/T_video);

subplot(2,1,2)
plot(l_vp/fs, acf_video_pack/max(acf_video_pack), 'LineWidth', 1.7);
grid on; title('Нормированная АКФ пачки видеоимпульсов');
xlabel('Нормированное время задержки'); ylabel('Нормированная АКФ');
xlim([-T_pulse_pack T_pulse_pack]/T_pulse_pack);

% Радиоимпульс + пачка радиоимпульсов
figure;
subplot(2,1,1)
plot(l_radio/fs, acf_radio/max(acf_radio), 'LineWidth', 1.7);
grid on; title('Нормированная АКФ одиночного радиоимпульса');
xlabel('Нормированное время задержки'); ylabel('Нормированная АКФ');
xlim([-T_video T_video]/T_video);

subplot(2,1,2)
plot(l_rp/fs, acf_radio_pack/max(acf_radio_pack), 'LineWidth', 1.7);
grid on; title('Нормированная АКФ пачки радиоимпульсов');
xlabel('Нормированное время задержки'); ylabel('Нормированная АКФ');
xlim([-T_pulse_pack T_pulse_pack]/T_pulse_pack);

% М-последовательность
figure; 
plot(l_m/fs_m, acf_m/max(acf_m)); grid on; 
title('Нормированная АКФ М-последовательности');
xlabel('Нормированное время'); ylabel('Нормированная АКФ');
xlim([0 T_m]/T_m);

% Код Баркера
figure; 
plot(l_b/fs, acf_barker/max(acf_barker)); grid on; 
title('Нормированная АКФ кода Баркера 13');
xlabel('Нормированное время'); ylabel('Нормированная АКФ');
xlim([0 T_video]/T_video);

% ЛЧМ
figure; 
plot(l_lfm/fs, acf_lfm/max(acf_lfm)); grid on; 
title('Нормированная АКФ ЛЧМ-сигнала');
xlabel('Нормированное время'); ylabel('Нормированная АКФ');
xlim([0 T_lfm]/T_lfm);
