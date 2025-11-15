clear; clc; close all;

%% ------------------ ПАРАМЕТРЫ ИЗ ТАБЛИЦЫ ------------------
A = 1;

% Гармонический сигнал
f_car_harm = 500;
T_harm = 0.12;
fs = 2e8;                     % Частота дискретизации (общая)
dt = 1/fs;

% Видеоимпульс
T_video = 0.02;

% Пачка импульсов
T_imp = 0.04;                 % Период (скважность q=2)
N_rep = 3;
T_pulse_pack = N_rep*T_imp;   % 0.12

% Радиоимпульс
f_car_radio = 500;

% M-последовательность
f_car_m = 100e3;
T_m = 0.02;
fs_m = 2e6; 
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
t_video = 0:dt:T_video;
t_pack = 0:dt:T_pulse_pack;
t_m = 0:dt_m:T_m;
t_barker = t_video;
t_lfm = t_video;

%% ------------------ СИГНАЛЫ ------------------

% Гармонический
harm = A * sin(2*pi*f_car_harm*t_harm);

% Видеоимпульс
video = A * rectpuls(t_video - T_video/2, T_video);

% Пачка видеоимпульсов
starts = 0:T_imp:T_pulse_pack;
video_pack = pulstran(t_pack, starts', @(t) rectpuls(t, T_video));

% Радиоимпульс
radio = video .* sin(2*pi*f_car_radio*t_video);

% Пачка радиоимпульсов
radio_pack = pulstran(t_pack, starts', @(t) rectpuls(t, T_video).*sin(2*pi*f_car_radio*t));

% М-последовательность
m_seq = 2*(round(rand(1,n_M))-0.5);     % псевдо М-код (±1)
m_signal = zeros(size(t_m));
for i=1:n_M
    idx = t_m >= (i-1)*sym_T & t_m < i*sym_T;
    m_signal(idx) = m_seq(i);
end
m_signal = m_signal .* sin(2*pi*f_car_m*t_m);

% Баркера 13
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

%% ------------------ ОТОБРАЖЕНИЕ ------------------

figure; plot(l_harm/fs, acf_harm); grid on; title('АКФ гармонического сигнала');
figure; plot(l_video/fs, acf_video); grid on; title('АКФ видеоимпульса');
figure; plot(l_vp/fs, acf_video_pack); grid on; title('АКФ пачки видеоимпульсов');
figure; plot(l_radio/fs, acf_radio); grid on; title('АКФ радиоимпульса');
figure; plot(l_rp/fs, acf_radio_pack); grid on; title('АКФ пачки радиоимпульсов');

figure; plot(l_m/fs_m, acf_m); grid on; title('АКФ М-последовательности');
figure; plot(l_b/fs, acf_barker); grid on; title('АКФ кода Баркера 13');
figure; plot(l_lfm/fs, acf_lfm); grid on; title('АКФ ЛЧМ-сигнала');

