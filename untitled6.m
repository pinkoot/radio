% Параметры
n_rep = 3; % Число повторений импульса
T_total = 120e-3;       % общая длительность 120 мс
fs = 100e6;                % частота дискретизации 100 кГц
t = 0:1/fs:T_total;     

% Параметры импульсов
tau_imp = 20e-3;             % длительность одного импульса 20 мс
Ts = 40e-3;              % период следования импульсов 40 мс
n_pulses = floor(T_total/Ts);
q = 2; % Скважность импульса
T_imp = q*tau_imp; % Период импульса
f_diskr = 2e8; % Частота дискретизации
step_d = 1/f_diskr; % Шаг дискретизации
f_car = 500; % Несущая частота
time_of_modeling = 0:step_d:n_rep*T_imp; % Время моделирования
Periods = 0:T_imp:n_rep*T_imp; % Массив начал периодов
ampl = 1/512; % Амплитуда несущей
t_pulses = (0:n_pulses-1)*Ts;


% Функции моделируемых сигналов
Car_of_imp = ampl*sin(2*pi*f_car.*time_of_modeling); % Функция несущей
Video_imp=ampl * double(abs(time_of_modeling - tau_imp/2) < tau_imp/2);
f_radio = 500; 
Radio_imp = Video_imp .* Car_of_imp; % Функция радиоимпульса
Pocket_video_imp = pulstran(time_of_modeling, Periods, @rectpuls, tau_imp); % Функция пачки видеоимпульсов
Pocket_radio_imp = Pocket_video_imp.*Car_of_imp; % Функция пачки радиоимпульсов

%Гауссовский шум
M = 0;
D = 1;
WGN = sqrt(D)*randn(1, length(time_of_modeling))+M;



Car_of_imp_GS = Car_of_imp+WGN;
[car_corr, lags_corr] = xcorr(Car_of_imp, Car_of_imp_GS);

Video_imp_GS = Video_imp+WGN;
[video_corr, lags_video] = xcorr(Video_imp, Video_imp_GS,'coef');

Radio_imp_GS = Radio_imp+WGN;
[radio_corr, lags_radio] = xcorr(Radio_imp, Radio_imp_GS,'coef');

Pocket_radio_imp_GS = Pocket_radio_imp+WGN;
[pocket_radio_corr, lags_radio_pock] = xcorr(Pocket_radio_imp, Pocket_radio_imp_GS);

Pocket_video_imp_GS = Pocket_video_imp+WGN;
[pocket_video_corr, lags_video_pock] = xcorr(Pocket_video_imp, Pocket_video_imp_GS);

figure(1);
subplot(2, 1, 1);
plot(time_of_modeling, Car_of_imp_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда несущей частоты с шумом');
title('График несущего сигнала с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Car_of_imp);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда несущей частоты без шума');
title('График несущего сигнала без шума');

figure(2);
subplot(2, 1, 1);
plot(time_of_modeling, Video_imp_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда видеоимпульса с шумом');
title('График видеоимпульса с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Video_imp);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда видеоимпульса без шума');
title('График видеоимпульса без шума');

figure(3);
subplot(2, 1, 1);
plot(time_of_modeling, Radio_imp_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда радиоимпульса с шумом');
title('График радиоимпульса с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Radio_imp);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда радиоимпульса с шумом');
title('График радиоимпульса без шума');

figure(4);
subplot(2, 1, 1);
plot(time_of_modeling, Pocket_radio_imp_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда пачки радиоимпульсов с шумом');
title('График пачки радиоимпульсов с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Pocket_radio_imp);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда пачки радиоимпульсов без шума');
title('График пачки радиоимпульсов без шума');

figure(5);
subplot(2, 1, 1);
plot(time_of_modeling, Pocket_video_imp_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда пачки видеоимпульсов с шумом');
title('График пачки видеоимпульсов с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Pocket_video_imp);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда пачки радиоимпульсов без шума');
title('График пачки видеоимпульсов без шума');

figure(6);
plot(lags_corr, car_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции несущего сигнала с шумом и без');

figure(7);
plot(lags_video, video_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции видеосигнала с шумом и без');

figure(8);
plot(lags_radio, radio_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции радиосигнала с шумом и без');

figure(9);
plot(lags_radio_pock, pocket_radio_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции пачки радиоимпульсов с шумом и без');

figure(10);
plot(lags_video_pock, pocket_video_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции пачки видеоимпульсов с шумом и без');

% М-последовательность
f_diskr = 2e8; % Частота дискретизации
step_d = 1/f_diskr; % Шаг дискретизации
f_car = 50e3; % Несущая частота
time_of_modeling = 0:step_d:tau_imp; % Время моделирования
Car_of_imp = sin(2*pi*f_car.*time_of_modeling)/128; % Функция несущей

n_Mcode = 511; % Число символов в М-последовательности
t_M = tau_imp/n_Mcode; % Длительность одного символа в коде баркера

result = function_test();
counter = 1;
M_code=ones(1, length(time_of_modeling));
M_model=ones(1, length(time_of_modeling));
for i=1:length(time_of_modeling)
    if time_of_modeling(1, i) <= tau_imp
        if time_of_modeling(1, i) >= t_M*counter 
            if counter < 511
                counter = counter+1;
            end
        end
        M_code(1, i) = Car_of_imp(1, i)*result(1, counter); % кодом М-последовательности
        M_model(1, i) = result(1, counter);
    else
        M_code(1, i) = Car_of_imp(1, i); % несущей (после оконачания модуляции М-последовательностью)
    end
end

M = 0;
D = 1;
WGN = sqrt(D)*randn(1, length(time_of_modeling))+M;

M_code_GS = M_code+WGN;
[m_code_corr, lags_m_corr] = xcorr(M_code, M_code_GS);

figure(11);
subplot(2, 1, 1);
plot(time_of_modeling, M_code_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда М-последовательности с шумом');
title('График М-последовательности с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, M_model);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда М-последовательности без шума');
title('График М-последовательности без шума');

figure(12);
plot(lags_m_corr, m_code_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции М-последовательности с шумом и без');

%Код Баркера 13
Bbase_code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; % Код Баркера 13
t_bark = tau_imp/13; % Длительность одного символа в коде баркера
f_car = 6.5e3; % Несущая частота
f_diskr = f_car*10; % Частота дискретизации
step_d = 1/f_diskr; % Шаг дискретизации
time_of_modeling = 0:step_d:0.12; % Время моделирования
Car_of_imp = sin(2*pi*f_car.*time_of_modeling); % Функция несущей

counter = 1;
Barker_code=ones(1, length(time_of_modeling));
Barker_model=ones(1, length(time_of_modeling));
for i=1:length(time_of_modeling)
    if time_of_modeling(1, i) <= tau_imp
        if time_of_modeling(1, i) >= t_bark*counter % Условие перехода к следующему значению кода Баркера
            if counter < 13
                counter = counter+1;
            end
        end
        phase_shift = (1 - Bbase_code(counter)) * pi/2; % 0 или π
        Barker_code(i) = ampl * sin(2*pi*f_car*time_of_modeling(i) + phase_shift);
        Barker_model(i) = Bbase_code(counter);
    else
        Barker_code(i) = Car_of_imp(i);
        Barker_model(i) = 0;
    end
end

M = 0;
D = 1;
WGN = sqrt(D)*randn(1, length(time_of_modeling))+M;
Barker_code_GS = Barker_code+WGN;
[barker_code_corr, barker_code_lags] = xcorr(Barker_code, Barker_code_GS,'coef');


figure(13);
subplot(2, 1, 1);
plot(time_of_modeling, Barker_code_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда кода Баркера 13 с шумом');
title('График кода Баркера 13 с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Barker_model);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда кода Баркера 13 без шума');
title('График кода Баркера 13 без шума');

figure(14);
plot(barker_code_lags, barker_code_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции кода Баркера 13 с шумом и без');


% Одиночный импульс ЛЧМ-сигнала
f_diskr = 2e8; % Частота дискретизации
step_d = 1/f_diskr; % Шаг дискретизации
time_of_modeling = 0:step_d:tau_imp; % Время моделирования

f_growth = (1000/((length(time_of_modeling)))); % Установка коэффициента скорости роста частоты
freq=ones(1, length(time_of_modeling));
for i = 1:length(time_of_modeling)
        freq(i) = (i*f_growth+1000);
end
Lchm_sign = sin(2*pi*freq.*time_of_modeling); % Функция ЛЧМ-сигнала

%Гауссовский шум
M = 0;
D = 1;
WGN = sqrt(D)*randn(1, length(time_of_modeling))+M;
Lchm_sign_GS = Lchm_sign+WGN;
[Lchm_sign_corr, lags_lchm_sign] = xcorr(Lchm_sign, Lchm_sign_GS);


figure(15);
subplot(2, 1, 1);
plot(time_of_modeling, Lchm_sign_GS);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда кода ЛЧМ-сигнала с шумом');
title('График кода ЛЧМ-сигнала с шумом');

subplot(2, 1, 2);
plot(time_of_modeling, Lchm_sign);
grid;
xlabel('Время моделирования, с'); ylabel('Амплитуда кода ЛЧМ-сигнала без шума');
title('График кода ЛЧМ-сигнала без шума');

figure(16);
plot(lags_lchm_sign, Lchm_sign_corr);
grid;
xlabel('Время задержки, с'); ylabel('Взаимокорреляционная функция');
title('График взаимокорреляционной функции ЛЧМ сигнала с шумом и без');

figure(17);
plot(time_of_modeling, WGN);
grid;
xlabel('Время, с'); ylabel('Амплитуда шума');
title('График БГШ');

figure(18);
H1 = histogram(WGN, 100);
[corr, lags] = xcorr(WGN, 'coef');
title('Гистограмма БГШ');

figure(19);
plot(lags, corr);
xlabel('Время задержки, с'); ylabel('Автокорреляционная функция БГШ');
title('Автокорреляционная функция БГШ');

