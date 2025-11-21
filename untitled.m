clc; clear; close all;

%% 1. Одиночные сигналы (пункт а)

% Параметры
T_total = 70e-3;      % общая длительность сигнала 70 мс
T_pulse = 20e-3;      % длительность импульса 20 мс
fs = 1e5;             % частота дискретизации 1 кГц (достаточно для низкочастотных сигналов)
t = 0:1/fs:T_total-1/fs;   % временная ось

% 1. Гармоническое колебание
f_harm = 500;   % частота гармоники 500 Гц
u_harm = zeros(size(t));
pulse_start_harm = (T_total - T_pulse)/2;  
pulse_end_harm = pulse_start_harm + T_pulse;    
harm_indices = t >= pulse_start_harm & t < pulse_end_harm;
u_harm(harm_indices) = cos(2*pi*f_harm*(t(harm_indices) - pulse_start_harm));

% 2. Одиночный прямоугольный видеоимпульс
u_video = zeros(size(t));
pulse_start = (T_total - T_pulse)/2;  
pulse_end = pulse_start + T_pulse;    
u_video(t >= pulse_start & t < pulse_end) = 1;

% 3. Одиночный прямоугольный радиоимпульс
f_radio = 1e3;     % 1 кГц     % частота радиоимпульса 1 кГц
u_radio = zeros(size(t));
idx = t >= pulse_start & t < pulse_end;
u_radio(idx) = cos(2*pi*f_radio*(t(idx) - pulse_start));

% Построение графиков
figure('Name','Одиночные сигналы');
subplot(3,1,1)
plot(t*1e3, u_harm, 'b','LineWidth',1.5)
title('Гармоническое колебание')
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 T_total*1e3]); ylim([-1.5 1.5])

subplot(3,1,2)
plot(t*1e3, u_video, 'r','LineWidth',2)
title('Одиночный прямоугольный видеоимпульс')
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 T_total*1e3]); ylim([-0.3 1.3])

subplot(3,1,3)
plot(t*1e3, u_radio, 'k','LineWidth',1.5)
title('Одиночный прямоугольный радиоимпульс')
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 T_total*1e3]); ylim([-1.5 1.5])

%% 2. Пачка прямоугольных видео- и радиоимпульсов

T_total = 120e-3;       % общая длительность 120 мс
fs = 1e5;                % частота дискретизации 100 кГц
t = 0:1/fs:T_total;     

% Параметры импульсов
tau = 20e-3;             % длительность одного импульса 20 мс
Ts = 40e-3;              % период следования импульсов 40 мс
n_pulses = floor(T_total/Ts);
t_pulses = (0:n_pulses-1)*Ts;  % моменты появления импульсов

% 1. Несущее колебание
f_carrier = 50;         % частота несущей 50 Гц
u_carrier = cos(2*pi*f_carrier*t);

% 2. Пачка видеоимпульсов с pulstran
video_pulse = @(t) double(abs(t)<tau/2);
u_video = pulstran(t, t_pulses, video_pulse);

% 3. Пачка радиоимпульсов с pulstran
f_radio = 500; 
radio_pulse = @(t) cos(2*pi*f_radio*t) .* (abs(t)<tau/2);
u_radio = pulstran(t, t_pulses, radio_pulse);

% Построение графиков
figure('Name','Пачка видео- и радиоимпульсов');
subplot(3,1,1)
plot(t*1e3, u_carrier, 'b','LineWidth',1.5)
title('Несущее колебание'); xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-1.2 1.2])

subplot(3,1,2)
plot(t*1e3, u_video, 'r','LineWidth',1.5)
title('Пачка видеоимпульсов'); xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-0.2 1.2])

subplot(3,1,3)
plot(t*1e3, u_radio, 'g','LineWidth',1.5)
title('Пачка радиоимпульсов'); xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-1.2 1.2])


%% 3. Одиночный ФКМ импульс (М-последовательность, код Баркера)

% Параметры
tau_imp = 0.002; % Длительность импульса 2 мс
f_car = 500; 
f_diskr = 1e5; % Частота дискретизации 100 кГц
step_d = 1/f_diskr; 
time_of_modeling = 0:step_d:tau_imp;
ampl = 1;

% Несущее колебание
Car_of_imp = ampl * sin(2*pi*f_car.*time_of_modeling);

% Генерация М-последовательности
function PSP = function_test()
    v = [-1 -1 -1 -1 -1 -1 -1 -1 -1];
    for i = 1:511
        PSP(i) = v(7);
        f = v(5) * v(9);
        v(2:9) = v(1:8);
        v(1) = f;
    end
    PSP = -PSP;
end

M_sequence = function_test();
n_Mcode = length(M_sequence);
t_M = tau_imp / n_Mcode;

% Формирование ФКМ сигнала с М-последовательностью
M_code = zeros(1, length(time_of_modeling));
M_model = ones(1, length(time_of_modeling));
counter = 1;
for i = 1:length(time_of_modeling)
    if time_of_modeling(i) <= tau_imp
        if time_of_modeling(i) >= t_M * counter
            if counter < n_Mcode
                counter = counter + 1;
            end
        end
        phase_shift = (1 - M_sequence(counter)) * pi/2; % 0 или π
        M_code(i) = ampl * sin(2*pi*f_car*time_of_modeling(i) + phase_shift);
        M_model(i) = M_sequence(counter);
    else
        M_code(i) = Car_of_imp(i);
        M_model(i) = 1;
    end
end       
% Код Баркера длиной 13
Bbase_code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
t_bark = tau_imp / length(Bbase_code);
Barker_code = zeros(1, length(time_of_modeling));
Barker_model = ones(1, length(time_of_modeling));
counter = 1;
for i = 1:length(time_of_modeling)
    if time_of_modeling(i) <= tau_imp
        if time_of_modeling(i) >= t_bark * counter
            if counter < length(Bbase_code)
                counter = counter + 1;
            end
        end
        phase_shift = (1 - Bbase_code(counter)) * pi/2; % 0 или π
        Barker_code(i) = ampl * sin(2*pi*f_car*time_of_modeling(i) + phase_shift);
        Barker_model(i) = Bbase_code(counter);
    else
        Barker_code(i) = Car_of_imp(i);
        Barker_model(i) = 1;
    end
end

% Построение графиков ФКМ сигналов
figure('Name','ФКМ сигналы');
subplot(3,2,1); plot(time_of_modeling, Car_of_imp,'LineWidth',1.5); title('Несущее колебание'); grid on;
subplot(3,2,3); plot(time_of_modeling, M_model,'LineWidth',1.5); title('М-последовательность'); grid on;
subplot(3,2,5); plot(time_of_modeling, M_code,'LineWidth',1.5); title('ФКМ с М-последовательностью'); grid on;
subplot(3,2,2); plot(time_of_modeling, Car_of_imp,'LineWidth',1.5); title('Несущее колебание'); grid on;
subplot(3,2,4); plot(time_of_modeling, Barker_model,'LineWidth',2); title('Код Баркера'); grid on;
subplot(3,2,6); plot(time_of_modeling, Barker_code,'LineWidth',1.5); title('ФКМ с кодом Баркера'); grid on;


%% 4. ЛЧМ импульс

% Параметры ЛЧМ сигнала
f0 = 1000;              % Начальная частота, Гц
f1 = 5000;              % Конечная частота, Гц
T = tau_imp;             % Длительность равна длительности ФКМ импульса
k = (f1 - f0)/T;        % Градиент частоты
Lchm_sign = cos(2*pi*(f0*time_of_modeling + 0.5*k*time_of_modeling.^2));

% Построение графиков ЛЧМ сигнала
figure(4); % ЛЧМ сигнал
subplot(2,1,1) 
plot(time_of_modeling, Car_of_imp, 'LineWidth', 1.5); % График несущей
grid on;
xlabel('Время, с'); ylabel('Амплитуда');
title('График несущей сигнала');
axis([0 max(time_of_modeling) -1.2 1.2])

subplot(2,1,2) 
plot(time_of_modeling, Lchm_sign, 'LineWidth', 1.5); % График ЛЧМ
grid on;
xlabel('Время, с'); ylabel('Амплитуда');
title('График одиночного импульса ЛЧМ сигнала');
axis([0 max(time_of_modeling) -1.2 1.2])
% Код Баркера
% Параметры сигнала
tau_imp = 0.02; % Длительность импульса 20 мс
f_car = 500; % Несущая частота 500 Гц
f_diskr = 10e4; % Частота дискретизации 100 кГц
step_d = 1/f_diskr; % Шаг дискретизации
ampl = 1; % Амплитуда

time_of_modeling = 0:step_d:0.12; % 120 мс
t = time_of_modeling;

% Несущее колебание
u_carrier = ampl * sin(2*pi*f_car.*t);

% Код Баркера длиной 13 + нули после окончания
Bbase_code = [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
% Добавляем нули после кода Баркера
extended_code = [Bbase_code, 0, 0, 0, 0, 0]; % 5 нулей после кода

t_symbol = tau_imp / length(Bbase_code); % Длительность одного символа Баркера
total_symbols = length(extended_code); % Общее количество символов

% Формирование модулирующей функции (с нулями после кода Баркера)
Barker_model = zeros(1, length(t));
for i = 1:length(t)
    symbol_index = min(floor(t(i) / t_symbol) + 1, total_symbols);
    if symbol_index <= total_symbols
        Barker_model(i) = extended_code(symbol_index);
    else
        Barker_model(i) = 0; % После всех символов - ноль
    end
end

% Формирование ФКМ сигнала с возвратом к несущей при нуле
Barker_code = zeros(1, length(t));
for i = 1:length(t)
    if Barker_model(i) == -1
        Barker_code(i) = -u_carrier(i); % Сдвиг фазы на 180°
    elseif Barker_model(i) == 0
        Barker_code(i) = 0; % НУЛЬ после окончания кода Баркера
    else
        Barker_code(i) = u_carrier(i); % Фаза без изменений
    end
end

% Построение графиков
figure('Name','ФКМ сигнал с кодом Баркера (исправленный)');
subplot(3,1,1)
plot(t*1e3, u_carrier, 'b','LineWidth',1.5)
title('Несущее колебание (500 Гц)'); 
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-1.2 1.2])

subplot(3,1,2)
stairs(t*1e3, Barker_model, 'r','LineWidth',2)
title('Модулирующая функция (код Баркера N=13 + нули)'); 
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-1.5 1.5])
yticks([-1, 0, 1])
yticklabels({'-1', '0', '1'})

subplot(3,1,3)
plot(t*1e3, Barker_code, 'g','LineWidth',1.5)
title('ФКМ сигнал (фазовая манипуляция с возвратом в ноль)'); 
xlabel('Время, мс'); ylabel('Амплитуда'); grid on
xlim([0 120]); ylim([-1.2 1.2])

% Добавляем вертикальные линии для наглядности переключений фаз
hold on
for i = 1:length(extended_code)
    x_line = (i-1) * t_symbol * 1e3;
    if extended_code(i) == -1
        plot([x_line x_line], [-1.2 1.2], 'r:', 'LineWidth', 0.5, 'Color', [1 0 0 0.3])
    elseif extended_code(i) == 0
        plot([x_line x_line], [-1.2 1.2], 'k:', 'LineWidth', 1, 'Color', [0 0 0 0.5])
    else
        plot([x_line x_line], [-1.2 1.2], 'b:', 'LineWidth', 0.5, 'Color', [0 0 1 0.3])
    end
end
hold off

% Вывод информации о сигнале
fprintf('\n=== Параметры исправленного ФКМ сигнала ===\n');
fprintf('Длительность импульса: %.0f мс\n', tau_imp*1000);
fprintf('Несущая частота: %.0f Гц\n', f_car);
fprintf('Длина кода Баркера: %d символов\n', length(Bbase_code));
fprintf('Добавлено нулей после кода: %d\n', length(extended_code) - length(Bbase_code));
fprintf('Длительность символа: %.3f мс\n', t_symbol*1000);
fprintf('Полный код: [');
for i = 1:length(extended_code)
    if extended_code(i) == 1
        fprintf('+1 ');
    elseif extended_code(i) == -1
        fprintf('-1 ');
    else
        fprintf(' 0 ');
    end
end
fprintf(']\n');
% Спектр гармонического сигнала
N_harm = 2^nextpow2(length(u_harm));
spec_harmo = fft(u_harm, N_harm);
amp_harmo = abs(spec_harmo) ./ length(u_harm); % Нормировка на длину сигнала!

% Построение графика спектра гармонического колебания
figure(6);
fs = 1e5;
F_range = (0:N_harm/2-1) * (fs/N_harm);

plot(F_range, amp_harmo(1:N_harm/2), 'LineWidth', 2, 'Color', 'blue');
title('Спектр гармонического сигнала', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Амплитуда', 'FontSize', 12);
grid on;

% Умное центрирование вокруг частоты сигнала
[peak_amp, peak_idx] = max(amp_harmo(1:N_harm/2));
peak_freq = F_range(peak_idx);

% Центрируем вокруг пика с запасом ±25%
f_span = peak_freq * 0.25;
xlim([peak_freq - f_span, peak_freq + f_span]);
ylim([0, peak_amp * 1.1]);

% Вывод информации в командное окно
fprintf('\n=== СПЕКТРАЛЬНЫЙ АНАЛИЗ ГАРМОНИЧЕСКОГО СИГНАЛА ===\n');
fprintf('Расчетная частота: %.2f МГц\n', f_harm/1e6);
fprintf('Фактическая частота пика: %.2f МГц\n', peak_freq/1e6);
fprintf('Амплитуда пика: %.6f\n', peak_amp);
fprintf('Отклонение: %.4f%%\n', abs(peak_freq - f_harm)/f_harm * 100);
fprintf('Частотное разрешение: %.2f кГц\n', fs/N_harm/1000);

%% Аналогично для других сигналов

%% Спектр прямоугольного видеоимпульса
T_total = 0.07;      % общая длительность сигнала 70 мс
T_pulse = 0.02;      % длительность импульса 20 мс
fs = 1e5;            % частота дискретизации 100 кГц
t = 0:1/fs:T_total-1/fs;   % временная ось

% Формирование видеоимпульса
u_video = zeros(size(t));
pulse_start = (T_total - T_pulse)/2;  
pulse_end = pulse_start + T_pulse;    
u_video(t >= pulse_start & t < pulse_end) = 1;

% Расчет спектра
N_video = 2^nextpow2(length(u_video)); 
spec_video = fft(u_video, N_video);
amp_video = abs(spec_video) ./ length(u_video); % нормировка на длину сигнала
amp_video_shifted = fftshift(amp_video);

% Ось частот
F_range_video = (-N_video/2:N_video/2-1) * (fs/N_video);

% Построение спектра
figure;
plot(F_range_video/1e3, amp_video_shifted, 'b', 'LineWidth', 2);
title('спектр видеоимпульса', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('f, Гц', 'FontSize', 12);
ylabel('A', 'FontSize', 12);
grid on;

% Диапазон, как на изображении
xlim([-5/T_pulse 5/T_pulse]); % ±5 главных лепестков
ylim([0 max(amp_video_shifted)*1.1]);

%% Спектр прямоугольного радиоимпульса (симметрия относительно 500 Гц)
% Параметры
tau_imp = 0.02; % Время одного импульса
n_rep = 3; % Число повторений импульса
q = 2; % Скважность импульса
T_imp = q*tau_imp; % Период импульса
f_diskr = 10e4; % Частота дискретизации
step_d = 1/f_diskr; % Шаг дискретизации
f_car = 500; % Несущая частота
time_of_modeling = 0:step_d:n_rep*T_imp; % Время моделирования
Periods = 0:T_imp:n_rep*T_imp; % Массив начал периодов
ampl = 1; % Амплитуда несущей

% Функции моделируемых сигналов
Car_of_imp = sin(2*pi*f_car.*time_of_modeling); % Функция несущей
Video_imp = ampl*rectpuls(time_of_modeling, q*tau_imp); % Функция видеоимпульса
Radio_imp = Video_imp.*Car_of_imp; % Функция радиоимпульса
n_count_for_fft = 2^nextpow2(length(time_of_modeling));
L = max(time_of_modeling)*f_diskr;
fft_rad_imp = fft(Radio_imp, n_count_for_fft); % БПФ для радиоимпульса
fft_rad_imp = fft_rad_imp(1:n_count_for_fft/2); % Отсечение половины частот
Real_fft_rad_imp = real(fft_rad_imp); % Получение действительная части
Imagine_fft_rad_imp = imag(fft_rad_imp); % Получение мнимой части
Ampl_fft_rad_imp = abs(fft_rad_imp)/L; % Получение модулей комплексных чисел
Frequencies = (0:n_count_for_fft-1)*f_diskr/n_count_for_fft; % Получение частоты
Frequencies = Frequencies(1:n_count_for_fft/2); % Отсечение половины частот
figure(8);
subplot(2, 1, 1);
plot(Frequencies, Ampl_fft_rad_imp);
xlabel('Частота, Гц'); ylabel('Амплитуда');
title('График амплитудного спектра радиоимпульса');
axis([0 2*f_car -0.2*max(Ampl_fft_rad_imp) 1.2*max(Ampl_fft_rad_imp)]);


% Пачка прямоугольных видео и радиоимпульсов
T_total = 120e-3;       % общая длительность 120 мс
fs = 1e5;                % частота дискретизации 100 кГц
t = 0:1/fs:T_total;     

% Параметры импульсов
tau = 20e-3;             % длительность одного импульса 20 мс
Ts = 40e-3;              % период следования импульсов 40 мс
n_pulses = floor(T_total/Ts);
t_pulses = (0:n_pulses-1)*Ts;  % моменты появления импульсов

% 1. Несущее колебание
f_carrier = 500;         % частота несущей 500 Гц
u_carrier = cos(2*pi*f_carrier*t);

% 2. Пачка видеоимпульсов 
video_pulse = @(t) double(abs(t)<tau/2);
u_video = pulstran(t, t_pulses, video_pulse);
N_video_pack = 2^nextpow2(length(u_video));
spec_video_pack = fft(u_video, N_video_pack);
amp_video_pack = abs(spec_video_pack) ./ N_video_pack;
amp_video_pack_shifted = fftshift(amp_video_pack);

F_range_video_pack = (-N_video_pack/2:N_video_pack/2-1) * (fs/N_video_pack);

figure(9);
plot(F_range_video_pack, amp_video_pack_shifted, 'LineWidth', 2, 'Color', 'magenta');
title('Спектр пачки прямоугольных видеоимпульсов', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Амплитуда', 'FontSize', 12);
grid on;

% Покажем низкие частоты (основной лепесток)
xlim([-500, 500]);
fprintf('\n=== СПЕКТР ПАЧКИ ВИДЕОИМПУЛЬСОВ ===\n');
fprintf('Период следования импульсов: %.1f мс\n', Ts*1e3);
fprintf('Длительность импульса: %.1f мс\n', tau*1e3);
fprintf('Теоретическая частота повторения импульсов: %.1f Гц\n', 1/Ts);
fprintf('Главный лепесток около 0 Гц (видеоимпульсы)\n');


% Пачка прямоугольных радиоимпульсов
f_radio = 500; 
radio_pulse = @(t) cos(2*pi*f_radio*t) .* (abs(t)<tau/2);
u_radio = pulstran(t, t_pulses, radio_pulse);
N_radio_pack = 2^nextpow2(length(u_radio));
spec_radio_pack = fft(u_radio, N_radio_pack);
amp_radio_pack = abs(spec_radio_pack) ./ N_radio_pack;

F_range_radio_pack = (0:N_radio_pack/2-1) * (fs/N_radio_pack);

figure(10);
plot(F_range_radio_pack, amp_radio_pack(1:N_radio_pack/2), 'LineWidth', 2, 'Color', 'green');
title('Спектр пачки прямоугольных радиоимпульсов', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Амплитуда', 'FontSize', 12);
grid on;

% Поиск пика несущей частоты
[peak_amp_pack, peak_idx_pack] = max(amp_radio_pack(1:N_radio_pack/2));
peak_freq_pack = F_range_radio_pack(peak_idx_pack);

% Окно вокруг несущей
f_span_pack = 200; % +/- 200 Гц достаточно, т.к. несущая = 500 Гц
xlim([f_radio - f_span_pack, f_radio + f_span_pack]);
ylim([0, peak_amp_pack*1.1]);

fprintf('\n=== СПЕКТР ПАЧКИ РАДИОИМПУЛЬСОВ ===\n');
fprintf('Несущая частота: %.1f Гц\n', f_radio);
fprintf('Период повторения импульсов: %.1f мс\n', Ts*1e3);
fprintf('Длительность импульса: %.1f мс\n', tau*1e3);
fprintf('Фактическая несущая в спектре: %.2f Гц\n', peak_freq_pack);
fprintf('Частотное разрешение: %.3f Гц\n', fs/N_radio_pack);




%% Создание сводной таблицы спектральных характеристик
fprintf('\n=== СВОДНАЯ ТАБЛИЦА СПЕКТРАЛЬНЫХ ХАРАКТЕРИСТИК ===\n');
fprintf('%-25s %-12s %-12s %-15s\n', 'Сигнал', 'Частота, МГц', 'Амплитуда', 'Отклонение');
fprintf('%-25s %-12.2f %-12.6f %-15s\n', 'Гармонический', peak_freq/1e6, peak_amp, ...
    sprintf('%.4f%%', abs(peak_freq - f_harm)/f_harm * 100));
fprintf('%-25s %-12s %-12.6f %-15s\n', 'Видеоимпульс', 'DC', max(amp_video), 'N/A');
%% ===========================================================
%% 6. Фазовые спектры всех сигналов
%% ===========================================================

%% Фазовый спектр гармонического сигнала 
fs = 1e5;
N_harm = length(spec_harmo);
F_phase_harmo = (-N_harm/2:N_harm/2-1) * (fs/N_harm);

% Расчет фазы
phase_harmo = angle(spec_harmo);
phase_harmo_shifted = fftshift(phase_harmo);

% Построение графика
figure(11);
plot(F_phase_harmo, phase_harmo_shifted, 'LineWidth', 1.5, 'Color', 'b');
title('Фазовый спектр гармонического сигнала', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Фаза, рад', 'FontSize', 12);
grid on;

% Умное центрирование вокруг частоты сигнала
[~, peak_idx] = max(abs(spec_harmo));
peak_freq_shifted = F_phase_harmo(peak_idx);
f_span = 2e2; % ±0.2 кГц вокруг пика

xlim([peak_freq_shifted - f_span, peak_freq_shifted + f_span]);
ylim([-pi pi]);


%% Фазовый спектр прямоугольного видеоимпульса 
fs = 1e5;

N_video = length(u_video);
spec_video = fft(u_video);
F_phase_video = (-N_video/2:N_video/2-1) * (fs/N_video);

phase_video = angle(spec_video);
phase_video_shifted = fftshift(phase_video);

figure(12);
plot(F_phase_video/1e3, phase_video_shifted, 'LineWidth', 2, 'Color', 'r');
title('Фазовый спектр прямоугольного видеоимпульса', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, кГц', 'FontSize', 12);
ylabel('Фаза, рад', 'FontSize', 12);
grid on;

xlim([-2, 2]); 
ylim([-pi, pi]);

Phase_rad_imp = atan2(Imagine_fft_rad_imp, Real_fft_rad_imp);
figure(13);
subplot(2, 1, 2);
plot(Frequencies, Phase_rad_imp);
xlabel('Частота, Гц'); ylabel('Фаза');
title('График фазового спектра радиоимпульса');
axis([0 2*f_car -pi-0.5 pi+0.5]);


% --- Фазовый спектр пачки прямоугольных видеоимпульсов ---
phase_video_pack = angle(spec_video_pack);
phase_video_pack_shifted = fftshift(phase_video_pack);

figure(14);
plot(F_range_video_pack, phase_video_pack_shifted, 'LineWidth', 1.5, 'Color', 'm');
title('Фазовый спектр пачки прямоугольных видеоимпульсов', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Фаза, рад', 'FontSize', 12);
grid on;
xlim([-500, 500]);
fprintf('\n=== ФАЗОВЫЙ СПЕКТР ПАЧКИ ВИДЕОИМПУЛЬСОВ ===\n');
fprintf('Фаза показывает периодическую структуру, связанную с повторением импульсов.\n');


% --- Фазовый спектр пачки радиоимпульсов ---
f_radio = 500;
phase_radio_pack = angle(spec_radio_pack);
figure(15);
plot(F_range_radio_pack, phase_radio_pack(1:N_radio_pack/2), 'LineWidth', 1.5, 'Color', 'g');
title('Фазовый спектр пачки прямоугольных радиоимпульсов', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Частота, Гц', 'FontSize', 12);
ylabel('Фаза, рад', 'FontSize', 12);
grid on;
xlim([f_radio - f_span_pack, f_radio + f_span_pack]);
fprintf('\n=== ФАЗОВЫЙ СПЕКТР ПАЧКИ РАДИОИМПУЛЬСОВ ===\n');
fprintf('Фазовый спектр содержит колебания вокруг несущей %.1f Гц.\n', f_radio);

% --- Спектр ЛЧМ сигнала 
f_growth = (2/(tau_imp*f_diskr)); % Установка коэффициента скорости роста частоты
freq=ones(1, length(time_of_modeling));
for i = 1:length(time_of_modeling)
    if i <= tau_imp*f_diskr
        freq(i) = i*f_growth;
    end
end
Lchm_sign = sin(2*pi*f_car*freq.*time_of_modeling); 
fft_LCHM = fft(Lchm_sign, n_count_for_fft); % БПФ для ЛЧМ сигнала
fft_LCHM = fft_LCHM(1:n_count_for_fft/2); % Отсечение половины частот
Real_fft_LCHM = real(fft_LCHM); % Получение действительная части
Imagine_fft_LCHM = imag(fft_LCHM); % Получение мнимой части
Ampl_fft_LCHM = abs(fft_LCHM)/L; % Получение модулей комплексных чисел
Phase_LCHM = atan2(Imagine_fft_LCHM, Real_fft_LCHM); % Получение фаз
figure(16);
subplot(2, 1, 1);
plot(Frequencies, Ampl_fft_LCHM);
xlabel('Частота, Гц'); ylabel('Амплитуда');
title('График амплитудного спектра ЛЧМ сигнала');
axis([0 2*f_car -0.2*max(Ampl_fft_LCHM) 1.2*max(Ampl_fft_LCHM)]);

subplot(2, 1, 2);
plot(Frequencies, Phase_LCHM);
xlabel('Частота, Гц'); ylabel('Фаза');
title('График фазового спектра ЛЧМ сигнала');
axis([0 2*f_car -pi-0.5 pi+0.5]);


% -- Вычисление спектра сигнала М-последовательности
n_Mcode = 511; % Число символов в М-последовательности
t_M = tau_imp/n_Mcode; % Длительность одного символа в коде баркера

result = function_test();
counter = 1;
M_code=ones(1, length(time_of_modeling));
M_model=ones(1, length(time_of_modeling));
for i=1:length(time_of_modeling) % Формирование модулированного в соответствии с
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

fft_m_posl = fft(M_code, n_count_for_fft); % БПФ для М-последовательности
fft_m_posl = fft_m_posl(1:n_count_for_fft/2); % Отсечение половины частот
Real_fft_m_posl = real(fft_m_posl); % Получение действительная части
Imagine_fft_m_posl = imag(fft_m_posl); % Получение мнимой части
Ampl_fft_m_posl = abs(fft_m_posl)/L; % Получение модулей комплексных чисел
Phase_m_posl = atan2(Imagine_fft_m_posl, Real_fft_m_posl); % Получение фаз
figure(17);
subplot(2, 1, 1);
plot(Frequencies, Ampl_fft_m_posl);
xlabel('Частота, Гц'); ylabel('Амплитуда');
title('График амплитудного спектра сигнала М-последовательности');
axis([0 2*f_car -0.2*max(Ampl_fft_m_posl) 1.2*max(Ampl_fft_m_posl)]);

subplot(2, 1, 2);
plot(Frequencies, Phase_m_posl);
xlabel('Частота, Гц'); ylabel('Фаза');
title('График фазового спектра сигнала М-последовательности');
axis([0 2*f_car -pi-0.5 pi+0.5]);


% Вычисление спектра сигнала кода Баркера 13
Bbase_code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; % Код Баркера 13
t_bark = tau_imp/13; % Длительность одного символа в коде баркера

counter = 1;
Barker_code=ones(1, length(time_of_modeling));
Barker_model=ones(1, length(time_of_modeling));
for i=1:length(time_of_modeling) % Формирование модулированного в соответствии с
    if time_of_modeling(1, i) <= tau_imp
        if time_of_modeling(1, i) >= t_bark*counter % Условие перехода к следующему значению кода Баркера
            if counter < 13
                counter = counter+1;
            end
        
        end
        Barker_code(1, i) = Bbase_code(1, counter)*Car_of_imp(1, i); % кодом Баркера 13 сигнала
        Barker_model(1, i) = Bbase_code(1, counter);
    else
        Barker_code(1, i) = Car_of_imp(1, i); % кодом Баркера 13 сигнала
        Barker_model(1, i) = 1;
    end
end

fft_bark_code = fft(Barker_code, n_count_for_fft); % БПФ для сигнала, модулированного кодом Баркера 13
fft_bark_code = fft_bark_code(1:n_count_for_fft/2); % Отсечение половины частот
Real_fft_barker_code = real(fft_bark_code); % Получение действительная части
Imagine_fft_barker_code = imag(fft_bark_code); % Получение мнимой части
Ampl_fft_barker_code = abs(fft_bark_code)/L; % Получение модулей комплексных чисел
Phase_barker_code = atan2(Imagine_fft_barker_code, Real_fft_barker_code); % Получение фаз
figure(18);
subplot(2, 1, 1);
plot(Frequencies, Ampl_fft_barker_code);
xlabel('Частота, Гц'); ylabel('Амплитуда');
title('График амплитудного спектра сигнала кода Баркера 13');
axis([0 2*f_car -0.2*max(Ampl_fft_barker_code) 1.2*max(Ampl_fft_barker_code)]);

subplot(2, 1, 2);
plot(Frequencies, Phase_barker_code);
xlabel('Частота, Гц'); ylabel('Фаза');
title('График фазового спектра кода сигнала кода Баркера 13');
axis([0 2*f_car -pi-0.5 pi+0.5]);


%% =======================  ФУНКЦИИ ДЛЯ СПЕКТРОВ ==========================
function plot_amp_spectrum(signal, fs, title_text)
    N = length(signal);
    Nfft = 2^nextpow2(N);
    S = fft(signal, Nfft);
    A = abs(S)/N;
    F = (0:Nfft/2-1) * (fs/Nfft);

    figure; stem(F, A(1:Nfft/2), 'LineWidth', 1.3);
    grid on; xlabel('Частота, Гц'); ylabel('Амплитуда');
    title(title_text); axis([0 max(F) -inf inf]);
end

function plot_phase_spectrum(signal, fs, title_text)
    N = length(signal);
    Nfft = 2^nextpow2(N);
    S = fft(signal, Nfft);
    PH = angle(S);
    F = (0:Nfft/2-1) * (fs/Nfft);

    figure; stem(F, PH(1:Nfft/2), 'LineWidth', 1.3);
    grid on; xlabel('Частота, Гц'); ylabel('Фаза, рад');
    title(title_text); axis([0 max(F) -pi pi]);
end
%% =======================================================================



%% =======================================================================
%                          ВХОДНЫЕ ПАРАМЕТРЫ
% ========================================================================

fs = 5e5;          % Частота дискретизации
t_end = 10e-3;     % Длительность анализа
t = 0 : 1/fs : t_end-1/fs;



%% =======================================================================
%                         1. ГАРМОНИЧЕСКИЙ СИГНАЛ
% ========================================================================

A = 1;
F0 = 10e3;  % 10 kHz гармонический

u_harm = A*cos(2*pi*F0*t);

figure; plot(t, u_harm);
title('Гармонический сигнал'); xlabel('t, c'); ylabel('u');
grid on;

plot_amp_spectrum(u_harm, fs, 'Амплитудный спектр гармонического сигнала');
plot_phase_spectrum(u_harm, fs, 'Фазовый спектр гармонического сигнала');



%% =======================================================================
%                           2. ВИДЕОИМПУЛЬС
% ========================================================================

tau = 300e-6;
u_video = double(t < tau);

figure; plot(t, u_video);
title('Видеоимпульс'); xlabel('t, c'); ylabel('u');
grid on;

plot_amp_spectrum(u_video, fs, 'Амплитудный спектр видеоимпульса');
plot_phase_spectrum(u_video, fs, 'Фазовый спектр видеоимпульса');



%% =======================================================================
%                      3. РАДИОИМПУЛЬС (несущая + окно)
% ========================================================================

fc = 40e3;
u_radio = u_video .* cos(2*pi*fc*t);

figure; plot(t, u_radio);
title('Радиоимпульс'); grid on;

plot_amp_spectrum(u_radio, fs, 'Амплитудный спектр радиоимпульса');
plot_phase_spectrum(u_radio, fs, 'Фазовый спектр радиоимпульса');



%% ====== Фазовый спектр М-последовательности  ======

% FFT
Nfft = 2^nextpow2(length(M_code));
SpecM = fft(M_code, Nfft);

% Центрирование спектра
SpecM_shift = fftshift(SpecM);

% Амплитуда и фаза
AmpM = abs(SpecM_shift);
PhaseM = angle(SpecM_shift);

% Развёртка фазы (устранение скачков ±π)
PhaseM_unwrap = unwrap(PhaseM);

% Частоты
F = (-Nfft/2:Nfft/2-1) * (f_diskr / Nfft);

figure;
subplot(2,1,1);
plot(F, AmpM, 'LineWidth', 1.5);
title('Амплитудный спектр М-последовательности');
xlabel('Частота, Гц'); ylabel('|S(f)|'); grid on;
xlim([f_car-800 f_car+800]);

subplot(2,1,2);
plot(F, PhaseM_unwrap, 'LineWidth', 1.5);
title('Фазовый спектр М-последовательности (развёрнутый)');
xlabel('Частота, Гц'); ylabel('Фаза, рад'); grid on;
xlim([f_car-800 f_car+800]);


%% =======================================================================
%                           5. КОД БАРКЕРА
% ========================================================================

Barker = [1 1 1 -1 -1 1 -1];
B_signal = repelem(Barker, round(fs/20e3));  % код на несущей 20 кГц

figure; plot(B_signal); title('Код Баркера'); grid on;

plot_phase_spectrum(B_signal, fs, 'Фазовый спектр кода Баркера');



%% =======================================================================
%                   6. ЛИНЕЙНО-ЧАСТОТНАЯ МОДУЛЯЦИЯ (ЛЧМ)
% ========================================================================

f1 = 10e3;
f2 = 50e3;
T_chirp = 1e-3;

t_ch = 0:1/fs:T_chirp-1/fs;

LFM = cos(2*pi*(f1*t_ch + (f2-f1)/(2*T_chirp)*t_ch.^2));

figure; plot(t_ch, LFM); title('ЛЧМ сигнал'); grid on;

plot_phase_spectrum(LFM, fs, 'Фазовый спектр ЛЧМ');
