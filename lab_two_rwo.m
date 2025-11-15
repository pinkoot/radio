%% --- ПАРАМЕТРЫ ---
fs = 1e5;       % Частота дискретизации (достаточно высокая)
dt = 1/fs;

T_video = 0.2;  % Длительность видеоимпульса
T_imp = 0.4;    % Период пачки
N_rep = 3;      % 3 импульса в пачке

t = -1:dt:1;    % Окно для визуализации (двустороннее)

%% --- ВИДЕОИМПУЛЬС (ПРЯМОУГОЛЬНЫЙ) ---
video = double( abs(t) <= T_video/2 );

% АКФ видеоимпульса → правильный треугольник
[acf_video, lag] = xcorr(video, 'coeff');
tau = lag * dt;

%% --- ПАЧКА ВИДЕОИМПУЛЬСОВ ---
starts = (-floor(N_rep/2)*T_imp) : T_imp : (floor(N_rep/2)*T_imp);

video_pack = zeros(size(t));
for k = 1:length(starts)
    video_pack = video_pack + double( abs(t - starts(k)) <= T_video/2 );
end

[acf_video_pack, lag_pack] = xcorr(video_pack, 'coeff');
tau_pack = lag_pack * dt;

%% === ГРАФИК: ВИДЕОИМПУЛЬС + ПАЧКА ===
figure;
subplot(2,1,1)
plot(tau, acf_video, 'LineWidth', 1.8);
grid on;
title('АКФ одиночного видеоимпульса');
xlabel('Время задержки, с'); ylabel('АКФ');

subplot(2,1,2)
plot(tau_pack, acf_video_pack, 'LineWidth', 1.8);
grid on;
title('АКФ пачки видеоимпульсов');
xlabel('Время задержки, с'); ylabel('АКФ');



%% --- РАДИОИМПУЛЬС ---
f_car = 50;    % несущая
radio = video .* cos(2*pi*f_car*t);

[acf_radio, lag_r] = xcorr(radio, 'coeff');
tau_r = lag_r * dt;

%% --- ПАЧКА РАДИОИМПУЛЬСОВ ---
radio_pack = zeros(size(t));
for k = 1:length(starts)
    radio_pack = radio_pack + ...
        double(abs(t - starts(k)) <= T_video/2) .* cos(2*pi*f_car*(t - starts(k)));
end

[acf_radio_pack, lag_rp] = xcorr(radio_pack, 'coeff');
tau_rp = lag_rp * dt;

%% === ГРАФИК: РАДИОИМПУЛЬС + ПАЧКА ===
figure;
subplot(2,1,1)
plot(tau_r, acf_radio, 'LineWidth', 1.7);
grid on;
title('АКФ одиночного радиоимпульса');
xlabel('Время задержки, с'); ylabel('АКФ');

subplot(2,1,2)
plot(tau_rp, acf_radio_pack, 'LineWidth', 1.7);
grid on;
title('АКФ пачки радиоимпульсов');
xlabel('Время задержки, с'); ylabel('АКФ');
