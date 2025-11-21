figure(19);

subplot(3,1,1)
plot(time_of_modeling, Car_of_imp, 'LineWidth', 1.5); % Построение графика несущего сигнала
grid;
xlabel('Время, с'); ylabel('Амплитуда');
title('График несущего сигнала');
axis([0 max(time_of_modeling) -1.2 1.2])

subplot(3,1,2)
plot(time_of_modeling, M_model, 'LineWidth', 1.5); % Построение графика несущего сигнала
grid;
xlabel('Время, с'); ylabel('Амплитуда');
title('График кода М-последовательности');
axis([0 max(time_of_modeling) -1.2 1.2])

subplot(3,1,3) 
plot(time_of_modeling,M_code, 'LineWidth', 1.5); % Построение графика сигнала, модулированного М-последовательностью
grid;
xlabel('Время, с'); ylabel('Амплитуда');
title('График одиночного импульса сигнала M-последовательности');
axis([0 max(time_of_modeling) -1.2 1.2])