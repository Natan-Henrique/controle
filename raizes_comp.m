s = tf('s');
G = 1/(s+2);

P = 1/(s^2 + 1);

K_values = linspace(0.05, 0.2, 10);
a_values = linspace(1, 5, 10);
p_values = linspace(0.05, 0.5, 10);

best_cost = inf;
best_params = [];

beta = 10; % peso para amortecimento (ajuste conforme necessário)

for K = K_values
    for a = a_values
        for p = p_values
            C = K*(s + a) / (s*(s^2 + 1)*(s + p));
            T = feedback(C*G, 1);
            
            poles = pole(T);
            if all(real(poles) < 0)
                info = stepinfo(T);
                ts = info.SettlingTime;
                
                % Selecionar os polos dominantes (menor parte real)
                [~, idx] = sort(real(poles), 'descend'); 
                dominant = poles(idx(end-1:end)); % dois polos dominantes complexos
                
                sigma = real(dominant(1));
                omega_d = imag(dominant(1));
                
                % Calcular zeta
                zeta = -sigma / sqrt(sigma^2 + omega_d^2);
                
                % Evitar NaN ou valores inválidos
                if isnan(ts) || isnan(zeta) || ts <= 0 || zeta <= 0
                    continue;
                end
                
                % Custo: tempo de acomodação + peso para baixo zeta
                cost = ts + beta*(1/zeta);
                
                if cost < best_cost
                    best_cost = cost;
                    best_params = [K, a, p];
                    best_ts = ts;
                    best_zeta = zeta;
                end
            end
        end
    end
end

fprintf('Melhores parâmetros encontrados:\n');
fprintf('K = %.4f, a = %.4f, p = %.4f\n', best_params(1), best_params(2), best_params(3));
fprintf('Tempo de acomodação = %.2f s\n', best_ts);
fprintf('Coeficiente de amortecimento zeta = %.4f\n', best_zeta);

% Simulação com os melhores parâmetros
K = best_params(1);
a = best_params(2);
p = best_params(3);
C = K*(s + a) / (s*(s^2 + 1)*(s + p));
T = feedback(C*G, 1);

S = feedback(1, C*G);
T_p = S * P;

t = 0:0.1:200;
[y_r, t_out] = step(T, t);
u_p = sin(t);
y_p = lsim(T_p, u_p, t);
y_total = y_r + y_p;

figure;
plot(t_out, y_total, 'r', 'LineWidth', 1.5); hold on;
plot(t_out, ones(size(t_out)), 'k--', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('Saída Y(t)');
title('Resposta ao Degrau com Rejeição à Perturbação Senoidal');
legend('Saída Total', 'Referência');
grid on;
