% Configurar parámetros iniciales
filename = 'profilesvelocityRe10.txt';
cilindro_data = [];
eje_data = [];

% Abrir el archivo para lectura
fid = fopen(filename, 'r');
if fid == -1
    error('No se pudo abrir el archivo %s', filename);
end

% Leer el archivo línea por línea
current_section = '';
while ~feof(fid)
    line = strtrim(fgetl(fid));
    
    % Detectar secciones
    if startsWith(line, '((xy/key/label "cilindro")')
        current_section = 'cilindro';
        continue;
    elseif startsWith(line, '((xy/key/label "eje")')
        current_section = 'eje';
        continue;
    end
    
    % Saltar líneas no relevantes
    if isempty(line) || startsWith(line, '(') || ~contains(line, sprintf('\t'))
        continue;
    end
    
    % Procesar líneas de datos
    data = sscanf(line, '%f %f');
    if numel(data) == 2
        switch current_section
            case 'cilindro'
                cilindro_data = [cilindro_data; data'];
            case 'eje'
                eje_data = [eje_data; data'];
        end
    end
end
fclose(fid);

% Crear el gráfico
figure;
hold on;

% Graficar cilindro en rojo
if ~isempty(cilindro_data)
    plot(cilindro_data(:,1)+12+1.398-1, cilindro_data(:,2), 'r.', 'LineWidth', 1.5, 'DisplayName', 'Cilindro');
end

% Graficar eje en azul
if ~isempty(eje_data)
    plot(eje_data(:,1)+12+1.398-1, eje_data(:,2), 'b.', 'LineWidth', 1.5, 'DisplayName', 'Eje');
end

% Configurar el gráfico
title('Perfil de Velocidad X');
xlabel('Posición');
ylabel('Velocidad X');
legend('show');
grid on;
%hold off;
