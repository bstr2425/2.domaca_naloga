
path_vozlisc = './vozlisca_temperature_dn2.txt';
path_cells = './celice_dn2.txt';

% x y in temp
opts = detectImportOptions(path_vozlisc, 'NumHeaderLines', 4); 
opts.VariableNames = {'x', 'y', 'temperatura'};
opts.Delimiter = ',';
vozlisca_data = readtable(path_vozlisc, opts);



% pt1, pt2, pt3, pt4

opts_cells = detectImportOptions(path_cells, 'NumHeaderLines', 2); 
opts_cells.VariableNames = {'pt1', 'pt2', 'pt3', 'pt4'};
opts_cells.Delimiter = ',';
celice_data = readtable(path_cells, opts_cells);




% razdelitev na x y in temperatura vektorje
x = vozlisca_data.x;
y = vozlisca_data.y;
temperatura = vozlisca_data.temperatura;

% metoda ScatteredInterpolant 
tic;
F_scattered = scatteredInterpolant(x, y, temperatura, 'linear', 'none');
T1 = F_scattered(0.403, 0.503);
time_scattered = toc;

% GriddedInterpolant metoda
tic;
[x_unique, y_unique, T_grid] = pripraviMreznePodatke(x, y, temperatura);
[X_nd, Y_nd] = ndgrid(x_unique, y_unique); 
F_gridded = griddedInterpolant(X_nd, Y_nd, T_grid, 'linear', 'none');
T2 = F_gridded(0.403, 0.503);
time_gridded = toc;

% Ročna bilinearna interpolacija
tic;
T3 = manualBilinearInterpolation(x, y, temperatura, 0.403, 0.503, celice_data);
time_manual = toc;

% Najdi najbližjega soseda
tic;
T_nearest = nearestNeighbor(x, y, temperatura, 0.403, 0.503);
time_nearest = toc;

% Poišči največjo temperaturo
[max_temp, idx_max] = max(temperatura);
max_coord = [x(idx_max), y(idx_max)];

% Izpiši rezultate
fprintf('ScatteredInterpolant: T = %.2f, čas = %.4f s\n', T1, time_scattered);
fprintf('GriddedInterpolant: T = %.2f, čas = %.4f s\n', T2, time_gridded);
fprintf('Bilinearna interpolacija: T = %.2f, čas = %.4f s\n', T3, time_manual);
fprintf('Najbližji sosed: T = %.2f, čas = %.4f s\n', T_nearest, time_nearest);
fprintf('Največja temperatura: %.2f pri koordinatah (%.3f, %.3f)\n', max_temp, max_coord(1), max_coord(2));

% Funkcija za pripravo mrežnih podatkov
function [x_unique, y_unique, T_grid] = pripraviMreznePodatke(x, y, temperatura)
    x_unique = unique(x);
    y_unique = unique(y);
    T_grid = nan(length(x_unique), length(y_unique)); % naredimo prazno matriko T
    for i = 1:length(x)
    
        ix = find(x_unique == x(i));
        iy = find(y_unique == y(i));
        T_grid(ix, iy) = temperatura(i); % v T matriko dodamo vrednosti ki pašejo skupaj iti x in y združimo z ito temperaturo
    end
end

% Funkcija za najbližjega soseda (dodatno
function T = nearestNeighbor(x, y, temperature, x_target, y_target)
    distances = sqrt((x - x_target).^2 + (y - y_target).^2);
    [neki, idx] = min(distances);
    T = temperature(idx);
    
end

% Funkcija za ročno bilinearno interpolacijo
function T = manualBilinearInterpolation(x, y, temperature, x_target, y_target, celice_data)
    num_cells = height(celice_data);
    T = NaN;
    for i = 1:num_cells
        cell_indices = table2array(celice_data(i, :));
        cell_x = x(cell_indices);
        cell_y = y(cell_indices);
        cell_temps = temperature(cell_indices);

        if x_target >= min(cell_x) && x_target <= max(cell_x) && ...
           y_target >= min(cell_y) && y_target <= max(cell_y)
            xmin = min(cell_x);
            xmax = max(cell_x);
            ymin = min(cell_y);
            ymax = max(cell_y);
            
            T11 = cell_temps(1);
            T21 = cell_temps(2);
            T22 = cell_temps(3);
            T12 = cell_temps(4);

            K1 = (xmax - x_target) / (xmax - xmin) * T11 + (x_target - xmin) / (xmax - xmin) * T21;
            K2 = (xmax - x_target) / (xmax - xmin) * T12 + (x_target - xmin) / (xmax - xmin) * T22;
            T = (ymax - y_target) / (ymax - ymin) * K1 + (y_target - ymin) / (ymax - ymin) * K2;
            return;
        end
    end
end
