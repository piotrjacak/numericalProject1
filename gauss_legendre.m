function [suma] = gauss_legendre(f, p, m)
% Projekt 1, zadanie 37
% Piotr Jacak, 327354

% Funkcja pomocnicza obliczająca przybliżoną całkę z funkcji f stosując
% na m przedziałach metodę 2-punktowej kwadratury Gaussa-Legendre'a
% WEJŚCIE:
%   f - uchwyt do funkcji podcałkowej
%   p - (m+1) - elementowy wektor z podziałem danego w zadaniu przedziału 
%       na m podprzedziałow
%   m - ilość pożądanych podprzedziałów
% WYJŚCIE:
%   suma - suma przybliżonych wartości całek na m przedziałach

% Zapisanie do wektorów w i x odpowiednio wag i węzłów dla kwadratury
% Gaussa-Legendre'a
w = [1, 1];
x = [-sqrt(1/3), sqrt(1/3)];

suma = 0;

for i = 1: m
    % Stworzenie wektora x_transform z przekształconymi węzłami x, tak aby
    % pasowały do danego podprzedziału
    x_transform = ((p(i+1) - p(i)) .* x) / 2 + (p(i) + p(i+1)) / 2;
    % Obliczenie przybliżonej całki dla podprzedziału i sumowanie
    suma = suma + ((p(i+1) - p(i))/2) * ... 
            (w(1) * f(x_transform(1)) + w(2) * f(x_transform(2)));
end

end % function

