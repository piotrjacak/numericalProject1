function [wynik, errest, m] = P1Z37_PJA_gaussLegendre(f, a, b, tol, mMax)
% Projekt 1, zadanie 37
% Piotr Jacak, 327354
%
% Obliczanie całki z funkcji f na przedziale [a,b] metodą 2-punktowej
% złożonej kwadratury Gaussa-Legendre'a z zadaną dokładnością (funkcja
% dzieli przedział [a,b] na coraz mniejsze przedziały i na każdym z nich
% stosuje powyższą metodę, aż do otrzymania zadanej dokładności lub 
% przekroczenia zadanej liczby iteracji)
% WEJŚCIE:
%   f - uchwyt do funkcji podcałkowej f
%   a - początek przedziału całkowania (a należy do R)
%   b - koniec przedziału całkowania (b należy do R)
%   tol - pożądana dokładność obliczeń (domyślnie 1e-10)
%   mMax - maksymalna ilość iteracji pętli w funkcji (domyślnie 18)
% WYJŚCIE:
%   wynik - wyznaczona przybliżona wartość liczbowa całki
%   errest - ostatnia uzyskana wartość bezwzględna różnicy kolejnych
%       przybliżeń wartości całki
%   m - końcowa ilość podprzedziałów koniecznych do uzyskania zadanej
%       dokładności

% Zapisanie do zmiennej m początkowej ilości podprzedziałów
m = 1;
% Inicjalizacja parametrów wyjściowych
wynik = 0;
errest = NaN;
% Zmienne pomocnicze
blad = Inf;

% Sprawdzenie poprawności podanych argumentów i ustawienie argumentu
% domyślnego
if nargin < 3 || b < a
    wynik = NaN;
    return;
elseif nargin == 3
    tol = 1e-10;
    mMax = 18;
elseif nargin == 4
    mMax = 18;
end

while blad > tol
    % Zapisanie do wektora p wszystkich krańców m podprzedziałów
    podzial = linspace(a, b, m+1);
    % Obliczanie wartości całki dla danych m podprzedziałów
    suma = gauss_legendre(f, podzial, m);
    % Obliczenie błędu kolejnych przybliżeń i nadpisanie zmiennej wynik
    blad = abs(wynik - suma);
    wynik = suma;
    % Sprawdzenie czy pętla nie wykonała za dużo iteracji (jeśli funkcja
    % nie może osiągnąć zadanej dokładności lub całka nie jest zbieżna)
    if m >= 2^mMax
        errest = NaN;
        return;
    end
    % Podwojenie liczby przedziałów
    m = m*2;
end

errest = blad;
m = m/2;

end % function
