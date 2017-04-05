

A = csvread('sm_1.csv');
E = [2:10];
Tau = [1:24];

E = E';
Tau = Tau';
Sm = A(:, 3);
Smm = zeros(24, 9);
for i = 1 : 9
    for j = 1 : 24
        Smm(j, i) = Sm(7*(i-1) + j);
    end
end
            
            
            
[x y] = meshgrid(E, Tau);
mesh(x, y, Smm)