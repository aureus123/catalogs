clear
clc

# Read matrix with magnitudes
A = csvread('results/matrices.csv');
[m n] = size(A);
fprintf("Number of original rows = %d\n", m)

# Filter by number
num = 1;

row = 0;
for i = 1:m
    if A(i,1) == num
        row = row + 1;
        B(row,1) = A(i,2);
        B(row,2) = A(i,3);
    end
end
fprintf("Number of filtered rows (Tomo %d) = %d\n", num, row)

# Perform quadratic fit
F = [ones(1, row) ; B(:,1)' ; (B(:,1).*B(:,1))' ];
M = F*F';
y = F*B(:,2);
coef = M\y
ECM = sqrt((1/row) * sum((coef(1) + coef(2)*B(:,1) + coef(3)*B(:,1).*B(:,1) - B(:,2)).^2))

# Show fit
plot(B(:,1), B(:,2), 'r+')
hold on
xx = 3:0.1:11;
plot(xx, coef(1) + coef(2)*xx + coef(3)*xx.*xx, 'b')