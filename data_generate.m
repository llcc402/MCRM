%% data generation process
%
%     mu_0 = 1/2 N(0,1) + 1/2 N(5,1)
%     mu_1 = 1/2 N(-5,1) + 1/2 N(-10,1)
%     mu_2 = 1/2 N(10,1) + 1/2 N(15,1)
% 
%     p_1 = 0.1 * mu_0 + 0.2 * mu_1 + 0.7 * mu_3
%     p_2 = 0.7 * mu_0 + 0.2 * mu_2 + 0.1 * mu_3

function data = data_generate()

x1 = randn(1, 100);
x2 = randn(1, 100) + 5;
x3 = randn(1, 100) - 5;
x4 = randn(1, 100) - 10;
x5 = randn(1, 100) + 10;
x6 = randn(1, 100) + 15;

data = zeros(2, 400);
data(1,:) = [x1, x2, x3, x4];
data(2,:) = [x1, x2, x3, x4];

figure(1)
hist(data(1,:), 50)
title('The first group')

figure(2)
hist(data(2,:), 50)
title('The second group')
