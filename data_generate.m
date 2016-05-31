%% data generation process
%
%     mu_0 = 1/2 N(0,1) + 1/2 N(5,1)
%     mu_1 = 1/2 N(-5,1) + 1/2 N(-10,1)
%     mu_2 = 1/2 N(10,1) + 1/2 N(15,1)
% 
%     p_1 = mu_0 + mu_1
%     p_2 = mu_0 + mu_2

function data = data_generate()

x1 = randn(1, 200);
x2 = randn(1, 200) + 5;
