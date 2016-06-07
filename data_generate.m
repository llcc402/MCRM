%% data generation process
%
%     mu_0 = 1/2 N(0,1) + 1/2 N(5,1)
%     mu_1 = 1/2 N(-5,1) + 1/2 N(-10,1)
%     mu_2 = 1/2 N(10,1) + 1/2 N(15,1)
% 
%     p_1 = 0.2 * mu_0 + 0.3 * mu_1 + 0.5 * mu_2
%     p_2 = 0.5 * mu_0 + 0.3 * mu_1 + 0.2 * mu_2

function data = data_generate()

data = zeros(2, 300);
for i = 1:300
    u = rand(1);
    if u < 0.2
        u1 = rand(1);
        if u1 < 0.5
            data(1,i) = randn(1);
        else
            data(1,i) = randn(1) + 5;
        end
    elseif u < 0.5
        u1 = rand(1);
        if u1 < 0.5
            data(1,i) = randn(1) - 5;
        else
            data(1,i) = randn(1) -10;
        end
    else
        u1 = rand(1);
        if u1 < 0.5
            data(1,i) = randn(1) + 10;
        else
            data(1,i) = randn(1) + 15;
        end
    end   
end

for i = 1:300
    u = rand(1);
    if u < 0.5
        u1 = rand(1);
        if u1 < 0.5
            data(2,i) = randn(1);
        else
            data(2,i) = randn(1) + 5;
        end
    elseif u < 0.8
        u1 = rand(1);
        if u1 < 0.5
            data(2,i) = randn(1) - 5;
        else
            data(2,i) = randn(1) -10;
        end
    else
        u1 = rand(1);
        if u1 < 0.5
            data(2,i) = randn(1) + 10;
        else
            data(2,i) = randn(1) + 15;
        end
    end   
end

figure(1)
hist(data(1,:), 50)
title('The first group')

figure(2)
hist(data(2,:), 50)
title('The second group')
end
            

