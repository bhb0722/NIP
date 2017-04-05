
train_d = csvread('train.csv');
test_d = csvread('test.csv');
train_x = train_d(:,1);
train_y = train_d(:, 2);
test_x = test_d(:,1);
test_y = test_d(:, 2);


figure
plot(train_x);

hold on;

plot(train_y);

figure
plot(test_x);

hold on;

plot(test_y);