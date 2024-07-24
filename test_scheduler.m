test_num = 1;
savename = 'test_task_scheduler.mat';
save(savename, "test_num");

figure;
for n = [0.5, 1, 2, 3]
    x = 0: 0.1: 15;
    y = 1./(1+(15-x).^n);
    plot(x, y, 'DisplayName', ['n = ', num2str(n)]);
    hold on;
end
hold off;
legend();


figure;
variance = [0.29452/0.14108, 0.31744/0.22344, 0.79219/0.16804];
figure; plot([1, 1, 2], variance, 'o');


% plot the function f = 2xy/(x+y)
figure;
x = 0:0.1:1;
y = 0:0.1:1;
[X,Y] = meshgrid(x, y);
f = 2*X.*Y./(X+Y);
surf(X, Y, f);
xlabel('pturn(blue)'); ylabel('pturn(red)'); zlabel('pturn(blue+red)');
