%clear
clc
close all

EP_x_sum_bi = EP_x_sum;
EP_x_sum_bi(EP_x_sum_bi>0.01) = 1;
EP_p = EP_p.*EP_x_sum_bi;
EP_p(EP_p<0.01) = 0;

SAV_num = 4;
target = EP_y(:,(T-1)*(SAV_num-1)+1:(T-1)*(SAV_num-1)+T-1);
target = EP_x_sum;
%target = EP_y_sum*4-EP_x_sum;
%target = EP_e;
scale = 5;

%% Step1 各ノードの座標を決定
node_grid = N;
ox_node_grid = zeros(1,node_grid);
oy_node_grid = zeros(1,node_grid);

% 左上から数字をつける
num = 0;
for i=1:grid_num
    for j=1:grid_num
        num = num+1;
        ox_node_grid(num) = -(grid_num+1)/2+j;
        oy_node_grid(num) = (grid_num+1)/2-i;
    end
end

% 中心から数字をつける
%{
for i=1:5
    if i==1
        ox_node_grid(i) = 0;
        oy_node_grid(i) = 0;
    else
        for j=1:(i*2-2)*4
            if j>0 && j<=i*2-2
                ox_node_grid(((i-1)*2-1)^2+j) = -(i-1)+mod(j-1,i*2-2);
                oy_node_grid(((i-1)*2-1)^2+j) = i-1;
            elseif j>i*2-2 && j<=(i*2-2)*2
                ox_node_grid(((i-1)*2-1)^2+j) = i-1;
                oy_node_grid(((i-1)*2-1)^2+j) = i-1-mod(j-1,i*2-2);
            elseif j>(i*2-2)*2 && j<=(i*2-2)*3
                ox_node_grid(((i-1)*2-1)^2+j) = i-1-mod(j-1,i*2-2);
                oy_node_grid(((i-1)*2-1)^2+j) = -(i-1);
            else
                ox_node_grid(((i-1)*2-1)^2+j) = -(i-1);
                oy_node_grid(((i-1)*2-1)^2+j) = -(i-1)+mod(j-1,i*2-2);
            end
        end
    end
end
%}


%% Step2
figh = figure;
xlim([-4 4])
ylim([-4 4])

for i=1:T-1
    %Wipe the slate clean so we are plotting with a black figure
    clf
    
    %Plot the current location of the particle
    for j=1:L/2
        if target(j*2,i) > target(j*2-1,i)
            if target(j*2,i)>0.001 %緑（左向き）
                plot([ox_node_grid(road_grid(j*2,2)),ox_node_grid(road_grid(j*2,3))],[oy_node_grid(road_grid(j*2,2)),oy_node_grid(road_grid(j*2,3))], 'g-', 'LineWidth', target(j*2,i)/scale)
                xlim([-(grid_num-1)/2 (grid_num-1)/2])
                ylim([-(grid_num-1)/2 (grid_num-1)/2])
                hold on
            end
            if target(j*2-1,i)>0.001 %赤（右向き）
                plot([ox_node_grid(road_grid(j*2-1,2)),ox_node_grid(road_grid(j*2-1,3))],[oy_node_grid(road_grid(j*2-1,2)),oy_node_grid(road_grid(j*2-1,3))], 'r-', 'LineWidth', target(j*2-1,i)/scale)
                xlim([-(grid_num-1)/2 (grid_num-1)/2])
                ylim([-(grid_num-1)/2 (grid_num-1)/2])
                hold on
            end
        else
            if target(j*2-1,i)>0.001 %赤（下向き）
                plot([ox_node_grid(road_grid(j*2-1,2)),ox_node_grid(road_grid(j*2-1,3))],[oy_node_grid(road_grid(j*2-1,2)),oy_node_grid(road_grid(j*2-1,3))], 'r-', 'LineWidth', target(j*2-1,i)/scale)
                xlim([-(grid_num-1)/2 (grid_num-1)/2])
                ylim([-(grid_num-1)/2 (grid_num-1)/2])
                hold on
            end
            if target(j*2,i)>0.001 %緑（上向き）
                plot([ox_node_grid(road_grid(j*2,2)),ox_node_grid(road_grid(j*2,3))],[oy_node_grid(road_grid(j*2,2)),oy_node_grid(road_grid(j*2,3))], 'g-', 'LineWidth', target(j*2,i)/scale)
                xlim([-(grid_num-1)/2 (grid_num-1)/2])
                ylim([-(grid_num-1)/2 (grid_num-1)/2])
                hold on
            end
            if target(j*2,i)==0 && target(j*2-1,i)==0
                xlim([-(grid_num-1)/2 (grid_num-1)/2])
                ylim([-(grid_num-1)/2 (grid_num-1)/2])
                hold on
            end
        end
    end
    hold off

    
    grid on
    xlabel('x')
    ylabel('y')
    title(['Particle at t = ', num2str(i), ' seconds'])
    %view([30 35])
    
    movieVector(i) = getframe(figh, [10 10 520 400]);
end

writeAnimation('GRID.gif')

myWriter = VideoWriter('GRID', 'MPEG-4');
myWriter.FrameRate = 5;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}

