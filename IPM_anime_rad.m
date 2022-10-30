%clear
clc
close all

EP_x_sum_bi = EP_x_sum;
EP_x_sum_bi(EP_x_sum_bi>0.01) = 1;
EP_p = EP_p.*EP_x_sum_bi;
EP_p(EP_p<0.01) = 0;

target = EP_x_sum;
%target = EP_y(:,(T-1)*24+1:(T-1)*24+T-1);
%target = EP_e;
scale = 100;

%% Step1 各ノードの座標を決定
node_rad = N;
ox_node_rad = zeros(1,node_rad);
oy_node_rad = zeros(1,node_rad);

% 中心から数字(座標)をつける
num = 0;
for i=1:rad_num+1
    if i==1
        num=num+1;
        ox_node_rad(num) = 0;
        oy_node_rad(num) = 0;
    else
        for j=1:(i-1)*6
            if j==1
                num=num+1;
                ox_node_rad(num) = 0;
                oy_node_rad(num) = i-1;
            elseif j>1 && j<=(i-1)*1+1
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1)+cos(11*pi/6);
                oy_node_rad(num) = oy_node_rad(num-1)+sin(11*pi/6);
            elseif j>(i-1)*1+1 && j<=(i-1)*2+1
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1);
                oy_node_rad(num) = oy_node_rad(num-1)-1;
            elseif j>(i-1)*2+1 && j<=(i-1)*3+1
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1)+cos(7*pi/6);
                oy_node_rad(num) = oy_node_rad(num-1)+sin(7*pi/6);
            elseif j>(i-1)*3+1 && j<=(i-1)*4+1
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1)+cos(5*pi/6);
                oy_node_rad(num) = oy_node_rad(num-1)+sin(5*pi/6);
            elseif j>(i-1)*4+1 && j<=(i-1)*5+1
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1);
                oy_node_rad(num) = oy_node_rad(num-1)+1;
            else
                num=num+1;
                ox_node_rad(num) = ox_node_rad(num-1)+cos(pi/6);
                oy_node_rad(num) = oy_node_rad(num-1)+sin(pi/6);
            end
        end
    end
end
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


%% Step2 時間ごとに変化
figh = figure;
pbaspect([1 1 1])
xlim([-4 4])
ylim([-4 4])

for i=1:T-1
    %Wipe the slate clean so we are plotting with a black figure
    clf
    
    %Plot the current location of the particle
    for j=1:L/2
        if target(j*2,i) > target(j*2-1,i)
            if target(j*2,i)>0.001 %外側→内側、反時計回り（緑）
                plot([ox_node_rad(road_rad(j*2,2)),ox_node_rad(road_rad(j*2,3))],[oy_node_rad(road_rad(j*2,2)),oy_node_rad(road_rad(j*2,3))], 'g-', 'LineWidth', target(j*2,i)/scale)
                pbaspect([1 1 1])
                xlim([-4 4])
                ylim([-4 4])
                hold on
            end
            if target(j*2-1,i)>0.001 %内側→外側、時計回り（赤）
                plot([ox_node_rad(road_rad(j*2-1,2)),ox_node_rad(road_rad(j*2-1,3))],[oy_node_rad(road_rad(j*2-1,2)),oy_node_rad(road_rad(j*2-1,3))], 'r-', 'LineWidth', target(j*2-1,i)/scale)
                pbaspect([1 1 1])
                xlim([-4 4])
                ylim([-4 4])
                hold on
            end
        else
            if target(j*2-1,i)>0.001 %外側→内側、反時計回り（緑）
                plot([ox_node_rad(road_rad(j*2-1,2)),ox_node_rad(road_rad(j*2-1,3))],[oy_node_rad(road_rad(j*2-1,2)),oy_node_rad(road_rad(j*2-1,3))], 'r-', 'LineWidth', target(j*2-1,i)/scale)
                pbaspect([1 1 1])
                xlim([-4 4])
                ylim([-4 4])
                hold on
            end
            if target(j*2,i)>0.001 %外側→内側、反時計回り（緑）
                plot([ox_node_rad(road_rad(j*2,2)),ox_node_rad(road_rad(j*2,3))],[oy_node_rad(road_rad(j*2,2)),oy_node_rad(road_rad(j*2,3))], 'g-', 'LineWidth', target(j*2,i)/scale)
                pbaspect([1 1 1])
                xlim([-4 4])
                ylim([-4 4])
                hold on
            end
            if target(j*2,i)==0 && target(j*2-1,i)==0
                pbaspect([1 1 1])
                xlim([-4 4])
                ylim([-4 4])
                hold on
            end
        end
    end
    hold off

    
    grid off
    for j=-4:1:4
        hline = refline([tan(pi/6),j]);
        hline.Color = [0.5 0.5 0.5];
        hline = refline([tan(5*pi/6),j]);
        hline.Color = [0.5 0.5 0.5];
        hline = xline(0+j*cos(pi/6));
        hline.Color = [17 17 17]/255;
        pbaspect([1 1 1])
        xlim([-4 4])
        ylim([-4 4])
    end
    xlabel('x')
    ylabel('y')
    title(['Particle at t = ', num2str(i), ' seconds'])
    %view([30 35])
    
    movieVector(i) = getframe(figh, [10 10 520 400]);
end

writeAnimation('RAD3.gif')

myWriter = VideoWriter('RAD3', 'MPEG-4');
myWriter.FrameRate = 5;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}
