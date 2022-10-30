%clear
%clc
%close all

SAV_num = 5;
target = EP_y(:,(T-1)*(SAV_num-1)+1:(T-1)*(SAV_num-1)+T-1);
%target = EP_y_sum;
target_sum = sum(target,2);
scale = 20;


%% Step1 各ノードのxy座標を決定
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


%% Step2 各ノード・リンクのSAV拠点からの距離を距離別に分類
%{
node_far = zeros(N,1); %ノードの遠さを表現
link_far = zeros(L,1); %リンクの遠さを表現

for i=1:N
    node_far(i,1) = abs(ox_node_grid(SAV_grid(SAV_num,2)) - ox_node_grid(i));
    node_far(i,1) = node_far(i,1) + abs(oy_node_grid(SAV_grid(SAV_num,2)) - oy_node_grid(i));
end
for i=1:L
    link_far(i,1) = node_far(road_grid(i,2)) + node_far(road_grid(i,3));
end
type_len = (max(link_far)+1)/2; %リンクの遠さの種類数
link_type = zeros(type_len,2); %リンクの遠さの種類別SAV台数
target_sum = sum(target,2);
for i=1:type_len
    link_type(i,1) = 2*i-1; %T列目にリンクの種類を表現
end
for i=1:L
    link_type((link_far(i,1)+1)/2,2) = link_type((link_far(i,1)+1)/2,2) + target_sum(i);
end
%}



%% Step3 図にプロット
clf
hold on
for j=-(grid_num-1)/2:1:(grid_num-1)/2
    hline = refline([0,j]);
    hline.Color = [220 220 220]/255;
    hline = xline(j);
    hline.Color = [200 200 200]/255;
    pbaspect([1 1 1])
    xlim([-(grid_num-1)/2 (grid_num-1)/2])
    ylim([-(grid_num-1)/2 (grid_num-1)/2])
end
%c=gray(20);
for i=1:L/2
    if target_sum(2*i-1)+target_sum(2*i)>0.001 %これが0だとエラー
        plot([ox_node_grid(road_grid(i*2,2)),ox_node_grid(road_grid(i*2,3))],...
            [oy_node_grid(road_grid(i*2,2)),oy_node_grid(road_grid(i*2,3))],...
            'r-', 'LineWidth', (target_sum(2*i-1)+target_sum(2*i))/scale)
            %'r-', 'LineWidth', link_type((link_far(i*2,1)+1)/2,2)/scale)
        xlim([-(grid_num-1)/2 (grid_num-1)/2])
        ylim([-(grid_num-1)/2 (grid_num-1)/2])
        hold on
    end
    figure_list(i) = target_sum(2*i-1)+target_sum(2*i);
end

plot(ox_node_grid(SAV_grid(SAV_num,2)), oy_node_grid(SAV_grid(SAV_num,2)),...
    '-s','MarkerSize',18,'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
xlim([-(grid_num-1)/2 (grid_num-1)/2])
ylim([-(grid_num-1)/2 (grid_num-1)/2])
hold on

for i=1:4
    plot(ox_node_grid(corner_in(i)), oy_node_grid(corner_in(i)),...
        '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'green',...
        'MarkerFaceColor', [.4 1 .4])
    pbaspect([1 1 1])
    xlim([-(grid_num-1)/2 (grid_num-1)/2])
    ylim([-(grid_num-1)/2 (grid_num-1)/2])
end

for i=1:4
    plot(ox_node_grid(corner_out(i)), oy_node_grid(corner_out(i)),...
        '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'red',...
        'MarkerFaceColor', [1 .6 .6])
    pbaspect([1 1 1])
    xlim([-(grid_num-1)/2 (grid_num-1)/2])
    ylim([-(grid_num-1)/2 (grid_num-1)/2])
end

% plot(ox_node_grid((grid_num^2+1)/2), oy_node_grid((grid_num^2+1)/2),...
%     '-o', 'MarkerSize', 7, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', [.8 .8 .8])
% pbaspect([1 1 1])
% xlim([-(grid_num-1)/2 (grid_num-1)/2])
% ylim([-(grid_num-1)/2 (grid_num-1)/2])

grid off
% for j=-(grid_num-1)/2:1:(grid_num-1)/2
%     hline = refline([0,j]);
%     hline.Color = [220 220 220]/255;
%     hline = xline(j);
%     hline.Color = [200 200 200]/255;
%     pbaspect([1 1 1])
%     xlim([-(grid_num-1)/2 (grid_num-1)/2])
%     ylim([-(grid_num-1)/2 (grid_num-1)/2])
% end
xlabel('x')
ylabel('y')

