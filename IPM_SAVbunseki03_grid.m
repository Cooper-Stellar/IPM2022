%clear
%clc
%close all

SAV_num = 2;
target = EP_y(:,(T-1)*(SAV_num-1)+1:(T-1)*(SAV_num-1)+T-1);
%target = EP_y_sum;
target_sum = sum(target,2);
scale = 50;


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
%%{
node_far = zeros(N,2); %ノードの遠さを表現
link_far = zeros(L,5); %リンクの遠さを表現（リンク番号，遠さ，近づいている(0)or離れている(1)，SAV量，付番（並び替え後））

for i=1:N
    node_far(i,1) = i;
    node_far(i,2) = abs(ox_node_grid(SAV_grid(SAV_num,2)) - ox_node_grid(i));
    node_far(i,2) = node_far(i,2) + abs(oy_node_grid(SAV_grid(SAV_num,2)) - oy_node_grid(i));
end
for i=1:L
    link_far(i,1) = i;
    link_far(i,2) = node_far(road_grid(i,2),2) + node_far(road_grid(i,3),2);
    if node_far(road_grid(i,2),2)<node_far(road_grid(i,3),2)
        link_far(i,3) = 1;
    end
    link_far(i,4) = target_sum(i);
end
type_len = (max(link_far(:,2))+1)/2; %リンクの遠さの種類数
link_type = zeros(type_len,2); %リンクの遠さの種類別SAV台数
for i=1:type_len
    link_type(i,1) = 2*i-1; %T列目にリンクの種類を表現
end
for i=1:L
    link_type((link_far(i,2)+1)/2,2) = link_type((link_far(i,2)+1)/2,2) + target_sum(i);
end

link_far_near = topkrows(link_far,L,2:3,{'ascend','ascend'});
for i=1:L
    link_far_near(i,5) = i;
end
%%}



%% Step3 図にプロット
link_scatter = zeros(L*(T-1),5); %（x座標, y座標, SAV量, リンク番号, 近づいている(0)or離れている(1), 遠さ）
link_shuffle = zeros(L,2+(T-1));
link_shuffle(:,T) = link_far_near(:,1);
link_shuffle(:,T+1) = link_far_near(:,5);
link_shuffle = topkrows(link_shuffle,L,T,{'ascend'});
link_shuffle(:,1:T-1) = target;
link_shuffle = topkrows(link_shuffle,L,T+1,{'ascend'});

for i=1:T-1
    link_scatter(L*(i-1)+1:L*i,1) = link_far_near(:,5);
    link_scatter(L*(i-1)+1:L*i,2) = i;
    link_scatter(L*(i-1)+1:L*i,3) = link_shuffle(:,i);
    link_scatter(L*(i-1)+1:L*i,4) = link_far_near(:,1);
    link_scatter(L*(i-1)+1:L*i,5) = link_far_near(:,3);
    link_scatter(L*(i-1)+1:L*i,6) = link_far_near(:,2);
end

link_scatter_a = link_scatter;
num=0;
for i=1:L*(T-1)
    if link_scatter_a(i-num,3)<0.00001 || link_scatter_a(i-num,5)==1
        link_scatter_a(i-num,:) = [];
        num = num+1;
    end
end

link_scatter_b = link_scatter;
num=0;
for i=1:L*(T-1)
    if link_scatter_b(i-num,3)<0.00001 || link_scatter_b(i-num,5)==0
        link_scatter_b(i-num,:) = [];
        num = num+1;
    end
end

clf
scatter(link_scatter_a(:,1),link_scatter_a(:,2),link_scatter_a(:,3)*10,'red')
hold on
scatter(link_scatter_b(:,1),link_scatter_b(:,2),link_scatter_b(:,3)*10,'blue')

xline(0,'-.',{['遠さ ',num2str(1)]},'FontSize',13)
for i=1:L*(T-1)-1
    if link_scatter(i,6) < link_scatter(i+1,6)
        xline(link_scatter(i,1),'-.',{['遠さ ',num2str((link_scatter(i,6)+1)/2+1)]},'FontSize',13)
    end
end
yline(D_hope(1),'-.g','LineWidth',1.3)
%xlim([0 180])
ylim([0 25])
hold off
            
