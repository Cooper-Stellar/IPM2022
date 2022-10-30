%% スイッチ
which_SAVarrange = 0; %0=全体，1=外側


%% 与件のパラメータ
Q = 144*105; %総利用者数
%Q = 144*275; %総利用者数


%% 直行(グリッド)型ネットワーク①
% ネットワーク情報(リンク番号，発ノード，着ノード，自由旅行時間，容量)
D_hope = 15; %希望到着時間 %デフォルト15 %改18

grid_num = 5; %1辺に存在するノード数(グリッドネットワークの大きさ) %デフォルト9
free_travel_time = 1; %自由旅行時間 %デフォルト1
capacity_grid = 30; %容量 %デフォルト30
%capacity = Inf; %容量
road_grid = zeros((grid_num-1)*grid_num*4,5);
road_grid(:,4) = free_travel_time;
road_grid(:,5) = capacity_grid;
num=0;
for i=1:grid_num-1 
    for j=1:grid_num-1 
        num = num+1;
        road_grid(num,1:3) = [num,grid_num*(i-1)+j,grid_num*(i-1)+j+1]; %左→右
        num = num+1;
        road_grid(num,1:3) = [num,grid_num*(i-1)+j+1,grid_num*(i-1)+j]; %右→左
    end
    for j=1:grid_num
        num = num+1;
        road_grid(num,1:3) = [num,grid_num*(i-1)+j,grid_num*i+j]; %上→下
        num = num+1;
        road_grid(num,1:3) = [num,grid_num*i+j,grid_num*(i-1)+j]; %下→上
    end
end
for j=1:grid_num-1 
    num = num+1;
    road_grid(num,1:3) = [num,grid_num*(grid_num-1)+j,grid_num*(grid_num-1)+j+1]; %左→右
    num = num+1;
    road_grid(num,1:3) = [num,grid_num*(grid_num-1)+j+1,grid_num*(grid_num-1)+j]; %下→上
end

road_grid(:,4) = free_travel_time;
road_grid(:,5) = capacity_grid;

%%{
enlarge = 100; %デフォルト100

road_grid(1:2,5) = capacity_grid*enlarge;
road_grid((grid_num-1)*2-1:(grid_num-1)*2+2,5) = capacity_grid*enlarge;
road_grid((grid_num*2-1)*2-1:(grid_num*2-1)*2,5) = capacity_grid*enlarge;

road_grid((grid_num*2-1)*2*(grid_num-2)+(grid_num-1)*2+1:(grid_num*2-1)*2*(grid_num-2)+(grid_num-1)*2+2,5) = capacity_grid*enlarge; %(255:256)
road_grid((grid_num*2-1)*2*(grid_num-1)-1:(grid_num*2-1)*2*(grid_num-1)+2,5) = capacity_grid*enlarge; %(271:274)
road_grid((grid_num*2-1)*2*(grid_num-1)+(grid_num-1)*2-1:(grid_num*2-1)*2*(grid_num-1)+(grid_num-1)*2,5) = capacity_grid*enlarge; %(287:288)
%%}


% OD(利用者)情報(出発ノード，到着ノード，希望到着時間，利用者数)
user_grid_one = zeros(4*4,4); %一方向OD表
user_grid_bi = zeros(4*4*2,4); %双方向OD表
corner_out = [1, grid_num, grid_num*(grid_num-1)+1, grid_num^2];
corner_in = zeros(1,4);
corner_in(1) = grid_num*((grid_num-1)/2-1)+(grid_num+1)/2;
corner_in(2) = grid_num*((grid_num-1)/2)+(grid_num-1)/2;
corner_in(3) = corner_in(2)+2;
corner_in(4) = grid_num*((grid_num-1)/2+1)+(grid_num+1)/2;
for i=1:4
    for j=1:4
        user_grid_one((i-1)*4+j,:) = [corner_out(1,i), corner_in(1,j), D_hope, Q/(4*4)];
    end
end
for i=1:4
    for j=1:4
        user_grid_bi((i-1)*4+j,:) = [corner_out(1,i), corner_in(1,j), D_hope, Q/(4*4*2)];
        user_grid_bi((3+j)*4+i,:) = [corner_in(1,j), corner_out(1,i), D_hope, Q/(4*4*2)];
    end
end


% SAVプロバイダー情報(SAV番号，拠点)
if which_SAVarrange==0
    %SAV_grid = zeros(4*4+1,2);
    SAV_grid = [
        1 1
        2 (grid_num+1)/2
        3 grid_num
        4 21 %21%7
        5 25 %25%9
        6 (grid_num^2+1)/2-grid_num
        7 (grid_num^2+1)/2-(grid_num-1)/2
        8 (grid_num^2+1)/2-1
        9 (grid_num^2+1)/2
        10 (grid_num^2+1)/2+1
        11 (grid_num^2+1)/2+(grid_num-1)/2
        12 (grid_num^2+1)/2+grid_num
        13 57 %57%17
        14 61 %61%19
        15 grid_num^2-(grid_num-1)
        16 grid_num^2-(grid_num-1)/2
        17 grid_num^2
        ];
    %%{
    SAV_grid = [
        1 1
        2 grid_num
        3 grid_num^2-(grid_num-1)
        4 grid_num^2
        5 (grid_num^2+1)/2-grid_num
        6 (grid_num^2+1)/2-1
        7 (grid_num^2+1)/2+1
        8 (grid_num^2+1)/2+grid_num
        ];
    %%}
    if grid_num<6
        SAV_grid = zeros(grid_num^2,2);
        for i=1:grid_num^2
            SAV_grid(i,:) = [i, i];
        end
    end
elseif which_SAVarrange==1
    SAV_grid = [
        1 1
        2 grid_num
        3 grid_num^2-(grid_num-1)
        4 grid_num^2
        ];
end



%% 直行(グリッド)型ネットワーク②
% ネットワーク情報(リンク番号，発ノード，着ノード，自由旅行時間，容量)

%{
grid_num = 9; %1辺に存在するノード数(グリッドネットワークの大きさ)
free_travel_time = 1; %自由旅行時間
capacity = 30; %容量
road_grid = zeros((grid_num-1)*grid_num*4,5);
road_grid(:,4) = free_travel_time;
road_grid(:,5) = capacity;
num=0;
if mod(grid_num,2) == 1 %奇数
    for i=1:grid_num-1 
        if mod(i,2) == 1 %連結リンク(i=1辺のリンク数)
            if i==1
                for k=1:4
                    num = num+1;
                    road_grid(num,1:3) = [num,1,k*2+1]; %内側→外側
                    num = num+1;
                    road_grid(num,1:3) = [num,k*2+1,1]; %外側→内側
                end
            else
                for j=1:4
                    if j==1 %上辺
                        for k=1:i
                            num = num+1;
                            road_grid(num,1:3) = [num,(i-2)^2+k,i^2+k+1]; %内側→外側
                            num = num+1;
                            road_grid(num,1:3) = [num,i^2+k+1,(i-2)^2+k]; %外側→内側
                        end
                    end
                    if j==2 %右辺
                        for k=1:i
                            num = num+1;
                            road_grid(num,1:3) = [num,(i-2)^2+(i-1)*1+k,i^2+(i+1)*1+1+k];
                            num = num+1;
                            road_grid(num,1:3) = [num,i^2+(i+1)*1+1+k,(i-2)^2+(i-1)*1+k];
                        end
                    end
                    if j==3 %下辺
                        for k=1:i
                            num = num+1;
                            road_grid(num,1:3) = [num,(i-2)^2+(i-1)*2+k,i^2+(i+1)*2+1+k];
                            num = num+1;
                            road_grid(num,1:3) = [num,i^2+(i+1)*2+1+k,(i-2)^2+(i-1)*2+k];
                        end
                    end
                    if j==4 %下辺
                        for k=1:i-1
                            num = num+1;
                            road_grid(num,1:3) = [num,(i-2)^2+(i-1)*3+k,i^2+(i+1)*3+1+k];
                            num = num+1;
                            road_grid(num,1:3) = [num,i^2+(i+1)*3+1+k,(i-2)^2+(i-1)*3+k];
                        end
                        num = num+1;
                        road_grid(num,1:3) = [num,(i-2)^2+1,i^2+(i+1)*3+1+i];
                        num = num+1;
                        road_grid(num,1:3) = [num,i^2+(i+1)*3+1+i,(i-2)^2+1];
                    end
                end
            end
        end
        if mod(i,2) == 0 %外縁リンク(i=1辺のリンク数)
            for j=1:i*4-1
                num=num+1;
                road_grid(num,1:3) = [num,(i-1)^2+j,(i-1)^2+j+1]; %時計回り
                num=num+1;
                road_grid(num,1:3) = [num,(i-1)^2+j+1,(i-1)^2+j]; %反時計回り
            end
            num=num+1;
            road_grid(num,1:3) = [num,(i-1)^2+i*4,(i-1)^2+1];
            num=num+1;
            road_grid(num,1:3) = [num,(i-1)^2+1,(i-1)^2+i*4];
        end
    end
end
if mod(grid_num,2) == 0 %偶数
    for i=1:grid_num-1 
        if mod(i,2) == 0 %連結リンク(i=1辺のリンク数)
            for j=1:4
                if j==1 %上辺
                    for k=1:i
                        num = num+1;
                        road_grid(num,1:3) = [num,(i-2)^2+k,i^2+k+1]; %外側→内側
                        num = num+1;
                        road_grid(num,1:3) = [num,i^2+k+1,(i-2)^2+k]; %内側→外側
                    end
                end
                if j==2 %右辺
                    for k=1:i
                        num = num+1;
                        road_grid(num,1:3) = [num,(i-2)^2+(i-1)*1+k,i^2+(i+1)*1+1+k];
                        num = num+1;
                        road_grid(num,1:3) = [num,i^2+(i+1)*1+1+k,(i-2)^2+(i-1)*1+k];
                    end
                end
                if j==3 %下辺
                    for k=1:i
                        num = num+1;
                        road_grid(num,1:3) = [num,(i-2)^2+(i-1)*2+k,i^2+(i+1)*2+1+k];
                        num = num+1;
                        road_grid(num,1:3) = [num,i^2+(i+1)*2+1+k,(i-2)^2+(i-1)*2+k];
                    end
                end
                if j==4 %下辺
                    for k=1:i
                        num = num+1;
                        road_grid(num,1:3) = [num,(i-2)^2+(i-1)*3+k,i^2+(i+1)*3+1+k];
                        num = num+1;
                        road_grid(num,1:3) = [num,i^2+(i+1)*3+1+k,(i-2)^2+(i-1)*3+k];
                    end
                end
            end
        end
        if mod(i,2) == 1 %外縁リンク(i=1辺のリンク数)
            for j=1:i*4
                num=num+1;
                road_grid(num,1:3) = [num,(i-1)^2+j,(i-1)^2+j+1]; %時計回り
                num=num+1;
                road_grid(num,1:3) = [num,(i-1)^2+j+1,(i-1)^2+j]; %反時計回り
            end
        end
    end
end
road_grid(:,4) = free_travel_time;
road_grid(:,5) = capacity;


% OD(利用者)情報(出発ノード，到着ノード，希望到着時間，利用者数)
user_grid_one = zeros(4*4,4); %一方向OD表
user_grid_bi = zeros(4*4*2,4); %双方向OD表
corner_out = [(grid_num-2)^2+1, (grid_num-2)^2+grid_num, (grid_num-2)^2+(grid_num-1)*2+1, (grid_num-2)^2+(grid_num-1)*3+1];
corner_in = [2,4,6,8];
for i=1:4
    for j=1:4
        user_grid_one((i-1)*4+j,:) = [corner_out(1,i), corner_in(1,j), D_hope, Q/(4*4)];
    end
end
for i=1:4
    for j=1:4
        user_grid_bi((i-1)*4+j,:) = [corner_out(1,i), corner_in(1,j), D_hope, Q/(4*4*2)];
        user_grid_bi((3+j)*4+i,:) = [corner_in(1,j), corner_out(1,i), D_hope, Q/(4*4*2)];
    end
end


% SAVプロバイダー情報(SAV番号，拠点)
%SAV_grid = zeros(4*4+1,2);
SAV_grid = [
    1 1
    2 3
    3 5
    4 7
    5 9
    6 10
    7 14
    8 18
    9 22
    10 50
    11 54
    12 58
    13 62
    14 66
    15 70
    16 74
    17 78
    ];
%}



%% 放射型ネットワーク
% ネットワーク情報(リンク番号，発ノード，着ノード，自由旅行時間，容量)
D_hope = 12; %希望到着時間12.0889

rad_num = 3; %1辺に存在するリンク数(グリッドネットワークの大きさ)
time_b = sqrt((8*sqrt(3))/9);
free_travel_time = 1; %自由旅行時間
%capacity_rad = 25; %容量24.8161(単位時間あたりの流出可能？)
capacity_rad = 37; %容量37.2242(単位時間あたりの流出可能？)
%capacity_rad = 20; %(単位時間あたりの流出可能)
%capacity_rad = Inf; %無限大
road_rad = zeros((rad_num-1)*rad_num*4,5);
road_rad(:,4) = free_travel_time;
road_rad(:,5) = capacity_rad;
rad_head = [1,2,8,20,38,62,92,128,170];
num=0;
for i=1:rad_num
    % 連結リンク
    if i==1
        for j=1:6
            num=num+1;
            road_rad(num,1:3) = [num,1,1+j]; %内側→外側
            num=num+1;
            road_rad(num,1:3) = [num,1+j,1]; %外側→内側
        end
    else
        inner_list = zeros(1,(i*2-1)*6);
        outer_list = zeros(1,(i*2-1)*6);
        a=0;
        for j=1:(i-1)*6
            if j==1
                inner_list(1,(i*2-1)*6) = rad_head(i)+j-1;
                inner_list(1,1:2) = rad_head(i)+j-1;
                a=a+2;
            elseif rem((j-1),i-1)==0 && j>1
                inner_list(1,a+1:a+3) = rad_head(i)+j-1;
                a=a+3;
            else
                inner_list(1,a+1:a+2) = rad_head(i)+j-1;
                a=a+2;
            end
        end
        a=0;
        for j=1:i*6
%             if j==1
%                 outer_list(1,(i*2-1)*6) = radical_head(i+1)+j-1;
%                 outer_list(1,1) = radical_head(i+1)+j-1;
%                 a=a+1;
            if mod(j,i)==1
                outer_list(1,a+1) = rad_head(i+1)+j-1;
                a=a+1;
            else
                outer_list(1,a+1:a+2) = rad_head(i+1)+j-1;
                a=a+2;
            end
        end
        for j=1:(i*2-1)*6
            num=num+1;
            road_rad(num,1:3) = [num,inner_list(j),outer_list(j)]; %内側→外側
            num=num+1;
            road_rad(num,1:3) = [num,outer_list(j),inner_list(j)]; %外側→内側
        end
    end
        
    % 外縁リンク
    for j=1:i*6-1
        num=num+1;
        road_rad(num,1:3) = [num,rad_head(i+1)+j-1,rad_head(i+1)+j]; %時計回り
        num=num+1;
        road_rad(num,1:3) = [num,rad_head(i+1)+j,rad_head(i+1)+j-1]; %反時計回り
    end
    num=num+1; 
    road_rad(num,1:3) = [num,rad_head(i+1)+i*6-1,rad_head(i+1)]; %時計回り
    num=num+1; 
    road_rad(num,1:3) = [num,rad_head(i+1),rad_head(i+1)+i*6-1]; %反時計回り
end
road_rad(:,4) = free_travel_time;
road_rad(:,5) = capacity_rad;

for i=1:6
    %road_rad(265+(i-1)*rad_num*2:266+(i-1)*rad_num*2,5) = capacity*enlarge;
    %road_rad(271+(i-1)*rad_num*2:272+(i-1)*rad_num*2,5) = capacity*enlarge;
    road_rad(145+(i-1)*rad_num*2:146+(i-1)*rad_num*2,5) = capacity_rad*enlarge;
    road_rad(149+(i-1)*rad_num*2:150+(i-1)*rad_num*2,5) = capacity_rad*enlarge;
end


% OD(利用者)情報(出発ノード，到着ノード，希望到着時間，利用者数)
user_rad_one = zeros(6*6,4); %一方向OD表
user_rad_bi = zeros(6*6*2,4); %双方向OD表
corner_rad_out = [rad_head(rad_num+1), rad_head(rad_num+1)+rad_num, rad_head(rad_num+1)+rad_num*2, rad_head(rad_num+1)+rad_num*3, rad_head(rad_num+1)+rad_num*4, rad_head(rad_num+1)+rad_num*5];
corner_rad_in = [2,3,4,5,6,7];
for i=1:6
    for j=1:6
        user_rad_one((i-1)*6+j,:) = [corner_rad_out(1,i), corner_rad_in(1,j), D_hope, Q/(6*6)];
    end
end
for i=1:6 %外側
    for j=1:6 %内側
        user_rad_bi((i-1)*6+j,:) = [corner_rad_out(1,i), corner_rad_in(1,j), D_hope, Q/(6*6*2)];
        user_rad_bi((5+j)*6+i,:) = [corner_rad_in(1,j), corner_rad_out(1,i), D_hope, Q/(6*6*2)];
    end
end

% SAVプロバイダー情報(SAV番号，拠点)
if which_SAVarrange==0
    %SAV_grid = zeros(4*4+1,2);
    SAV_rad = [
        1 1
        2 2
        3 3
        4 4
        5 5
        6 6
        7 7
        8 21
        9 22
        10 24
        11 25
        12 27
        13 28
        14 30
        15 31
        16 33
        17 34
        18 36
        19 37
        20 rad_head(rad_num+1)
        21 rad_head(rad_num+1)+rad_num
        22 rad_head(rad_num+1)+rad_num*2
        23 rad_head(rad_num+1)+rad_num*3
        24 rad_head(rad_num+1)+rad_num*4
        25 rad_head(rad_num+1)+rad_num*5
        ];
    if rad_num<4
        SAV_rad = zeros(rad_head(rad_num+2)-1,2);
        for i=1:rad_head(rad_num+2)-1
            SAV_rad(i,:) = [i, i];
        end
    end
elseif which_SAVarrange==1
    SAV_rad = [
        1 rad_head(rad_num+1)
        2 rad_head(rad_num+1)+rad_num
        3 rad_head(rad_num+1)+rad_num*2
        4 rad_head(rad_num+1)+rad_num*3
        5 rad_head(rad_num+1)+rad_num*4
        6 rad_head(rad_num+1)+rad_num*5
        ];
end



% % ネットワーク情報(リンク番号，発ノード，着ノード，自由旅行時間，容量)
% grid_num = 9; %グリッドネットワークの大きさ
% free_travel_time = 1; %自由旅行時間
% capacity = 50; %容量
% road_grid = zeros((grid_num-1)*grid_num*4,5);
% road_grid(:,4) = free_travel_time;
% road_grid(:,5) = capacity;
% for i=1:grid_num-1
%     for j=1:grid_num-1
%         road_grid((grid_num*4-2)*(i-1)+j,1:3) = [(grid_num*4-2)*(i-1)+j,grid_num*(i-1)+j,grid_num*(i-1)+j+1]; %右リンク
%         road_grid((grid_num*4-2)*(i-1)+grid_num-1+j,1:3) = [(grid_num*4-2)*(i-1)+grid_num-1+j,grid_num*(i-1)+j+1,grid_num*(i-1)+j]; %左リンク
%     end
%     for j=1:grid_num
%         road_grid((grid_num*4-2)*(i-1)+grid_num*2-2+j,1:3) = [(grid_num*4-2)*(i-1)+grid_num*2-2+j,grid_num*(i-1)+j,grid_num*i+j]; %下リンク
%         road_grid((grid_num*4-2)*(i-1)+grid_num*3-2+j,1:3) = [(grid_num*4-2)*(i-1)+grid_num*3-2+j,grid_num*i+j,grid_num*(i-1)+j]; %上リンク
%     end
% end
% for j=1:grid_num-1 %横リンク
%     road_grid((grid_num*4-2)*(grid_num-1)+j,1:3) = [(grid_num*4-2)*(grid_num-1)+j,grid_num*(grid_num-1)+j,grid_num*(grid_num-1)+j+1];
%     road_grid((grid_num*4-2)*(grid_num-1)+grid_num-1+j,1:3) = [(grid_num*4-2)*(grid_num-1)+grid_num-1+j,grid_num*(grid_num-1)+j+1,grid_num*(grid_num-1)+j];
% end