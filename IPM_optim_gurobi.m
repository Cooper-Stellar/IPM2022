%% スイッチ
which_network = 4; %4=グリッド，6=放射
which_OD = 2; %1=一方向，2=双方向
which_short = 1; %0=全部，1=ショートカット


%% 与件のパラメータ
alfa = 40; %移動時間価値(円/min)
alfa_e = 25; %早着時間価値(円/min)
alfa_l = 50; %遅着時間価値(円/min)
C = 4; % 車両容量(人/台)

if which_network==4
    road = road_grid;
    SAV = SAV_grid;
    if which_OD==1
        user = user_grid_one;
        if which_short==1
            node_flow_user = node_flow_user_grid_one;
        end
    elseif which_OD==2
        user = user_grid_bi;
        if which_short==1
            node_flow_user = node_flow_user_grid_bi;
        end
    end

    T = 25; %時間ステップ数25or30
    alfa = alfa*1;

    N = max([max(road(:,2)) max(road(:,3))]); %ノード数
    L = height(road); %リンク数
    t = repmat(road(:,4),1,T-1); %(L,(T-1))に対応した自由旅行時間
    limit = repmat(road(:,5),1,T-1); %(L,(T-1))に対応した道路容量

    OD = height(user); %ODトリップ数
    O_node = user(:,1)'; %各経路の出発ノード
    D_node = user(:,2)'; %各経路の到着ノード
    D_hope = user(:,3)'; %各経路の希望到着時間
    Q = user(:,4)'; %各経路の利用者数

    R = height(SAV); %SAVプロバイダーの数
    R_base = SAV(:,2)';% 拠点
    d = 20*t; %リンクサービス費用(円/time)
    beta = 50; %アクティブ車両数１台あたりのコスト係数(円/台・time)※デフォルト50
    ganma = 100; %SAV車両１台あたりの保有コスト係数(円/台)※デフォルト100

elseif which_network==6
    road = road_rad;
    SAV = SAV_rad;
    if which_OD==1
        user = user_rad_one;
        if which_short==1
            node_flow_user = node_flow_user_rad_one;
        end
    elseif which_OD==2
        user = user_rad_bi;
        if which_short==1
            node_flow_user = node_flow_user_rad_bi;
        end
    end

    time_b = sqrt((8*sqrt(3))/9);
    T = 20; %時間ステップ数25/time_b=20.1482
    alfa = alfa*time_b;
    alfa_e = alfa_e*time_b;
    alfa_l = alfa_l*time_b;

    N = max([max(road(:,2)) max(road(:,3))]); % ノード数
    L = height(road); % リンク数
    t = repmat(road(:,4),1,T-1); % (L,(T-1))に対応した自由旅行時間
    limit = repmat(road(:,5),1,T-1); % (L,(T-1))に対応した道路容量

    OD = height(user); % ODトリップ数
    O_node = user(:,1)'; % 各経路の出発ノード
    D_node = user(:,2)'; % 各経路の到着ノード
    D_hope = user(:,3)'; % 各経路の希望到着時間
    Q = user(:,4)'; % 各経路の利用者数

    R = height(SAV); %SAVプロバイダーの数
    R_base = SAV(:,2)'; %拠点
    d = 20*t*time_b; %リンクサービス費用(円/time)
    beta = 50*time_b; %アクティブ車両数１台あたりのコスト係数(円/台・time)
    ganma = 100; %SAV車両１台あたりの保有コスト係数(円/台)
end


% スケジュール費用
w = zeros(1,T*OD);
for i=1:OD
    for j=1:D_hope(1,i)
        w(1,T*(i-1)+j) = (D_hope(1,i)-j)*alfa_e;
    end
    for j=D_hope(1,i)+1:T
        w(1,T*(i-1)+j) = (j-D_hope(1,i))*alfa_l;
    end
end


%% 変数の準備
% 目的変数の係数
t_ij = repmat(reshape(t,1,[]),1,OD);
x_co = ones(1,L*(T-1)*OD).*alfa.*t_ij;
f_co = ones(1,T*OD).*w;
d_ij = repmat(reshape(d,1,[]),1,R);
weight_lambda = 0;
y_co = ones(1,L*(T-1)*R).*d_ij + weight_lambda;
g_co = zeros(1,T*R);
h_co = ones(1,T*R).*beta;
S_co = ones(1,R).*ganma;
coefficient = [x_co f_co y_co g_co h_co S_co]; % 係数のlinprog
column_long = L*(T-1)*(OD+R)+T*(OD+R*2)+R;


% 各主体のノードのフロー保存則のフローブロックを作成
flow_block = zeros(N*T,L*(T-1));
for i=1:L %リンク
    j=1; %時間
    while N*(j-1+road(i,4))+road(i,3)<=N*T
        flow_block(N*(j-1)+road(i,2),L*(j-1)+i) = 1; % 流出
        flow_block(N*(j-1+road(i,4))+road(i,3),L*(j-1)+i) = -1; % 流入
        j=j+1;
    end
end


% 時空間ネットワークの関係ないリンクの容量を0にする(Aeq)
useless_link = zeros(1,L*(T-1));
for i=1:L
    if road(i,4)>1
        for j=2:road(i,4)
            useless_link(1,L*(T-j)+i) = 1;
        end
    end
end
useless_block = [repmat(useless_link,1,OD) zeros(1,T*OD) repmat(useless_link,1,R) zeros(1,T*R*2+R)];


% 利用者のノードのフロー保存則
if which_short==0
    node_flow_user = sparse(zeros((N-1)*T*OD,1));
    node_flow_user = repmat(node_flow_user,1,column_long);
    user_flow_block = sparse(zeros(N*T*OD,1));
    user_flow_block = repmat(user_flow_block,1,L*(T-1)*OD);
    for i=1:OD
        user_flow_block(N*T*(i-1)+1:N*T*i,L*(T-1)*(i-1)+1:L*(T-1)*i) = flow_block;
    end
    for i=1:OD
        %i
        num=0;
        for j=1:N*T
            if mod(j,N) ~= mod(O_node(1,i),N)
                num=num+1;
                % 下の行のコードめっちゃ時間かかってる
                node_flow_user((N-1)*T*(i-1)+num,1:L*(T-1)*OD) = user_flow_block(N*T*(i-1)+j,:);
            end
        end
        if O_node(1,i)<D_node(1,i)
            D_node(1,i)=D_node(1,i)-1;
        end
        for j=1:T
            node_flow_user((N-1)*T*(i-1)+(N-1)*(j-1)+D_node(1,i),L*(T-1)*OD+T*(i-1)+j)=1; % fのところ
        end
    end
end
if which_OD==1
    node_flow_user_grid_one = node_flow_user;
elseif which_OD==2
    node_flow_user_grid_bi = node_flow_user;
end
%%}



% 利用者のODフロー保存則
OD_flow_user = sparse(zeros(OD,column_long));
for i=1:OD
    OD_flow_user(i,L*(T-1)*OD+T*(i-1)+1:L*(T-1)*OD+T*(i-1)+T) = ones(1,T)*(-1);
end

Text = "ここまで利用者制約準備"


% SAV車両のノードのフロー保存則
node_flow_sav = sparse(zeros(N*T*R,1));
node_flow_sav = repmat(node_flow_sav,1,column_long);
for i=1:R
    node_flow_sav(N*T*(i-1)+1:N*T*i,L*(T-1)*OD+T*OD+L*(T-1)*(i-1)+1:L*(T-1)*OD+T*OD+L*(T-1)*i) = flow_block;
    for j=1:T
        %node_flow_sav(N*T*(i-1)+N*(j-1)+i,L*(T-1)*(OD+R)+T*OD+T*(i-1)+j)=1; % gのところ
        node_flow_sav(N*T*(i-1)+N*(j-1)+R_base(1,i),L*(T-1)*(OD+R)+T*OD+T*(i-1)+j)=1; % gのところ
    end
end
node_flow_sav = node_flow_sav*(-1);


% SAV車両のODフロー保存則(gとhであらわされるもの)(Aeq,beq)
OD_flow_sav = sparse(zeros(T*R,1));
OD_flow_sav = repmat(OD_flow_sav,1,column_long);
for i=1:T
    for j=0:R-1
        OD_flow_sav(T*j+i,L*(T-1)*(OD+R)+T*OD+T*j+1:L*(T-1)*(OD+R)+T*OD+T*j+i) = ones(1,i); % gのところ
        OD_flow_sav(T*j+i,L*(T-1)*(OD+R)+T*(OD+R)+T*j+i) = 1; % hのところ
    end
end


% SAV車両の上限(h<=Sであらわされるもの)(A,b)
SAV_ub = sparse(zeros(T*R,1));
SAV_ub = repmat(SAV_ub,1,column_long);
for i=1:T
    for j=1:R
        %SAV_ub(T*j+i,L*(T-1)*(OD+R)+T*OD+T*j+1:L*(T-1)*(OD+R)+T*OD+T*j+i) = ones(1,i); % gのところ
        SAV_ub(T*(j-1)+i,L*(T-1)*(OD+R)+T*(OD+R)+T*(j-1)+i) = 1; % hのところ
        SAV_ub(T*(j-1)+i,L*(T-1)*(OD+R)+T*(OD+R*2)+j) = -1; % Sのところ
    end
end


% 時刻TでSAV車両は0
sav_end = sparse(zeros(R,1));
sav_end = repmat(sav_end,1,column_long);
for i=R-1:-1:0
    sav_end((i-R)*(-1),column_long-R-T*i) = 1;
end


% 道路容量制約(A,b)
road_matrix = sparse(zeros(L*(T-1),1));
road_matrix = repmat(road_matrix,1,column_long);
% for i=1:T-1
%     for j=1:L
%         for k=0:R-1
%             road_matrix(L*(i-1)+j,L*(T-1)*OD+T*OD+L*(T-1)*k+L*(i-1)+j) = 1; % yのところ
%         end
%     end
% end
road_matrix(:,L*(T-1)*OD+T*OD+1:L*(T-1)*OD+T*OD+L*(T-1)*R) = repmat(sparse(eye(L*(T-1))),1,R);


% MSサービス供給制約(A,b)
service_matrix = sparse(zeros(L*(T-1),1));
service_matrix = repmat(service_matrix,1,column_long);
% for i=1:T-1
%     for j=1:L
%         for k=0:OD-1
%             service_matrix(L*(i-1)+j,L*(T-1)*k+L*(i-1)+j) = 1;% xのところ
%         end
%         for k=0:R-1
%             service_matrix(L*(i-1)+j,L*(T-1)*OD+T*OD+L*(T-1)*k+L*(i-1)+j) = C*(-1); %yのところ
%         end
%     end
% end
eye_sub = sparse(eye(L*(T-1)));
service_matrix(:,1:L*(T-1)*OD) = repmat(eye_sub,1,OD); %xのところ
service_matrix(:,L*(T-1)*OD+T*OD+1:L*(T-1)*OD+T*OD+L*(T-1)*R) = repmat(sparse(eye(L*(T-1))),1,R)*C*(-1); %yのところ

Text = "ここまで制約準備"



%% linprogを回す
A = sparse([road_matrix; service_matrix; SAV_ub]);
%A_S = sparse(A);
b = sparse([ones(1,L*(T-1)).*reshape(limit,1,[]) zeros(1,L*(T-1)+T*R)]);
Aeq = sparse([node_flow_user; OD_flow_user; node_flow_sav; OD_flow_sav; sav_end; useless_block]);
%Aeq_S = sparse(Aeq);
beq = sparse([zeros(1,(N-1)*T*OD) ones(1,OD).*Q.*(-1) zeros(1,N*T*R) zeros(1,T*R) zeros(1,R) 0]);
% lb_SAVnet = ones(1,T*R);
% for i=1:R
%     lb_SAVnet(1,T*(i-1)+1:T*i) = S(1,i);
% end
lb = [zeros(1,(L*(T-1)*(OD+R)+T*OD)) Inf(1,T*R).*(-1) zeros(1,T*R+R)];
% ub_ODflow = ones(1,T*OD);
% for i=1:OD
%     ub_ODflow(1,T*(i-1)+1:T*i) = Q(1,i);
% end
ub = [Inf(1,column_long)];
%ub = sparse([C.*repmat(reshape(limit,1,[]),1,OD) ub_ODflow repmat(reshape(limit,1,[]),1,R) lb_SAVnet lb_SAVnet]);

Text = "ここまで最適化問題準備"
%% 最適化問題の実行
%options = optimoptions('linprog','Algorithm','interior-point-legacy');
%options = optimoptions('linprog','Algorithm','interior-point');
%[EPX,EP,exitflag,output,lambda] = linprog(coefficient, A, b, Aeq, beq, lb, ub);
%[EPX,EP,exitflag,output,lambda] = linprog(coefficient, A_S, b, Aeq_S, beq, lb, ub);

% gurobiの用意
model.obj = coefficient;
model.A = [A; Aeq]; % A must be sparse
%model.modelsense = 'Max';
model.rhs = full([b(:); beq(:)]); % rhs must be dense
model.lb = lb;
model.ub = ub;
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];

% gurobiの実行
result = gurobi(model);

% gurobiの結果
EP = result.objval;
%disp(result.x);
lambda.lower = [];
lambda.upper = [];
lambda.ineqlin = [];
lambda.eqlin = [];
if isfield(result,'rc')
    lambda.lower = max(0,result.rc);
    lambda.upper = -min(0,result.rc);
end
if isfield(result,'pi')
    if ~isempty(A)
        lambda.ineqlin = -result.pi(1:size(A,1));
    end
    if ~isempty(Aeq)
        lambda.eqlin = -result.pi((size(A,1)+1):end);
    end
end
%disp(lambda.lower);
%disp(lambda.upper);
%disp(lambda.ineqlin);
%disp(lambda.eqlin);
%lambda.lower



%% 均衡状態の数値
EPX = result.x;
EP_x = reshape(EPX(1:L*(T-1)*OD,1),L,(T-1)*OD);
EP_x_3dim = zeros(L,T-1,OD);
EP_x_sum = zeros(L,T-1);
for i=1:OD
    EP_x_3dim(:,:,i) = EP_x(:,(T-1)*(i-1)+1:(T-1)*i);
    EP_x_sum = EP_x_sum + EP_x_3dim(:,:,i);
end
if OD==1
    EP_x;
elseif OD==2
    EP_x1 = EP_x(:,1:T-1);
    EP_x2 = EP_x(:,(T-1)+1:(T-1)*2);
elseif OD==3
    EP_x1 = EP_x(:,1:T-1);
    EP_x2 = EP_x(:,(T-1)+1:(T-1)*2);
    EP_x3 = EP_x(:,(T-1)*2+1:(T-1)*3);
elseif OD==4
    EP_x1 = EP_x(:,1:T-1);
    EP_x2 = EP_x(:,(T-1)+1:(T-1)*2);
    EP_x3 = EP_x(:,(T-1)*2+1:(T-1)*3);
    EP_x4 = EP_x(:,(T-1)*3+1:(T-1)*4);
end
EP_f = reshape(EPX(L*(T-1)*OD+1:L*(T-1)*OD+T*OD,1),T,OD)';
EP_y = reshape(EPX(L*(T-1)*OD+T*OD+1:L*(T-1)*(OD+R)+T*OD,1),L,(T-1)*R);
EP_y_3dim = zeros(L,T-1,R);
EP_y_sum = zeros(L,T-1);
for i=1:R
    EP_y_3dim(:,:,i) = EP_y(:,(T-1)*(i-1)+1:(T-1)*i);
    EP_y_sum = EP_y_sum + EP_y_3dim(:,:,i);
end
%EP_y_exist = (EP_y_sum > 0.001);
if R==1
    EP_y1 = EP_y(:,1:T-1);
elseif R==2
    EP_y1 = EP_y(:,1:T-1);
    EP_y2 = EP_y(:,(T-1)+1:(T-1)*2);
elseif R==3
    EP_y1 = EP_y(:,1:T-1);
    EP_y2 = EP_y(:,(T-1)+1:(T-1)*2);
    EP_y3 = EP_y(:,(T-1)*2+1:(T-1)*3);
elseif R==4
    EP_y1 = EP_y(:,1:T-1);
    EP_y2 = EP_y(:,(T-1)+1:(T-1)*2);
    EP_y3 = EP_y(:,(T-1)*2+1:(T-1)*3);
    EP_y4 = EP_y(:,(T-1)*3+1:(T-1)*4);
end
EP_g = reshape(EPX(L*(T-1)*(OD+R)+T*OD+1:L*(T-1)*(OD+R)+T*(OD+R),1),T,R)';
EP_h = reshape(EPX(L*(T-1)*(OD+R)+T*(OD+R)+1:L*(T-1)*(OD+R)+T*(OD+R*2),1),T,R)';
EP_S = reshape(EPX(L*(T-1)*(OD+R)+T*(OD+R*2)+1:L*(T-1)*(OD+R)+T*(OD+R*2)+R,1),1,R)';

%EP_SPU = (x_co+(lambda.ineqlin(L*(T-1)+1:L*(T-1)*2,1))')*EPX(1:L*(T-1)*OD,1) + f_co*EPX(L*(T-1)*OD+1:L*(T-1)*OD+T*OD,1);
EPX_traf = x_co*EPX(1:L*(T-1)*OD,1);
EPX_sche = f_co*EPX(L*(T-1)*OD+1:L*(T-1)*OD+T*OD,1);
EPY_linkcost = y_co*EPX(L*(T-1)*OD+T*OD+1:L*(T-1)*(OD+R)+T*OD,1);
EPY_activecost = h_co*EPX(L*(T-1)*(OD+R)+T*(OD+R)+1:L*(T-1)*(OD+R)+T*(OD+R*2),1);
EPY_ownership = S_co*EPX(L*(T-1)*(OD+R)+T*(OD+R*2)+1:L*(T-1)*(OD+R)+T*(OD+R*2)+R,1);



% ラグランジュ乗数の生成
EP_e = reshape(lambda.ineqlin(1:L*(T-1),1),L,(T-1));
EP_p = reshape(lambda.ineqlin(L*(T-1)+1:L*(T-1)*2,1),L,(T-1));
%EP_p = EP_p.*EP_y_exist;
EP_p_round = round(EP_p);
EP_u = reshape(lambda.eqlin(1:(N-1)*T*OD,1),N-1,T*OD);
if OD==1
    EP_u1 = EP_u(:,1:T);
elseif OD==2
    EP_u1 = EP_u(:,1:T);
    EP_u2 = EP_u(:,T+1:T*2);
elseif OD==3
    EP_u1 = EP_u(:,1:T);
    EP_u2 = EP_u(:,T+1:T*2);
    EP_u3 = EP_u(:,T*2+1:T*3);
end
EP_rho = lambda.eqlin((N-1)*T*OD+1:(N-1)*T*OD+OD,1);
EP_v = reshape(lambda.eqlin((N-1)*T*OD+OD+1:(N-1)*T*OD+OD+N*T*R,1),N,T*R);
if R==1
    EP_v1 = EP_v(:,1:T);
elseif R==2
    EP_v1 = EP_v(:,1:T);
    EP_v2 = EP_v(:,T+1:T*2);
elseif R==3
    EP_v1 = EP_v(:,1:T);
    EP_v2 = EP_v(:,T+1:T*2);
    EP_v3 = EP_v(:,T*2+1:T*3);
end
EP_phi = reshape(lambda.eqlin((N-1)*T*OD+OD+N*T*R+1:(N-1)*T*OD+OD+N*T*R+T*R,1),T,R)';
EP_v_dif = (-1)*EP_phi;
EP_pi = reshape(lambda.upper((T-1)*L*(OD+R)+T*(OD+R)+1:(T-1)*L*(OD+R)+T*(OD+R*2),1),T,R)';

result_show = [EP EPX_traf EPX_sche EPY_linkcost EPY_activecost EPY_ownership];
EP_S_yoko = EP_S';
EP_S_sum = sum(EP_S)
result_show(2)/(sum(Q)*alfa)