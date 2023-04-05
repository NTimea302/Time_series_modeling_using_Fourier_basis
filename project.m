clear variables
close all
clc
load('product_16.mat');
P = 12; % folosit mai mult pt formula la matrice, puteam sa in locuiesc direct cu 12, dar arata mai frumos asa :)
%% Split the data into identification and validation
t_id = time(1:88);    % timpul pt primele 80% 
t_val = time(89:109); % timpul pt ultimele 20% 
y_id = y(1:88);       % datele pt primele 80% 
y_val = y(89:109);    % datele pt ultimele 20% 

%% 
for m=1:7
    figure;
    %Crearea matricii ID
    for k=1:length(t_id) % adica i = 1->88
        for j=1:2+2*m
            phi_id(k,1)=1;        
            phi_id(k,2)=k;   
            r = mod(j,2);
            if r==1
                phi_id(k,j)=cos(2*pi*k*(j-1)/2/P);
            else
                phi_id(k,j)=sin(2*pi*k*(j-2)/2/P);
            end
        end
    end

    theta=phi_id\y_id; %teta e aproximatorul, se calculeaza doar pt primele 80% din valori
    yhat_id=phi_id*theta; %aproximarea pt primele 80%

    %Crearea matricii val
    for k=1:length(t_val) % adica i = 1->21
        for j=1:2+2*m
            phi_val(k,1)=1;        
            phi_val(k,2)=t_val(k);   
            r = mod(j,2);
            if r==1
                phi_val(k,j)=cos(2*pi*t_val(k)*(j-1)/2/P);
            else
                phi_val(k,j)=sin(2*pi*t_val(k)*(j-2)/2/P);
            end
        end
    end
    
    %theta nu se mai calculeaza pt ultimele 20%, ci se foloseste cea de
    % dinainte pt a calcula aproximarea pe ultimele 20%
    yhat_val=phi_val*theta; % aproximarea pt ultimele 20%
    
    % MSE pt primele 80%
    s=0;
    for k=1:length(t_id)
        s=s+(y_id(k)-yhat_id(k)).^2;
    end
    MSE_id(m)=s/length(t_id);

    % MSE pt ultimele 20%
    s=0;
    for k=1:length(t_val)
        s=s+(y_val(k)-yhat_val(k)).^2;
    end
    MSE_val(m)=s/length(t_val);
    
    % Ploturile (albastru sunt datele initale, verde e aprox pt primele 80%,
    % rosu e aprox pt ultimele 20%)
    plot(time,y),title("m="+ m + "MSEval = "+ MSE_val+ "MSEid = "+ MSE_id);
    hold on;
    plot(t_val,yhat_val,'r');
    plot(t_id,yhat_id,'g');legend('Real values','Aprox values','Aprox ID');
end

figure,
m=1:7;
subplot(211);plot(m, MSE_id); title("MSE ID");
hold on
subplot(212);plot(m, MSE_val); title("MSE VAL");

%% Creating the matrices for m optim
figure,
m = 3;
close all
% load('product_16.mat');
% P = 12; %folosit mai mult pt formula la matrice, puteam sa in locuiesc direct cu 12, dar arata mai frumos asa :)
% 
% t_id = time(1:88);  %timpul pt primele 80% 
% t_val = time(89:109); %timpul pt ultimele 20% 
% 
% y_id = y(1:88); %datele pt primele 80% 
% y_val = y(89:109); %datele pt ultimele 20% 

%Crearea matricii ID
for k=1:length(t_id) % adica i = 1->88
    for j=1:2+2*m
        phi_id_optim(k,1)=1;        
        phi_id_optim(k,2)=k;   
        r = mod(j,2);
        if r==1
            phi_id_optim(k,j)=cos(2*pi*k*(j-1)/2/P);
        else
            phi_id_optim(k,j)=sin(2*pi*k*(j-2)/2/P);
        end
    end
end
theta=phi_id_optim\y_id; %teta e aproximatorul, se calculeaza doar pt primele 80% din valori
yhat_id=phi_id_optim*theta; %aproximarea pt primele 80%

%Crearea matricii ID
for k=1:length(t_val) % adica i = 1->21
    for j=1:2+2*m
        phi_val_optim(k,1)=1;        
        phi_val_optim(k,2)=t_val(k);   
        r = mod(j,2);
        if r==1
            phi_val_optim(k,j)=cos(2*pi*t_val(k)*(j-1)/2/P);
        else
            phi_val_optim(k,j)=sin(2*pi*t_val(k)*(j-2)/2/P);
        end
    end
end

%teta nu se mai calculeaza pt ultimele 20%, ci se foloseste cea de
%dinainte pt a calcula aproximarea pe ultimele 20%
yhat_val = phi_val_optim * theta; %aproximarea pt ultimele 20%

%Ploturile (albastru sunt datele initale, verde e aprox pt primele 80%,
%rosu e aprox pt ultimele 20%)
plot(time,y),title("m_o_p_t_i_m = 3");
hold on;
plot(t_val, yhat_val, 'r');
plot(t_id, yhat_id, 'g');legend('Real values', 'Aprox values', 'Aprox ID');

figure;
plot(time(89:109), y(89:109)),title("Validation data for m_o_p_t_i_m = 3 ");
hold on;
plot(t_val, yhat_val,'r');
