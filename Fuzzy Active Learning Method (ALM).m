clear
close all
clc
x1 = linspace(-2,2,500);
figure,plot(x1,sinc(x1))
hold on, plot(x1,tanh(x1))
x2 = linspace(-3,4,500);
y_sinc = sinc(x1 + 3*x2);
y_saddle = tanh(0.3*x1 + 0.9*x2);
figure, plot(x1,y_sinc,x1,y_saddle);
figure, plot(x2,y_sinc,x2,y_saddle);
%% Take the training data
x1_train = x1;
x2_train = x2;
y_train = y_sinc;
%% Create a IDS Gridded plane
% takes the max and min of training data;

x1_min = min(x1_train);
x1_max = max(x1_train);

x2_min = min(x2_train);
x2_max = max(x2_train);
y_min = min(y_train);
y_max = max(y_train);

R = 1/1000; %Resolution
F1 = floor((x1_train - x1_min)/((x1_max - x1_min)*R));
F2 = floor((x2_train - x2_min)/((x2_max - x2_min)*R));
H = floor((y_train - y_min)/((y_max - y_min)*R));
%% 
x1_hat = x1_min:(x1_max - x1_min)/1000:x1_max;
x2_hat = x2_min:(x2_max - x2_min)/1000:x2_max;
y_hat = y_min:(y_max - y_min)/1000:y_max;
[X1,Y] = meshgrid(x1_hat,y_hat);
[X2,Y] = meshgrid(x2_hat,y_hat);
%% IDS plane
I0 = 1;
I1 = 1e-2;
CD1 = zeros(1001,1001);
for i = 1:length(F1)
    for j = -10:10
        for k = -10:10
            if((F1(i)+j<1)||(F1(i)+j>1000)||(H(i)+k>1000)||(H(i)+k<1))
                continue
            end
            CD1(F1(i)+j,H(i)+k) = CD1(F1(i)+j,H(i)+k) + I0*exp(-I1*j^2)*exp(-I1*k^2);
        end
    end
end    
CD2 = zeros(1001,1001);
for i = 1:length(F1)
    for j = -10:10
        for k = -10:10
            if((F2(i)+j<1)||(F2(i)+j>1000)||(H(i)+k>1000)||(H(i)+k<1))
                continue
            end
            CD2(F2(i)+j,H(i)+k) = CD2(F2(i)+j,H(i)+k) + I0*exp(-I1*j^2)*exp(-I1*k^2);
        end
    end
end    

%% 
close all
%CD_norm = CD1/max(CD1,[],'all');
figure,surf(Y,X1,CD1','edgecolor','none');
figure,surf(Y,X2,CD2','edgecolor','none');
%% 
% F_hat = linspace(x_min,x_max,1000);
%close all
% figure
% hold on
% for i=1:1000
%     for j=1:1000
%         if ~(CD1(i,j) == 0)
% %             temp=CD(i,j);
% %             disp(i);
% %             disp(j);
%             plot3(i,j,CD1(i,j));
%         end
%     end
% end
%% Narrow Path
clc
clear temp tempy tempol firstOccurence lastOccurence saay Spread a ind
for i = 1:size(CD1,1)
    tempy=CD1(i,:)>1.2;
    if isempty(find(tempy~=0,1,'first'))
        firstOccurence(i) = 1;
    else
    firstOccurence(i)=find(tempy~=0,1,'first');
    end
    if isempty(find(tempy~=0,1,'last'))
        lastOccurence(i) = size(CD1,2);
    else
        lastOccurence(i)=find(tempy~=0,1,'last');
    end
    [a ind]=max(CD1(i,:));
    y1_narrow(i) =  Y(ind);
end
Spread1 = lastOccurence - firstOccurence;
figure,p=plot(x1_hat,Spread1,'MarkerSize',20)
clear temp tempy tempol firstOccurence lastOccurence saay Spread
for i = 1:size(CD2,1)
    %temp = CD(i,:)>1
   % firstOccurence = find(temp~=0,1,'first');
    %lastOccurence = find(temp~=0,1,'last')
    tempy=CD2(i,:)>1.2;
    if isempty(find(tempy~=0,1,'first'))
        firstOccurence(i) = 1;
    else
    firstOccurence(i)=find(tempy~=0,1,'first');
    end
    if isempty(find(tempy~=0,1,'last'))
        lastOccurence(i) = size(CD2,2);
    else
        lastOccurence(i)=find(tempy~=0,1,'last');
    end
    [a ind]=max(CD2(i,:));
    y2_narrow(i) =  Y(ind);
end
Spread2 = lastOccurence - firstOccurence;
figure,p=plot(x2_hat,Spread2,'MarkerSize',20)
%%
%close all

% figure, plot(x1_hat,y1_narrow)
% figure, plot(x2_hat,y2_narrow)
% figure, plot(y2_narrow)
%% Plot Results
% close all
net = 1./Spread1 + 1./Spread2;
beta1 = (1./Spread1)./net;
beta2 = (1./Spread2)./net;
y = beta1.*y1_narrow + beta2.*y2_narrow;
figure,plot(0.3*x1_hat + 0.9*x2_hat,y)
hold on
plot(0.3*x1_hat + 0.9*x2_hat,tanh(0.3*x1_hat + 0.9*x2_hat))
figure,plot(x1_hat + 3*x2_hat,y)
hold on
plot(x1_hat + 3*x2_hat,sinc(x1_hat + 3*x2_hat))