close all;
clear all;
clc;
filename = '0118_7.jpg';
A  = imread(filename);
I = A;
[m,n,l] = size(A); %��óߴ�
W_cuts = 3;
H_cuts = 3;
BB = zeros(m, n);
AAG = zeros(W_cuts, H_cuts);

%% �����ɰ����
mask = zeros(m, n);
for p = 1 : m
    for q = 1 : n
        max_AAg = max([abs(A(p, q, 1) - A(p, q, 2)), abs(A(p, q, 3) - A(p, q, 2)), abs(A(p, q, 1) - A(p, q, 3))]);
        if max_AAg >= 30
            mask(p, q) = 1;
        end
    end
end
se1 = strel('disk',7);
mask = imdilate(mask, se1.Neighborhood);%ͼ��A1���ṹԪ��B����
mask = ~mask;% �����ɰ棬

%% ��ֵ��
Ag = im2double(rgb2gray(A));%�Ҷ�ͼ
BB = filtervalley(Ag); %ʹ��log�������ư�ɫ����
BB = BB(3:m+2, 3:n+2); %����Ϊԭʼ�ߴ�
BB = BB .* mask; %���ɰ�ɾ������
BB = edge(BB); %��ȡ��Ե��ת��Ϊ��ֵͼ
BW = BB;
se1 = strel('disk',3); %�����ȡ��Ե���пյ�����
BW = imdilate(BW, se1.Neighborhood);

BW(:, 1:20) = 0; %ɾ����߿�����
BW(:, end - 20: end) = 0;
BW(1:20, :) = 0;
BW(end - 20: end, :) = 0;

figure;
subplot(1,3,1);
imshow(BW * 100);
L = bwlabel(BW);
subplot(1,3,2);
imshow(L);

%% ɾ�����Ե�޹����򣬱���ֱ������
S = regionprops(L, 'Area');
W = regionprops(L, 'MajorAxisLength');
H = regionprops(L, 'MinorAxisLength');
BB = regionprops(L, 'BoundingBox');

S = [S.Area];
W = [W.MajorAxisLength];
H = [H.MinorAxisLength];
% W_H = W ./ H;
BB = [BB.BoundingBox];
W_H = zeros(length(S), 1);
WH = zeros(length(S), 1);
for i = 1 : length(W)
    W_H(i) = BB(:,4*i-1) / BB(:,4*i);
    WH(i) = S(i) / (BB(:,4*i-1) * BB(:,4*i));
end

threshold_div = 3;
L1 = ismember(L, find(W_H >= threshold_div)); % ֱ���жϱ�׼1������������������
L2 = ismember(L, find(W_H <= 1 / threshold_div)); % ֱ���жϱ�׼1������������������
L3 = ismember(L, intersect(union(find(W > 500), find(H > 500)), find(WH < 0.05))); % ֱ���жϱ�׼2��ռ������ϴ󣬵������ص�ϡ��
L = L1 + L2 + L3;
imshow(L  * 1000);

%% ����任��ȡ�߶�
[H,T,R] = hough(L);%�����ֵͼ��ı�׼����任��HΪ����任����I,RΪ�������任�ĽǶȺͰ뾶ֵ
P  = houghpeaks(H,30);%��ȡ3����ֵ��
lines=houghlines(BW,T,R,P);%��ȡ�߶�

%% �ϲ�����ͬһ��ֱ���ϵ��߶�
subplot(1,3,3);
imshow(I);
hold on;
xys = zeros(length(lines), 4); % �����߶ε���ֹ����
count = 1; % �����߶θ���
xy = [lines(1).point1; lines(1).point2];
xys(count, :) = reshape(xy, [4, 1]);
for k = 2:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    xy = reshape(xy, [4, 1]);
    Q1 = xy(1:2:3)';
    Q2 = xy(2:2:4)';
    flag = false;
    for l = 1 : count
        P1 = xys(l, 1:2:3);
        P2 = xys(l, 2:2:4);
        d1 = abs(det([Q2-Q1;P1-Q1]))/norm(Q2-Q1);% �����ĸ����㵽�߶εľ���
        d2 = abs(det([Q2-Q1;P2-Q1]))/norm(Q2-Q1);
        d3 = abs(det([P2-P1;Q1-P1]))/norm(P2-P1);
        d4 = abs(det([P2-P1;Q2-P1]))/norm(P2-P1);
        if max([d1, d2, d3, d4]) < 100
            P = [Q1; Q2; P1; P2];
            P = sortrows(P); % �ϲ��߶�
            xys(l, :) = [P(1,1), P(4,1), P(1,2), P(4,2)];
            flag = true;
            break;
        end
    end
    
    if flag == false
        count = count + 1;
        xys(count, :) = [Q1(1), Q2(1), Q1(2), Q2(2)];
    end
end

xys = xys(1:count, :);

%% �����ͼ
for k = 1:count
    xy = [xys(k, 1), xys(k, 3); xys(k, 2), xys(k, 4)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');%�����߶�
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');%���
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');%�յ�
end
hold off

%% ����б�ʣ�ɾ��б�ʹ���ĵ�

k = (xys(:, 3) - xys(:, 4)) ./ (xys(:, 1) - xys(:, 2));
v_h = length(find(abs(k) < 10)) / count;
if v_h >= 0.5
    xys = sortrows(xys, 4);
    for i = 1 : length(find(abs(k) < 10)) - 1
        up = min([xys(i,4), xys(i, 3)]);
        down = max([xys(i + 1,4), xys(i + 1, 3)]);
        sub_image =  I(up : down, :, : );
        imwrite(sub_image, [filename(1: end - 4), '_', num2str(i), '.jpg']);
    end
else
    xys = sortrows(xys, 2);
    for i = 1 : length(find(abs(k) >= 10)) - 1
        left = min([xys(i,1), xys(i, 2)]);
        right = max([xys(i + 1,1), xys(i + 1, 2)]);        
        sub_image =  I(:, left: right, : );
        imwrite(sub_image, [filename(1: end - 4), '_', num2str(i), '.jpg']);
    end
end


function int_image_range2 = filtervalley(int_image_range)
% use 5 * 5 Laplace of Gaussion(LoG) to calc edge
[m, n] = size(int_image_range);
L_o_G = [-2 -4 -4 -4 -2;
    -4 0 8 0 -4;
    -4 8 24 8 -4;
    -4 0 8 0 -4;
    -2 -4 -4 -4 -2];
int_image_range2 = conv2(int_image_range, L_o_G);
end
