close all;
hold on;
w = 0.1;

%%(1) Lateral Inhibition
I1 = [1 1 1 1 1 0 0 0 0 0];
A1 = LateralInhib(I1,w)

% using a threshold of 0.85 enlight the edge
T = Threshold(A,0.85)

%% (2) Optical Illusion
I2 = [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];

% with lateral inhibition
A2 = LateralInhib(I2,w)

% with convolution
conv(I2,[-w 1 -w],'same')

%% (3) Mona Lisa
Mona = im2double(imread('MonaLisa.jpg')); %read the image
colormap('gray'); %set the colormap to gray

% test the filter for various filter sizes
for n= 1:12
    F_n = Filter(w,n);
    mona_f = conv2(Mona,F_n,'same');
    mona_bin = imbinarize(mona_f); % automatic threshold function 
    subplot(3,4,n);
    imagesc(mona_bin)
end

figure;
subplot(1,3,1);
colormap('gray')
imagesc(Mona)
subplot(1,3,2);
imagesc(imbinarize(conv2(Mona,Filter(w,10),'same'))) % Show the final Mona with only the edges

% thresholding of the original image
subplot(1,3,3);
imagesc(imbinarize(Mona))

%% (4) Hermann Grid
figure;
colormap('gray')
H = im2double(imread('hermann.jpg'));
for n= 1:6
    F_n = Filter(w,n);
    H_f = conv2(H,F_n,'same');
    subplot(2,3,n);
    imagesc(H_f)
end


%% FUNCTIONS

% for Ex 1 and 2, return the Activation vector for an input I and a
% constant w
function [A] = LateralInhib(I,w)
    N =length(I); A = zeros(1,N); % Initialization
    A(1) = A(1) - w*I(2); A(N) = A(N) - w*I(N-1); % Cells on the edges
    A(2:N-1) = A(2:N-1) - w*(I(1:N-2)+I(3:N)); % Inner cells
end

% Return a binary vector for an Activation vector A and a threshold tr
function [T] = Threshold(A,tr)
    T = A*0;
    T(find(A>=tr))=1;
end

% Return a filter F_n for a constant w and a size n
function [F_n] = Filter(w,n)
    w_n = ones(n,n)*w;
    F_n = [-w_n -w_n -w_n;
        -w_n ones(n,n) -w_n;
        -w_n -w_n -w_n];
end
