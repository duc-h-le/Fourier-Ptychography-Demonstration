%% Retrieving phase and HR image from LR intensity images
close all
clc
clear
%% Define target amplitude and phase
ampl = double(imread('goldhill.png'))/255;
sz = 501;
ampl = ampl(1:sz,1:sz);
phase = rgb2hsv(imread('fruits.png'));
phase = phase(1:sz,1:sz,3);

%% Fourier transform
f = sqrt(ampl).*exp(i*(phase-0.5)*2*pi);
F = myfft(f);
figure
imshow(img2color(f))
figure
imshow(img2color(F/500))
imwrite(img2color(f), 'phase+ampl.png')
%% Pick regions in Fourier transform
bw = 90; % cutoff frequency (radius of the captured region in Fourier domain)

sx = -120:40:120; % List of shifts in Fourier domain
sy = -120:40:120;

data = zeros(length(sx),length(sy), sz, sz);
section0 = zeros(length(sx),length(sy),sz,sz);  
section = zeros(length(sx),length(sy),sz,sz);  

for i = 1:1:length(sx)
    for j = 1:1:length(sy)              
        section0(i,j,(sz+1)/2-sx(i)-bw:(sz+1)/2-sx(i)+bw,(sz+1)/2-sy(j)-bw:(sz+1)/2-sy(j)+bw) = mycircle(bw);
        %section(i,j,:,:)bw = bw+20;
        delta = 0;
        section(i,j,(sz+1)/2-sx(i)-bw-delta:(sz+1)/2-sx(i)+bw+delta,(sz+1)/2-sy(j)-bw-delta:(sz+1)/2-sy(j)+bw+delta) = mycircle(bw+delta);
        flow = myifft(F.*squeeze(section(i,j,:,:)));
        data(i,j,:,:) = abs(flow).^2;
        %imshow(abs(myfft(squeeze(data(i,j,:,:))))/500)
        %imshow(squeeze(data(i,j,:,:)))
        %imshow(squeeze(section(i,j,:,:)))
        imwrite(abs(flow).^2, ['Intensity_LR\', num2str(i), num2str(j),'.png'])
        imwrite(abs(F.*squeeze(section(i,j,:,:)))/500, ['Fourier_LR\', num2str(i), num2str(j),'.png'])
        imwrite(abs(myfft(squeeze(data(i,j,:,:))))/500, ['Fourier_LR_real\', num2str(i), num2str(j),'.png'])
    end
end

%% Gerchberg-Saxton loop

I = squeeze(data((length(sx)+1)/2,(length(sy)+1)/2,:,:));
Ph = ones(sz, sz)*pi/2;

img = sqrt(I).*exp(1i*Ph);
fou = fftshift(fft2(img));

% Iterations
figure
num = length(sx);
M = spiral(length(sx));
[~,idx] = sort(M(:));
[n_list,m_list] = ind2sub([num,num],idx);
figure

dist = 0;

for iter = 1:1:5
for m = 1:1:length(sx)*length(sy)
    i = n_list(length(sx)*length(sy)+1-m);
    j = m_list(length(sx)*length(sy)+1-m);
    F_LR = fou.*squeeze(section0(i,j,:,:)); % pick the section in freq domain
    ampl_LR = myifft(F_LR) ; % inverse fourier tranform
    ampl_LR = sqrt(squeeze(data(i,j,:,:))).*exp(1i*angle(ampl_LR)); % update the intensity with the captured LR images, keep phase
    F_LR = myfft(ampl_LR); % Fourier transform 
    
    fou = fou.*(1-squeeze(section0(i,j,:,:))); 
    fou = fou+F_LR.*squeeze(section0(i,j,:,:)); % Update the section in freq domain
    
    img = myifft(fou);

    imshow(img2color(fou)/500) % See how the frequency domain is updated
    %imshow(img2color(img)) % See how intensity+phase is updated   
    %imshow(phase2color(img)); % See how phase is updated   
    %imshow(abs(img).^2)    % See how intensity is updated   
end
end
 
%% Postprocessing
figure
imshow((abs(img).^2))
imwrite(abs(img).^2, 'HR_reconstructed.png')
ps = mod((angle(img)+pi)/2/pi,1);
figure
imshow(phase2color(img))
imwrite(phase2color(img), 'phase_reconstructed.png')

%% Functions
function F = myfft(f)
F = fftshift(fft2(f));
end
function F = myifft(f)
F = ifft2(ifftshift(f));
end

function y = img2color(f)
hsv = ones(size(f,1), size(f,2),3);
hsv(:,:,1) = 1-(angle(f)+pi)/(2*pi) ;
hsv(:,:,3) = abs(f);
y = hsv2rgb(hsv);
end

function y = phase2color(f)
hsv = ones(size(f,1), size(f,2),3);
hsv(:,:,1) = 1-(angle(f)+pi)/(2*pi) ;
y = hsv2rgb(hsv);
end

function y = angle2color(f)
hsv = ones(size(f,1), size(f,2),3);
hsv(:,:,1) = 1-(f+pi)/(2*pi) ;
y = hsv2rgb(hsv);
end

function circle = mycircle(bw)
circle = zeros(2*bw+1,2*bw+1);
for i = 1:1:2*bw+1
    circle(i,bw+1-round(sqrt(bw^2-(bw+1-i)^2)):bw+1+round(sqrt(bw^2-(bw+1-i)^2))) = 1;
end
end

function y = imgaussfiltcmplx(x, s)
y = imgaussfilt(real(x),s) + 1i*imgaussfilt(imag(x),s);
end

