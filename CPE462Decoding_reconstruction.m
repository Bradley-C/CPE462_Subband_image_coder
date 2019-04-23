% define the wavelet type and how many level of reconstruction performed
wavelet = 'haar';
N = 1;
[LoD,HiD] = wfilters (wavelet,'d');

% Input ethropy encoded image


%Decode image

%Reconstruct image
% depending on how you decode into use idwt2 if you have the subbands
% seperated or waverec2 if they are in a single vectors with a bookkeeping
% matrix
rX1 = uint8(idwt2(cA1, cH1, cV1, cD1, wavelet));
rX2 = uint8(idwt2(cA2, cH2, cV2, cD2, wavelet)); 
rX3 = uint8(idwt2(cA3, cH3, cV3, cD3, wavelet)); 
rX1 = uint8(waverec2(CXC1,LXC1,LoD,HiD));
rX2 = uint8(waverec2(CXC2,LXC1,LoD,HiD));
rX3 = uint8(waverec2(CXC3,LXC1,LoD,HiD));

rX(:,:,1)=rX1;
rX(:,:,2)=rX2;
rX(:,:,3)=rX3;

figure;
imagesc (rX);
imwrite(rX,'out.jpg');

% To find the compression ratio
infoin = imfinfo('in.jpg');
infoout = imfinfo('out.jpg');
bin = infoin.FileSize;
bout = infoout.FileSize;
disp(bin);
disp(bout);
compression_ratio= bin/bout;
disp('Compression Ratio =');
disp(compression_ratio);

% Reconstruct the original image from the I level decomposition. 
% To reconstruct the original image from the wavelet decomposition structure
% Calculate mean square error and power signal to noise ratio
% X = double(X);
% X0 = double(X0);
% mse=(sum(sum(sum((X-X0).*(X-X0)))))/(r*c*p);
% PSNR=20*log10(255/sqrt(mse));
% disp('Mean square error = ');
% disp(mse);
% disp('PSNR =');
% disp(PSNR);
