clc;
close all;

% Input the image file name from a given folder
folderpath='Test_Images/';
filename='colorwheel';
filetype='.png';
X=uint8(imread(strcat(folderpath,filename,filetype)));

% Seperate image by RGB colors
X1= X(:,:,1);
X2= X(:,:,2);
X3= X(:,:,3);

% Define the decomposition wavelet type
wavelet = 'db1';

% Get levels to quantize image by (can be 1 to 20) from user
N = input("Enter desired quanitization level between 1 and 20 (lower level = more compression): ");

% Perform a level 1 decomposition of the image along the RGB channels
[cA1,cH1,cV1,cD1] = dwt2(X1,wavelet);
[cA2,cH2,cV2,cD2] = dwt2(X2,wavelet);
[cA3,cH3,cV3,cD3] = dwt2(X3,wavelet);
[height, width] = size(cA1);

% Find N threshhold levels (only need to do one set since all sets share
% range of values)
threshA1=multithresh(cA1,N);
threshH1=multithresh(cH1,N);
threshV1=multithresh(cV1,N);
threshD1=multithresh(cD1,N);

% Constructs a vector such that the min value in each quantization interval
% is assigned to the N levels of the output image
valuesMinA1= [min(cA1(:)) threshA1];
valuesMinH1= [min(cH1(:)) threshH1];
valuesMinV1= [min(cV1(:)) threshD1];
valuesMinD1= [min(cD1(:)) threshV1];

% Quantize the components of the decomposed image using N threshhold levels
[qcA1, indexA1] = imquantize(cA1,threshA1,valuesMinA1);
[qcH1, indexH1] = imquantize(cH1,threshH1,valuesMinH1);
[qcV1, indexV1] = imquantize(cV1,threshV1,valuesMinV1);
[qcD1, indexD1] = imquantize(cD1,threshD1,valuesMinD1);

[qcA2, indexA2] = imquantize(cA2,threshA1,valuesMinA1);
[qcH2, indexH2] = imquantize(cH2,threshH1,valuesMinH1);
[qcV2, indexV2] = imquantize(cV2,threshV1,valuesMinV1);
[qcD2, indexD2] = imquantize(cD2,threshD1,valuesMinD1);

[qcA3, indexA3] = imquantize(cA3,threshA1,valuesMinA1);
[qcH3, indexH3] = imquantize(cH3,threshH1,valuesMinH1);
[qcV3, indexV3] = imquantize(cV3,threshV1,valuesMinV1);
[qcD3, indexD3] = imquantize(cD3,threshD1,valuesMinD1);


%---------------- HUFFMAN ENCODING ----------------------------%   
%-------------------- RED CHANNEL -----------------------------%
% Encode each of the red channel components and turn them into vectors
A1Vec = qcA1(:);
H1Vec = qcH1(:);
V1Vec = qcV1(:);
D1Vec = qcD1(:);

% Find the probability of each of the symbols
% Find unique symbols and how many of them occur in the image
[probA1, symbolsA1] = hist(A1Vec,double(unique(A1Vec)));
[probH1, symbolsH1] = hist(H1Vec,double(unique(H1Vec)));
[probV1, symbolsV1] = hist(V1Vec,double(unique(V1Vec)));
[probD1, symbolsD1] = hist(D1Vec,double(unique(D1Vec)));

% Determine the probability of each symbol occuring
probA1=probA1/sum(probA1);
probH1=probH1/sum(probH1);
probV1=probV1/sum(probV1);
probD1=probD1/sum(probD1);

% Create a Huffman code dictionary using the symbols and probabilities we
% just found
[RedA1HuffDict,avglenA1]=huffmandict(symbolsA1,probA1);
[RedH1HuffDict,avglenH1]=huffmandict(symbolsH1,probH1);
[RedV1HuffDict,avglenV1]=huffmandict(symbolsV1,probV1);
[RedD1HuffDict,avglenD1]=huffmandict(symbolsD1,probD1);

% Encode the red components using the indivdual dictionaries for them
encodedA1 = huffmanenco(A1Vec, RedA1HuffDict);
encodedH1 = huffmanenco(H1Vec, RedH1HuffDict);
encodedV1 = huffmanenco(V1Vec, RedV1HuffDict);
encodedD1 = huffmanenco(D1Vec, RedD1HuffDict);
%-------------------- RED CHANNEL -----------------------------%

%-------------------- GREEN CHANNEL -----------------------------%
% Encode each of the green channel components and turn them into vectors
A2Vec = qcA2(:);
H2Vec = qcH2(:);
V2Vec = qcV2(:);
D2Vec = qcD2(:);

% Find the probability of each of the symbols
% Find unique symbols and how many of them occur in the image
[probA2, symbolsA2] = hist(A2Vec,double(unique(A2Vec)));
[probH2, symbolsH2] = hist(H2Vec,double(unique(H2Vec)));
[probV2, symbolsV2] = hist(V2Vec,double(unique(V2Vec)));
[probD2, symbolsD2] = hist(D2Vec,double(unique(D2Vec)));

% Determine the probability of each symbol occuring
probA2=probA2/sum(probA2);
probH2=probH2/sum(probH2);
probV2=probV2/sum(probV2);
probD2=probD2/sum(probD2);

% Create a Huffman code dictionary using the symbols and probabilities we
% just found
[GreenA2HuffDict,avglenA2]=huffmandict(symbolsA2,probA2);
[GreenH2HuffDict,avglenH2]=huffmandict(symbolsH2,probH2);
[GreenV2HuffDict,avglenV2]=huffmandict(symbolsV2,probV2);
[GreenD2HuffDict,avglenD2]=huffmandict(symbolsD2,probD2);

% Encode the green components using the indivdual dictionaries for them
encodedA2 = huffmanenco(A2Vec, GreenA2HuffDict);
encodedH2 = huffmanenco(H2Vec, GreenH2HuffDict);
encodedV2 = huffmanenco(V2Vec, GreenV2HuffDict);
encodedD2 = huffmanenco(D2Vec, GreenD2HuffDict);
%-------------------- GREEN CHANNEL -----------------------------%

%-------------------- BLUE CHANNEL -----------------------------%
% Encode each of the green channel components and turn them into vectors
A3Vec = qcA3(:);
H3Vec = qcH3(:);
V3Vec = qcV3(:);
D3Vec = qcD3(:);

% Find the probability of each of the symbols
% Find unique symbols and how many of them occur in the image
[probA3, symbolsA3] = hist(A3Vec,double(unique(A3Vec)));
[probH3, symbolsH3] = hist(H3Vec,double(unique(H3Vec)));
[probV3, symbolsV3] = hist(V3Vec,double(unique(V3Vec)));
[probD3, symbolsD3] = hist(D3Vec,double(unique(D3Vec)));

% Determine the probability of each symbol occuring
probA3=probA3/sum(probA3);
probH3=probH3/sum(probH3);
probV3=probV3/sum(probV3);
probD3=probD3/sum(probD3);

% Create a Huffman code dictionary using the symbols and probabilities we
% just found
[BlueA3HuffDict,avglenA3]=huffmandict(symbolsA3,probA3);
[BlueH3HuffDict,avglenH3]=huffmandict(symbolsH3,probH3);
[BlueV3HuffDict,avglenV3]=huffmandict(symbolsV3,probV3);
[BlueD3HuffDict,avglenD3]=huffmandict(symbolsD3,probD3);

% Encode the green components using the indivdual dictionaries for them
encodedA3 = huffmanenco(A3Vec, BlueA3HuffDict);
encodedH3 = huffmanenco(H3Vec, BlueH3HuffDict);
encodedV3 = huffmanenco(V3Vec, BlueV3HuffDict);
encodedD3 = huffmanenco(D3Vec, BlueD3HuffDict);
%-------------------- BLUE CHANNEL -----------------------------%
%---------------- HUFFMAN ENCODING ----------------------------%

% Combine the 12 encoded components by adding their vectors together end to
% send with a unique seperating character in between each component. The
% vectors will be added together in the same order they were encoded (A1,
% H1, V1, D1 and then repeat for green then blue)
encodedStream = [encodedA1; 2; encodedH1; 2; encodedV1; 2; encodedD1; 2; encodedA2; 2; encodedH2; 2; encodedV2; 2; encodedD2; 2; encodedA3; 2; encodedH3; 2; encodedV3; 2; encodedD3];

% Write the encoded stream to a file
binText='binary.txt';
fileID = fopen(binText, 'w');
fprintf(fileID, '%d\n', encodedStream);
fclose(fileID);


% This stuff will be in the decoding script
% Import the data from the file
fileData=importdata(binText);

%----------------- EXTRACT RED CHANNEL COMPONENTS ----------------------%
% Get A1 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dA1 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get H1 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dH1 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get V1 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dV1 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get D1 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dD1 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end
%----------------- EXTRACT RED CHANNEL COMPONENTS ----------------------%

%----------------- EXTRACT GREEN CHANNEL COMPONENTS ----------------------%
% Get A2 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dA2 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get H2 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dH2 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get V2 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dV2 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get D2 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dD2 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end
%----------------- EXTRACT GREEN CHANNEL COMPONENTS ----------------------%

%----------------- EXTRACT BLUE CHANNEL COMPONENTS ----------------------%
% Get A3 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dA3 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get H3 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dH3 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get V3 from encoded data stream
for i = 1:length(fileData)
   if fileData(i) == 2
       dV3 = fileData(1:i-1);
       fileData = fileData(i+1:end);
       break;
   end    
end

% Get D3 from encoded data stream
dD3 = fileData;
%----------------- EXTRACT BLUE CHANNEL COMPONENTS ----------------------%

%----------------- DECODE RED COMPONENTS --------------------------------%
% Decode each red channel component
decodedA1 = huffmandeco(dA1, RedA1HuffDict);
decodedH1 = huffmandeco(dH1, RedH1HuffDict);
decodedV1 = huffmandeco(dV1, RedV1HuffDict);
decodedD1 = huffmandeco(dD1, RedD1HuffDict);

% Turn the red components from vectors back into arrays
decodedA1 = reshape(decodedA1, [height, width]);
decodedH1 = reshape(decodedH1, [height, width]);
decodedV1 = reshape(decodedV1, [height, width]);
decodedD1 = reshape(decodedD1, [height, width]);
%----------------- DECODE RED COMPONENTS --------------------------------%

%----------------- DECODE GREEN COMPONENTS -------------------------------%
% Decode each green channel component
decodedA2 = huffmandeco(dA2, GreenA2HuffDict);
decodedH2 = huffmandeco(dH2, GreenH2HuffDict);
decodedV2 = huffmandeco(dV2, GreenV2HuffDict);
decodedD2 = huffmandeco(dD2, GreenD2HuffDict);

% Turn the green components from vectors back into arrays
decodedA2 = reshape(decodedA2, [height, width]);
decodedH2 = reshape(decodedH2, [height, width]);
decodedV2 = reshape(decodedV2, [height, width]);
decodedD2 = reshape(decodedD2, [height, width]);
%----------------- DECODE GREEN COMPONENTS -------------------------------%

%----------------- DECODE BLUE COMPONENTS --------------------------------%
% Decode each blue channel component
decodedA3 = huffmandeco(dA3, BlueA3HuffDict);
decodedH3 = huffmandeco(dH3, BlueH3HuffDict);
decodedV3 = huffmandeco(dV3, BlueV3HuffDict);
decodedD3 = huffmandeco(dD3, BlueD3HuffDict);

% Turn the blue components from vectors back into arrays
decodedA3 = reshape(decodedA3, [height, width]);
decodedH3 = reshape(decodedH3, [height, width]);
decodedV3 = reshape(decodedV3, [height, width]);
decodedD3 = reshape(decodedD3, [height, width]);
%----------------- DECODE BLUE COMPONENTS --------------------------------%

% Recompose image into color channels using inverse wavelet transform and
% approximation, horizantal, vertical, and diagonal components
decrX1 = uint8(idwt2(decodedA1,decodedH1,decodedV1,decodedD1,wavelet));
decrX2 = uint8(idwt2(decodedA2,decodedH2,decodedV2,decodedD2,wavelet));
decrX3 = uint8(idwt2(decodedA3,decodedH3,decodedV3,decodedD3,wavelet));

% Reconstruct image using the 3 color channels
decQ(:,:,1)=decrX1;
decQ(:,:,2)=decrX2;
decQ(:,:,3)=decrX3;

% Displaying the original image
figure; 
imshow(X);
title('Original image');

% Show reconstructed compressed image
figure;
fileOutput=strcat(folderpath,filename,'_decQ',filetype);
imwrite(decQ,fileOutput);
imshow(decQ);
title("Compressed image with encoding");

% Calculate peak noise to signal ratio to determine performance
peaksnr = psnr(decQ,X);
