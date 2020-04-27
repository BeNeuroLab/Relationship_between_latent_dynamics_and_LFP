function[] = CovMat(dataTensor,surrTensor)
    if length(size(dataTensor))~=3 || length(size(surrTensor))~=3
        quit
    end
    
    close all
    [dataT, dataN, dataC, ~] = extractFeatures(dataTensor);
    [surrT, surrN, surrC, ~] = extractFeatures(surrTensor);

    figure; subplot(1,2,1); imagesc(dataT); subplot(1,2,2); imagesc(surrT); suptitle('Time'); colormap 'jet';
    figure; subplot(1,2,1); imagesc(dataN); subplot(1,2,2); imagesc(surrN); suptitle('Neurons'); colormap 'jet';
    figure; subplot(1,2,1); imagesc(dataC); subplot(1,2,2); imagesc(surrC); suptitle('Trials'); colormap 'jet';
        

end