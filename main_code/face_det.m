% global_trflag stands for global track flag which is set if the input is a
% video so that tracking can be done. Otherwise, if it is just a bunch of
% random images, there is no need of tracking.

function face_det(global_trflag)
    
    %Load images
    imgs = dir('data_dir/*jpg');
    
    %Load model
    load dramanan_files/face_p146_small.mat;
    
    for i=1:numel(imgs)
        
    end;

end