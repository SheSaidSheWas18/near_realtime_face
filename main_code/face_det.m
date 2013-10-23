% global_trflag stands for global track flag which is set if the input is a
% video so that tracking can be done. Otherwise, if it is just a bunch of
% random images, there is no need of tracking.

function face_det(global_trflag)

%Load images
imgs = dir('data_dir/*jpg');

%Load model
load dramanan_files/face_p146_small.mat;
model.interval = 5;
model.thresh = min(-1, model.thresh);
thresh = model.thresh;

% define the mapping from view-specific mixture id to viewpoint
if length(model.components)==13
    posemap = 90:-15:-90;
elseif length(model.components)==18
    posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
else
    error('Can not recognize this model');
end

%Make a strucuture for previous detections that can be used for tracking
track_prev.s  = 0;
track_prev.c  = 0;
track_prev.xy = 0;
track_prev.level = 0;
track_prev.coords = [-1 -1 -1 -1];


for i=1:numel(imgs)
    %Load image
    im_path = ['data_dir/' imgs(i).name];
    im = imread(im_path);
    
    %Get viola-jones bounding box
    [n det] = runFacedet(im_path,save_flag);
    
    % No detection
    if n == 0
        %empty the previous detection strucuture because no tracking.
        continue;
    end;
    
    %For all the viola-jones detections, check if there is any overlap
    %with the previous detection. If not there is no tracking.
    
    for ii = 1:n
        
        if global_trflag == 1 && i ~= 1
            
            %Checking for overlap between the current vj detection and tracked.
            [check track_num] = check_overlap([det(ii,1),det(ii,3),det(ii,2),det(ii,4)],track_prev);
            
            if check == -1
                %Pass empty set because no tracking
            else
                
                %Tracking is there.
            end;
            
            
        else
            %Pass empty set to Face_detect.
        end;
        
        
        % Face_detect part
        if check == -1
            ww = det(ii,2) - det(ii,1) + 1;
            hh = det(ii,4) - det(ii,3) + 1;
            face_bb = [1,size(im,2),1,size(im,1)];
            face_bb(1) = max(face_bb(1),(round(det(ii,1)-0.3*ww)));
            face_bb(2) = min(face_bb(2),(round(det(ii,2)+0.3*ww)));
            face_bb(3) = max(face_bb(3),(round(det(ii,3)-0.3*hh)));
            face_bb(4) = min(face_bb(4),(round(det(ii,4)+0.3*hh)));
            face_box = im(face_bb(3):face_bb(4),face_bb(1):face_bb(2),:);
        else
            %Get the face_box of the previous tracked frame from the track
            %structure.
            
        end;
        
        %Get pyramid level
        if check == -1
            lvl = max(7,ceil(((log(size(face_box,1)/(34*4))/(log(2)))*5)+6));
        else
            lvl = max(7,bs.level);
        end;
        
        %Compute feature pyramid
        pyrm = featpyramid(face_box,model,lvl);
        
        
        
        
        
        
    end;
end;
