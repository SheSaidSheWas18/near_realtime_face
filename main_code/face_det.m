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

BOXCACHESIZE = 100000;

[components,filters]  = modelcomponents(model);

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
        
        st = 1;
        en = length(components);
        
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
            bs = track_prev(track_num);
            part_xy = bs.coords;
            face_box = im(part_xy(2):part_xy(4),part_xy(1):part_xy(3),:);
            bs.xy(:,1) = bs.xy(:,1) - part_xy(1);
            bs.xy(:,2) = bs.xy(:,2) - part_xy(2);
            bs.xy(:,3) = bs.xy(:,3) - part_xy(1);
            bs.xy(:,4) = bs.xy(:,4) - part_xy(2);
            
            st = max(1,bs.c- 1);
            en = min(13,bs.c+1);
            
        end;
        
        BOXCACHESIZE = 100000;
        cnt = 0;
        boxes.s  = 0;
        boxes.c  = 0;
        boxes.xy = 0;
        boxes.level = 0;
        boxes.coords = [-1 -1 -1 -1];
        boxes(BOXCACHESIZE) = boxes;
        
        %Get pyramid level
        if check == -1
            lvl = max(7,ceil(((log(size(face_box,1)/(34*4))/(log(2)))*5)+6));
        else
            lvl = max(7,bs.level);
        end;
        
        %Compute feature pyramid
        pyrmd = featpyramid(face_box,model,lvl);
        resp    = cell(length(pyrmd.feat),1);
        
        padx  = pyrmd.padx;
        pady  = pyrmd.pady;
        
        %For components
        for cmpn_num = st:en
            %On pyramid levels
            for rlvl = lvl-1:lvl+1
                parts    = components{cmpn_num};
                numparts = length(parts);
                
                %Convolution for each part
                for part_num = 1:numparts
                    
                    % ???
                    f = parts(part_num).filterid;
                    level = rlvl-parts(part_num).scale*model.interval;
                    
                    %Check if for this level convolution has already been
                    %done or not
                    if isempty(resp{level})
                        scale = pyrmd.scale(level);
                        
                        % ?? Which boxes are these
                        bbs = -1*ones(4,length(filters));
                        
                        for cmpn_in = st:en
                            parts_in = components{cmpn_curr_in};
                            numparts_in = length(parts_in);
                            
                            for part_num_in = 1:numparts_in
                                % ?? If else
                                if check == 1 && size(bs.xy,1) == numparts_in
                                    x1 = floor(((bs.xy(k1,1) - 1)/scale) + 1 + padx);
                                    y1 = floor(((bs.xy(k1,2) - 1)/scale) + 1 + pady);
                                    x2 = ceil(((bs.xy(k1,3) - 1)/scale) + 1 + padx);
                                    y2 = ceil(((bs.xy(k1,4) - 1)/scale) + 1 + pady);
                                else
                                    x1 = 3;
                                    y1 = 3;
                                    x2 = size(pyra.feat{level},2)-2;
                                    y2 = size(pyra.feat{level},1)-2;
                                end;
                                f_in  = parts_in(part_num_in).filterid;
                                bbs(:,f_in) = [x1-1 y1-1 x2-1 y2-1];
                            end;
                        end;
                        resp{level} = fconv(pyrmd.feat{level},filters,1, length(filters),bbs);
                    end;
                    
                    parts(part_num).score = resp{level}{f};
                    parts(part_num).level = level;
                end;
                
                %Shift distance transform
                for k=numparts:-1:2
                    child = parts(k);
                    par   = child.parent;
                    [Ny,Nx,~] = size(parts(par).score);
                    
                    % why this condition??
                    if check == 1 && size(bs.xy,1) == numparts
                        ccx1 = floor(((bs.xy(k,1) - 1)/scale) + 1 );%+ padx);
                        ccy1 = floor(((bs.xy(k,2) - 1)/scale) + 1 );%+ pady);
                        ccx2 = floor(((bs.xy(k,3) - 1)/scale) + 1 + padx);
                        ccy2 = floor(((bs.xy(k,4) - 1)/scale) + 1 + pady);
                        
                        ppx1 = floor(((bs.xy(par,1) - 1)/scale) + 1 );%+ padx);
                        ppy1 = floor(((bs.xy(par,2) - 1)/scale) + 1 );%+ pady);
                        ppx2 = floor(((bs.xy(par,3) - 1)/scale) + 1 + padx);
                        ppy2 = floor(((bs.xy(par,4) - 1)/scale) + 1 + pady);
                    else
                        ccx1 = 1;
                        ccx2 = size(child.score,2);
                        ccy1 = 1;
                        ccy2 = size(child.score,1);
                        
                        ppx1 = 1;
                        ppx2 = Nx;
                        ppy1 = 1;
                        ppy2 = Ny;
                    end;
                    
                    [msg,parts(k).Ix,parts(k).Iy] = shiftdt2(child.score, child.w(1),child.w(2),child.w(3),child.w(4), ...
                        child.startx, child.starty, Nx, Ny, child.step, ...
                        ppx1, ppy1, ppx2, ppy2, ccx1, ccy1, ccx2, ccy2);
                    
                    tmpmsg = msg;
                    msg = -1000*ones(size(tmpmsg));
                    msg(ppy1:ppy2,ppx1:ppx2) = tmpmsg(ppy1:ppy2,ppx1:ppx2);
                    parts(par).score = parts(par).score + msg;
                end;
                
                % Add bias to root score
                rscore = parts(1).score + parts(1).w;
                [Y,X] = find(rscore >= thresh);
                if ~isempty(X)
                    XY = backtrack( X, Y, parts, pyra);
                end;
                
                
                % Walk back down tree following pointers
                for i = 1:length(X)
                    x = X(i);
                    y = Y(i);
                    
                    if cnt == BOXCACHESIZE
                        b0 = nms_face(boxes,0.3);
                        clear boxes;
                        boxes.s  = 0;
                        boxes.c  = 0;
                        boxes.xy = 0;
                        boxes.level = 0;
                        boxes(BOXCACHESIZE) = boxes;
                        cnt = length(b0);
                        boxes(1:cnt) = b0;
                    end
                    
                    cnt = cnt + 1;
                    boxes(cnt).c = c;
                    boxes(cnt).s = rscore(y,x);
                    boxes(cnt).level = rlevel;
                    boxes(cnt).xy = XY(:,:,i);
                end
            end;
        end;
        
        % change the boxes coordinates back to originals. TODO.
    end;
end;


