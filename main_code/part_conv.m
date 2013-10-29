function [resp_out, parts_out, scale_out] =  part_conv(parts,rlvl,model,resp,pyrmd,filters,st,en,components,bs,check)

padx  = pyrmd.padx;
pady  = pyrmd.pady;

for part_num = 1:numparts
    % ???
    f = parts(part_num).filterid;
    level = rlvl-parts(part_num).scale*model.interval;
    %Check if for this level convolution has already been        %done or not
    if isempty(resp{level})
        scale = pyrmd.scale(level);
        % ?? Which boxes are these
        bbs = -1*ones(4,length(filters));
        for cmpn_in = st:en
            parts_in = components{cmpn_in};
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
                    x2 = size(pyrmd.feat{level},2)-2;
                    y2 = size(pyrmd.feat{level},1)-2;
                end;
                f_in  = parts_in(part_num_in).filterid;
                bbs(:,f_in) = [x1-1 y1-1 x2-1 y2-1];
            end;
        end;
        resp{level} = fconv(pyrmd.feat{level},filters,1, length(filters),bbs);
    end;
    
    parts(part_num).score = resp{level}{f};
    parts(part_num).level = level;
    resp_out = resp;
    parts_out = parts;
    scale_out = scale;
end;