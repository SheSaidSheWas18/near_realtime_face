function parts = part_shiftdt(numparts,parts,check,bs,scale,padx,pady)

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
