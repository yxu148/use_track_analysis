% Stitch after knowing what tracks to stitch. Adapted from Isabel

% Track‐ID groups to merge into letter tracks
track_groups = struct( ...
    'A',[3, 7], ...
    'B',[2], ...
    'C',[9], ...
    'D',[15], ...
    'E',[11, 16], ...
    'F',[12], ...
    'G',[13, 4, 8], ...
    'H',[14] );


stitched_tracks = struct();                     % output container
fps = 20;                                       % frames s⁻¹


group_labels = fieldnames(track_groups);
for g = 1:numel(group_labels)
    group_name = group_labels{g};
    ids        = track_groups.(group_name);     % numeric IDs …


    % Pre‑allocate dq struct for this stitched track
    stitched_tracks.(group_name).dq = struct( ...
        'iloc', [], 'shead', [], 'smid', [], ...
        'speed', [], 'eti', [] );


    % ---- copy directly if only one source track ------------------------
    if numel(ids) == 1
        src  = eset.expt.track(ids).dq;
        stitched_tracks.(group_name).dq.iloc  = src.iloc;
        stitched_tracks.(group_name).dq.shead = src.shead;
        stitched_tracks.(group_name).dq.smid  = src.smid;
        stitched_tracks.(group_name).dq.speed = src.speed(:);
        stitched_tracks.(group_name).dq.eti   = src.eti(:);
        continue
    end


    % ---- merge consecutive tracks -------------------------------------
    for k = 1:numel(ids)-1
        t1 = eset.expt.track(ids(k));
        t2 = eset.expt.track(ids(k+1));


        % Fallback if a field is unexpectedly empty
        if isempty(t1.dq.iloc), t1.dq = t2.dq; end
        if isempty(t2.dq.iloc), t2.dq = t1.dq; end


        % Frames missing between the two source tracks
        framesGap = (t2.startFrame - t1.endFrame) - 1;


        % Extract final & initial samples -------------------------------
        endPos      = t1.dq.iloc(:,end);   startPos   = t2.dq.iloc(:,1);
        endHead     = t1.dq.shead(:,end);  startHead  = t2.dq.shead(:,1);
        endMid      = t1.dq.smid(:,end);   startMid   = t2.dq.smid(:,1);
        endTime     = t1.dq.eti(end);      startTime  = t2.dq.eti(1);
        endSpeed    = t1.dq.speed(end);    startSpeed = t2.dq.speed(1);


        % -------- bridge if there is a temporal gap --------------------
        if framesGap > 0
            % positions (iloc)
            bridgeX   = linspace(endPos(1),  startPos(1),  framesGap+2);
            bridgeY   = linspace(endPos(2),  startPos(2),  framesGap+2);
            bridgeLoc = [bridgeX(2:end-1);   bridgeY(2:end-1)];


            % head & mid
            bridgeHeadX = linspace(endHead(1), startHead(1), framesGap+2);
            bridgeHeadY = linspace(endHead(2), startHead(2), framesGap+2);
            bridgeHead  = [bridgeHeadX(2:end-1); bridgeHeadY(2:end-1)];


            bridgeMidX  = linspace(endMid(1),  startMid(1),  framesGap+2);
            bridgeMidY  = linspace(endMid(2),  startMid(2),  framesGap+2);
            bridgeMid   = [bridgeMidX(2:end-1); bridgeMidY(2:end-1)];


            % uniform speed across the gap
            gapDist  = norm(startPos - endPos);
            gapTime  = framesGap / fps;
            bridgeSp = repmat(gapDist / gapTime, framesGap, 1);


            % time stamps
            bridgeT  = linspace(endTime, startTime, framesGap+2);
            bridgeT  = bridgeT(2:end-1);
        else
            bridgeLoc = []; bridgeHead = []; bridgeMid = [];
            bridgeSp = [];  bridgeT    = [];
        end


        % -------- concatenate into stitched_track ----------------------
        if k == 1
            S = stitched_tracks.(group_name).dq;   % shorthand
            S.iloc  = t1.dq.iloc;
            S.shead = t1.dq.shead;
            S.smid  = t1.dq.smid;
            S.speed = t1.dq.speed(:);
            S.eti   = t1.dq.eti(:);
        else
            S = stitched_tracks.(group_name).dq;
        end


        S.iloc  = [S.iloc,  bridgeLoc,  t2.dq.iloc];
        S.shead = [S.shead, bridgeHead, t2.dq.shead];
        S.smid  = [S.smid,  bridgeMid,  t2.dq.smid];
        S.speed = [S.speed; bridgeSp;   t2.dq.speed(:)];
        S.eti   = [S.eti;   bridgeT(:); t2.dq.eti(:)];


        stitched_tracks.(group_name).dq = S;       % write back
    end
end