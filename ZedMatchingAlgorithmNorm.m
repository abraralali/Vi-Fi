len = max(max(cellfun('size',Ad,1)));
MissingAd = NaN(len,1);
AllWeights = cell(size(Ad,1),size(Ad,2));
for r=1:size(AllWeights,1)
    for c=1:size(AllWeights,2)
        AllWeights{r,c} = MissingAd;
        AllWeights{r,c}(1:size(Ad{r,c},1))=Ad{r,c};
    end
end
Ad = AllWeights;
TP = 0;
FP = 0;
TN = 0;
FN = 0;

TPs=0;
FPs=0;
Accuracy = [];
HeadStd = [];
mulTimeSegment = 1;
subTP = nan(1,size(IsHere,2)-1);
subTP(gndtrackId) = 0;
CorrelationResult = zeros(size(IsHere,1),size(IsHere,2));
correlated = zeros(1,size(IsHere,2)-1);
% MinDist = inf(size(Ad,1),size(Ad,2));
for time_seg = 1 : len
    assrow=[];
    distM = cellfun(@(v)v(time_seg),Ad);
    %     if time_seg == 1
    %         MinDist = distM;
    %     end
    %     for row=1:size(distM,1)
    %         for col=1:size(distM,2)
    %             if distM(row,col) < MinDist(row,col)
    %                 MinDist(row,col) = distM(row,col);
    %             end
    %         end
    %     end
    MinDist = distM;
    if time_seg*global_windSize<size(IsHere,1)
        PeopleInFrame = sum(IsHere(mulTimeSegment:time_seg*global_windSize,2:end))>0;
    else
        PeopleInFrame = sum(IsHere(mulTimeSegment:end,2:end))>0;
    end
    
    %         FrameGNDAssignment = GNDAssignemnt.*PeopleInFrame;
    
    
    %     mdist = MinDist(assrow,asscol);
    %     for yy=1:numel(assrow)
    %         simTrackID(yy,:) = abs(MinDist(yy,:)-mdist(yy))<=10;
    %     end
    %     for y=1:numel(assrow)
    %         A = fix(assrow(y)/size(MinDist,1));
    %         assrow(y) = assrow(y)-size(MinDist,1)*A;
    %     end
    %     if time_seg*global_windSize<size(IsHere,1)
    %         for isHereInd = mulTimeSegment : time_seg*global_windSize
    %                 FrameMinDist = MinDist.*IsHere(isHereInd,2:end);
    %                 FrameMinDist(FrameMinDist==0)=Inf;
    %                 [assignment,cost] = munkres(FrameMinDist);% problem: matches only one BB with a phone, while we could have more than BB correlated with the phone
    %                 [assrow,asscol] = find(assignment);
    %                 CorrelationResult(isHereInd,1) = IsHere(isHereInd,1);
    %                 CorrelationResult(isHereInd,asscol+1) = assrow'.*IsHere(isHereInd,asscol+1);
    %         end
    %     else
    %         for isHereInd = mulTimeSegment : size(IsHere,1)
    FrameMinDist = MinDist.*IsHere(time_seg,2:end);
    FrameMinDist = FrameMinDist./max(max(FrameMinDist));
    FrameMinDist(FrameMinDist==0|isnan(FrameMinDist))=inf;
    Frames_FrameMinDist{time_seg,2} = FrameMinDist;
    %             [FrameMinDistrows,FrameMinDistcols] = find(~isnan(FrameMinDist));
    %             FrameMinDistGraph = FrameMinDist(unique(FrameMinDistrows),unique(FrameMinDistcols));
    %             [~,original_col] = rmmissing(FrameMinDist,2);
    %             [~,original_row] = rmmissing(FrameMinDist,1);
    %             original_col = ~original_col;
    %             original_row = ~original_row;
    FrameMinDistGraph = FrameMinDist;
    dim = size(FrameMinDistGraph);
    [max_nodes,max_dim] = max(dim);
    if max_dim == 2
        num_of_additional_phones = max_nodes - size(FrameMinDistGraph,1);
        Graph_FrameMinDist = [FrameMinDistGraph;inf(num_of_additional_phones,size(FrameMinDistGraph,2))];
    else % max_dim is 1
        num_of_additional_subjects = max_nodes - size(FrameMinDistGraph,2);
        Graph_FrameMinDist = [FrameMinDistGraph,inf(size(FrameMinDistGraph,1),num_of_additional_subjects)];
    end
    if sum(~isinf(Graph_FrameMinDist))>0
        m0 = min(Graph_FrameMinDist(isfinite(Graph_FrameMinDist(:))));   % get minimum in A
        m1 = max(Graph_FrameMinDist(isfinite(Graph_FrameMinDist(:))));   % get maximum in A
        Graph_FrameMinDistnorm = (Graph_FrameMinDist -m0)/(m1-m0 ); 
        [assignment,cost] = munkres(Graph_FrameMinDistnorm);
    else
        [assignment,cost] = munkres(Graph_FrameMinDist);
    end
    [assrow,asscol] = find(assignment);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1 = cellfun(@(v)v(time_seg,1),Adindiv(:,asscol));
    x2 = cellfun(@(v)v(time_seg,2),Adindiv(:,asscol));
    x3 = cellfun(@(v)v(time_seg,3),Adindiv(:,asscol));
    
%     Frames_FrameMinDist{time_seg,2} = x1;
    Frames_FrameMinDist{time_seg,3} = std(FrameMinDist(~isinf(FrameMinDist)));
%     Frames_FrameMinDist{time_seg,4} = x3;
   
    %%
    
    %             if time_seg==1
    %                 % Draw graphs:
    %                 [Distrows,Distcols] = find(~isnan(distM));
    %                 FrameDistGraph = distM(unique(Distrows),unique(Distcols));
    %                 dim = size(FrameDistGraph);
    %                 [max_nodes,max_dim] = max(dim);
    %                 if max_dim == 2
    %                     num_of_additional_phones = max_nodes - size(FrameMinDistGraph,1);
    %                     Graph_FrameDistGraph = [FrameDistGraph;inf(num_of_additional_phones,size(FrameDistGraph,2))];
    %                 else % max_dim is 1
    %                     num_of_additional_subjects = max_nodes - size(FrameDistGraph,2);
    %                     Graph_FrameDistGraph = [FrameDistGraph,inf(size(FrameDistGraph,1),num_of_additional_subjects)];
    %                 end
    %                 [as,~] = munkres(Graph_FrameDistGraph);
    %                 matching_graph = Graph_FrameDistGraph.*as;
    %                 matching_graph(matching_graph==Inf)=0;
    %                 n = size(matching_graph,2);
    %                 m = size(matching_graph,1);
    %                 big_matching_graph = [zeros(m,m), matching_graph;
    %                          matching_graph', zeros(n,n)];
    %                 g = graph(big_matching_graph);
    %                 g.Nodes.Name = {'phone_5' 'phone_4' 'phone_3' 'phone_2' 'phone_1' 'TrackID_4' 'TrackID_2' 'None1' 'None2' 'None3'}';
    %                 figure
    %                 h = plot(g,'EdgeLabel',g.Edges.Weight);%,'NodeLabel',g.Nodes.Name);
    %                 h.LineStyle = '-';
    %                 m = size(matching_graph,1);
    %                 n = size(matching_graph,2);
    %                 h.XData(1:m) = 1;
    %                 h.XData((m+1):end) = 2;
    %                 h.YData(1:m) = linspace(0,1,m);
    %                 h.YData((m+1):end) = linspace(0,1,n);
    %             end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             MatchedFrameDist = MinDist.*assignment;
    %             MatchedFrameDist = max(MatchedFrameDist);
    %             MatchedFrameDistSTD = [MatchedFrameDistSTD;nanstd(MatchedFrameDist)];
    %             existingPeople = find(IsHere(time_seg,2:end));
    %             NotCorrelatedPeopel = ~ismember(asscol,existingPeople);
    CorrelationResult(time_seg,1) = IsHere(time_seg,1);
    CorrelationResult(time_seg,asscol+1) = assrow'.*IsHere(time_seg,asscol+1);
    %             if NotCorrelatedPeopel==0
    %                 CorrelationResult(time_seg,NotCorrelatedPeopel+1) = -1.*IsHere(time_seg,asscol+1);
    %             end
    %         end
    %     end
    
    
    
    %         [gndrow,gndcol] = find(FrameGNDAssignment);
    %         gnd = [gndrow,gndcol];
    %         ass = [assrow,asscol]
    %         DistAss{time_seg} = [MinDist,assignment];
    
    %         MatchingResult{time_seg} = isequal(FrameGNDAssignment,assignment);
    
    %         if isempty(ass)% if none of the person is correlated
    %             FN = 3;
    %             AllFN = [AllFN;FN];
    %         else
    %             k = ismember(ass,gnd, 'rows');
    %             CorrectlyAssigned = ass(k,:)
    %             phoneid{ass(k,2)}
    %             kFN = ~ismember(gnd,ass, 'rows');
    %             negk = ~ismember(ass,gnd,'rows');
    %
    %             TP = sum(k)
    %             TPs = TPs + TP;
    %
    %
    %             FP = sum(negk)% the person should be correlated but it was not
    %             FPs = FPs + FP;
    %             AllFP = AllFP+FP;
    %             Accuracy = [Accuracy;TPs/(TPs+FPs)]
    %             HeadStd = [HeadStd;sum(PeopleInFrame)]
    %
    %             asTN = ~(assignment);
    %             ass_neg_Col = sum(asTN);
    %             ass_neg_idx = find(ass_neg_Col==length(gndrow));
    %
    %             gnTN = ~(GNDAssignemnt);
    %             gn_neg_Col = sum(gnTN);
    %             gn_neg_idx = find(gn_neg_Col==length(gndrow));
    %             TN_people = ismember(ass_neg_idx,gn_neg_idx);
    %             FN_People = ~ismember(ass_neg_idx,gn_neg_idx);
    %
    %             TN = TN + sum(TN_people);
    %             AllTN = AllTN+TN;
    %
    %             FN = FN + sum(FN_People);
    %             AllFN = AllFN +FN;
    %
    %             matched_elements = sum(k(:));   % total number of equal elements, I will get it as p (i fixed it earlier)
    %             matching_precentage(time_seg) = matched_elements/(length(gndcol));
    mulTimeSegment = mulTimeSegment+global_winShift;
    
    %         end
    %         precesions(time_seg,:) = [TP/length(gndrow),FP/length(gndrow),TN/length(gndrow),FN/length(gndrow)];
    %
    %         Precesion = sum(AllTP)/(sum(AllTP)+sum(AllFP));
    %         Recall = sum(AllTP)/(sum(AllTP)+sum(AllFN));
    %
    %         F1_score = 2*((Precesion*Recall)/(Precesion+Recall));
    %
    %         time_seg_precesions(time_seg) = F1_score;
    %         All_precesions(time_seg,:) = [sum(AllTP)/(3*time_seg),sum(AllFP)/(3*time_seg),sum(AllFN)/(3*time_seg)];
    %
    %         prev_distM = distM;
    
    %     if figFlag1
    % Draw graphs:
    %             matching_graph = ValidMinDist.*assignment;
    %             n = size(matching_graph,2);
    %             m = size(matching_graph,1);
    %             big_matching_graph = [zeros(m,m), matching_graph;
    %                      matching_graph', zeros(n,n)];
    %             g = digraph(big_matching_graph);
    %             phonesHolders = {'Mohamed','Hongyu','Sid','Murtadha','Hansi'};
    %             subjects = {'Hansi','Mohamed','Hongyu','Sid','Murtadha'};
    % %             g.Nodes.Name = {subjects,phonesHolders};
    %             h = plot(g,'EdgeLabel',g.Edges.Weight);
    %             h.LineStyle = '-';
    %             m = size(MinDist,1);
    %             n = size(MinDist,2);
    %             h.XData(1:m) = 1;
    %             h.XData((m+1):end) = 2;
    %             h.YData(1:m) = linspace(0,1,m);
    %             h.YData((m+1):end) = linspace(0,1,n);
    %         end
    
end
TotalPeople = 0;
matchID = 0;
% if k==15
%     phonesHolders = {'Mohamed','Hongyu','Sid','Murtadha','Hansi'};
% else
AllphonesHolders = {'Hansi','Mohamed','Hongyu','Sid','Murtadha'};
phonesHolders = AllphonesHolders(1:num_of_phones);
% end
subjects = {'Hansi','Mohamed','Hongyu','Sid','Murtadha'};

for ai = 1:size(CorrelationResult,1)
    phoneTP = [];
    phoneFP = [];
    phoneGND = [];
    fTP = 0;
    fFP = 0;
    fFN = 0;
    fTN = 0;
    [~,trackID]  = find(IsHere(ai,2:end)>0);
    peopleInFrame = sum(IsHere(ai,2:end)>0);
    TotalPeople = TotalPeople+peopleInFrame;
    frno = CorrelationResult(ai,1);
    for ti=1:length(trackID)
        TrackIDSub = ZedImgPosition{1,trackID(ti)}.SubName(1);
        if (CorrelationResult(ai,trackID(ti)+1)==0) && (IsHere(ai,trackID(ti)+1)==1)
            if ~any(find(TrackIDSub==phonesHolders))
                ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = "None";
                TN = TN+1;
                fTN = fTN+1;
            else
                ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = "None";
                FN = FN+1;
                fFN = fFN+1;
            end
        else
            phone = phonesHolders{CorrelationResult(ai,trackID(ti)+1)};
            phone = convertCharsToStrings(phone);
            
            ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = phone;
            %             ZedImgPosition{1,trackID(ti)}.Distance(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = confidenceMatrix(ai,trackID(ti)+1);
            
            %               TracktorgndPosition{1,trackID(ti)}.AssignedPhone = repmat(phone,size(TracktorgndPosition{1,trackID(ti)},1),1);
            ModDecision{end+1,1} = k;
            ModDecision{end,2} = CorrelationResult(ai,1);
            if strcmp(TrackIDSub,phone)
                %                 subjTP(trackID(ti)) = subjTP(trackID(ti))+1;
                TP = TP+1;
                fTP = fTP+1;
                sTP(CorrelationResult(ai,trackID(ti)+1))=sTP(CorrelationResult(ai,trackID(ti)+1))+1;
                phoneTP = [phoneTP,CorrelationResult(ai,trackID(ti)+1)];
            else
                FP = FP+1;
                sFP(CorrelationResult(ai,trackID(ti)+1))=sFP(CorrelationResult(ai,trackID(ti)+1))+1;
                phoneFP = [phoneFP,CorrelationResult(ai,trackID(ti)+1)];
                phoneGND = [phoneGND,find(phonesHolders==TrackIDSub)];
            end
            ModDecision{end,3} = TrackIDSub;
            ModDecision{end,4} = phone;
        end
    end
    %             TN = TN + peopleInFrame-(FP+TP+FN)
    frameAcc(ai) = fTP/peopleInFrame;
    frameTP(ai) = fTP;
    detected_people(ai)=peopleInFrame;
    frame_number(ai) = CorrelationResult(ai,1);
    seq_num(ai) = k;
    Frames_FrameMinDist{ai,1}=frameAcc(ai);
    
    Frames_FrameMinDist{ai,8}=detected_people(ai);
    Frames_FrameMinDist{ai,9}=phoneTP;
    Frames_FrameMinDist{ai,10}=phoneFP;
    Frames_FrameMinDist{ai,11}=phoneGND;
end
total_frameAcc = [total_frameAcc;seq_num',frame_number',frameAcc',detected_people',frameTP'];
%%
overallFrameAccuracy = nansum(frameAcc(1,:))/sum(~isnan(frameAcc(1,:)));
%%
% Acc = matchID/sum(~isnan(subjTP));
%%

for bb=1:size(ZedImgPosition,2)
    if ~isempty(ZedImgPosition{1,bb}) %&& size(ZedImgPosition{1,bb},1)>=Kframes
        totalbb = totalbb+size(ZedImgPosition{1,bb},1);
        [matchFrame,~] = find(ZedImgPosition{1,bb}.SubName==ZedImgPosition{1,bb}.AssignedPhone);
        numberOfFramesToCorrectMatch = size(matchFrame,1);
        TrackIdLengthAcc = numberOfFramesToCorrectMatch/size(ZedImgPosition{1,bb},1);
        TrackIdLength = [TrackIdLength;size(ZedImgPosition{1,bb},1),TrackIdLengthAcc];
    end
end
%%
%ACC = (TP+TN) / (TP+TN+FP+FN);
ACC = TP/sum(detected_people)
PR = TP / (TP + FP);
REC = TP / (TP + FN);
NPV = TN / (TN + FP +FN + TP);
FS = 2 * (PR * REC) / (PR + REC);
%%
fname = subFolders(k).name+"ZedSubPosition";
cd("C:\Users\aalali\CameraIMUCorrelation/Trials/NewTrials/"+subFolders(k).name)
save(fname,'ZedImgPosition');
cd('C:\Users\aalali\CameraIMUCorrelation')