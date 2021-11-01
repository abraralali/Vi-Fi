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
% subTP = nan(1,size(IsHere,2)-1);
% subTP(gndtrackId) = 0;
CorrelationResult = zeros(size(IsHere,1),size(IsHere,2));
correlated = zeros(1,size(IsHere,2)-1);
distM = cellfun(@(v)v(1),Ad);
% MinDist = inf(size(Ad,1),size(Ad,2));
for f = 1 : size(IsHere,1)
    assrow=[];
    [~,currentTrackids] = find(IsHere(f,2:end));
    FrameMinDist = distM.*IsHere(f,2:end);
    FrameMinDist(FrameMinDist==0|isnan(FrameMinDist))=inf;
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
    [assignment,cost] = munkres(Graph_FrameMinDist);
    [assrow,asscol] = find(assignment);
%     Frames_FrameMinDist{f,1} = distM.*IsHere(f,2:end);
    %Frames_FrameMinDist{f,3} = nanstd(distM.*IsHere(f,2:end),[],2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             if time_seg==1
    %                 % Draw graphs:
    %                 [Distrows,Distcols] = find(~isnan(distM));
    %                 FrameDistGraph = distM(unique(Distrows),unique(Distcols));
    % %                 dim = size(FrameDistGraph);
    % %                 [max_nodes,max_dim] = max(dim);
    % %                 if max_dim == 2
    % %                     num_of_additional_phones = max_nodes - size(FrameMinDistGraph,1);
    % %                     Graph_FrameDistGraph = [FrameDistGraph;inf(num_of_additional_phones,size(FrameDistGraph,2))];
    % %                 else % max_dim is 1
    % %                     num_of_additional_subjects = max_nodes - size(FrameDistGraph,2);
    % %                     Graph_FrameDistGraph = [FrameDistGraph,inf(size(FrameDistGraph,1),num_of_additional_subjects)];
    % %                 end
    %                 [as,~] = munkres(FrameDistGraph);
    %                 matching_graph = FrameDistGraph.*as;
    %                 n = size(matching_graph,2);
    %                 m = size(matching_graph,1);
    %                 big_matching_graph = [zeros(m,m), matching_graph;
    %                          matching_graph', zeros(n,n)];
    %                 g = digraph(big_matching_graph);
    %                 %g.Nodes.Name = {'Pixel 2' 'Pixel 3' 'Hongyou' 'Sid'}';
    %                 figure
    %                 h = plot(g,'EdgeLabel',g.Edges.Weight);%,'NodeLabel',g.Nodes.Name);
    %                 h.LineStyle = '-';
    %                 m = size(FrameDistGraph,1);
    %                 n = size(FrameDistGraph,2);
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
    CorrelationResult(f,1) = IsHere(f,1);
    CorrelationResult(f,asscol+1) = assrow'.*IsHere(f,asscol+1);
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
cd(sequences_path)
load(subFolders(k).name+"NoOfPohneHolsers")
% if k==15
%     phonesHolders = {'Mohamed','Hongyu','Sid','Murtadha','Hansi'};
% else
% AllphonesHolders = {'Hansi','Mohamed','Hongyu','Sid','Murtadha'};
% phonesHolders = AllphonesHolders(1:num_of_phones);
% % end
% subjects = {'Hansi','Mohamed','Hongyu','Sid','Murtadha'};
if day == "20211004"
    AllphonesHolders = {'Hansi','Nicholas'};
    subjects = {'Hansi','Nicholas'};
else
    AllphonesHolders = {'Hansi','Nicholas','Bo'};
    subjects = {'Hansi','Nicholas','Bo'};
end

for ai = 1:size(CorrelationResult,1)
    phoneTP = [];
    trackIDTP = [];
    phoneFP = [];
    trackIDFP = [];
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
        TrackIDSub = string(ZedImgPosition{1,trackID(ti)}.SubName(1));
        if (CorrelationResult(ai,trackID(ti)+1)==0) && (IsHere(ai,trackID(ti)+1)==1)
            if ~any(find(TrackIDSub==phonesHolders))
                ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = "None";
                ZedImgPosition{1,trackID(ti)}.AssPhoneID(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = -1;
                TN = TN+1;
                fTN = fTN+1;
            else
                ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = "None";
                ZedImgPosition{1,trackID(ti)}.AssPhoneID(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = -1;
                FN = FN+1;
                fFN = fFN+1;
            end
        else
            phone = phonesHolders{CorrelationResult(ai,trackID(ti)+1)};
            phone = convertCharsToStrings(phone);
            
            ZedImgPosition{1,trackID(ti)}.AssignedPhone(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = phone;
            ZedImgPosition{1,trackID(ti)}.AssPhoneID(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = CorrelationResult(ai,trackID(ti)+1);
            %             ZedImgPosition{1,trackID(ti)}.Distance(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = confidenceMatrix(ai,trackID(ti)+1);
            
            %               TracktorgndPosition{1,trackID(ti)}.AssignedPhone = repmat(phone,size(TracktorgndPosition{1,trackID(ti)},1),1);
            phoneGND = [phoneGND,find(phonesHolders==TrackIDSub)];
            ModDecision{end+1,1} = k;
            ModDecision{end,2} = CorrelationResult(ai,1);
            if strcmp(TrackIDSub,phone)
                %                 subjTP(trackID(ti)) = subjTP(trackID(ti))+1;
                TP = TP+1;
                fTP = fTP+1;
                sTP(CorrelationResult(ai,trackID(ti)+1))=sTP(CorrelationResult(ai,trackID(ti)+1))+1;
                phoneTP = [phoneTP,CorrelationResult(ai,trackID(ti)+1)];
                trackIDTP = [trackIDTP,trackID(ti)];
            else
                FP = FP+1;
                sFP(CorrelationResult(ai,trackID(ti)+1))=sFP(CorrelationResult(ai,trackID(ti)+1))+1;
                phoneFP = [phoneFP,CorrelationResult(ai,trackID(ti)+1)];
                
                trackIDFP = [trackIDFP,trackID(ti)];
            end
        end
    end
    %             TN = TN + peopleInFrame-(FP+TP+FN)
    frameAcc(ai) = fTP/peopleInFrame;
    frameTP(ai) = fTP;
    %detected_people(ai)=peopleInFrame;
    cd(sequences_path)
    detected_people(ai) = NoOfPohneHolsers(ai);
    frame_number(ai) = ai;
    seq_num(ai) = k;
    Frames_FrameMinDist{ai,1}=frameAcc(ai);
    x1 = cellfun(@(v)v(1),Adindiv(:,trackID));
    x2 = cellfun(@(v)v(2),Adindiv(:,trackID));
    if size(Adindiv(:,trackID),2)>2
        x3 = cellfun(@(v)v(3),Adindiv(:,trackID));
    else
        x3 = [];
    end
    Frames_FrameMinDist{ai,2} = x1;
    Frames_FrameMinDist{ai,3} = x2;
    Frames_FrameMinDist{ai,4} = x3;
    Frames_FrameMinDist{ai,5}=detected_people(ai);
    Frames_FrameMinDist{ai,6}=phoneTP;
    Frames_FrameMinDist{ai,7}=phoneFP;
    Frames_FrameMinDist{ai,8}=phoneGND;
    Frames_FrameMinDist{ai,9}=trackIDTP;
    Frames_FrameMinDist{ai,10}=trackIDFP;
end
total_frameAcc = [total_frameAcc;seq_num',frame_number',frameAcc',detected_people',frameTP'];
%%
overallFrameAccuracy = nansum(frameAcc(1,:))/sum(~isnan(frameAcc(1,:)));
%%
% Acc = matchID/sum(~isnan(subjTP));
%%

% for bb=1:size(ZedImgPosition,2)
%     if ~isempty(ZedImgPosition{1,bb}) %&& size(ZedImgPosition{1,bb},1)>=Kframes
%         totalbb = totalbb+size(ZedImgPosition{1,bb},1);
%         [matchFrame,~] = find(ZedImgPosition{1,bb}.SubName==ZedImgPosition{1,bb}.AssignedPhone);
%         numberOfFramesToCorrectMatch = size(matchFrame,1);
%         TrackIdLengthAcc = numberOfFramesToCorrectMatch/size(ZedImgPosition{1,bb},1);
%         TrackIdLength = [TrackIdLength;size(ZedImgPosition{1,bb},1),TrackIdLengthAcc];
%     end
% end
%%
ACC = TP/sum(detected_people)

%ACC = TP/sum(detected_people)
PR = TP / (TP + FP);
REC = TP / (TP + FN);
NPV = TN / (TN + FP +FN + TP);
FS = 2 * (PR * REC) / (PR + REC);

AllSeqAcc = [AllSeqAcc,ACC]
cd(sequences_path+"/"+subFolders(k).name)
fname = subFolders(k).name+"CorrelationResultOffline.mat";
save(fname,'CorrelationResult');
%%
% fname = subFolders(k).name+"ZedSubPosition";
% cd("C:\Users\aalali\CameraIMUCorrelation/Trials/NewTrials/"+subFolders(k).name)
% save(fname,'ZedImgPosition');
% cd('C:\Users\aalali\CameraIMUCorrelation')
% %%
% ZedImgPosition22 = {};
% ZedImgPosition22 = ZedImgPosition(find(~cellfun(@isempty,ZedImgPosition)));
% ZedImgPosition222 = {};
% TTSync = timetable();
% for m=1:size(ZedImgPosition22,2)
%     ZedImgPosition222{1,m} = timetable(ZedImgPosition22{1,m}.timestamp,ZedImgPosition22{1,m}.trackID,...
%         ZedImgPosition22{1,m}.AssPhoneID,ZedImgPosition22{1,m}.TimeString);
%     TTSync = [TTSync;ZedImgPosition222{1,m}];
% end
% %TTSync = synchronize(ZedImgPosition222{:,:},'union');
% TTSync1 = timetable2table(TTSync);
% uniTimes = unique(TTSync1.Time);
% uniTimesString = unique(TTSync1.Var3);
% uniTTSync = {};
% for ut=1:numel(uniTimes)
%     uniTTSyncrows = TTSync1(TTSync1.Time==uniTimes(ut),2:3);
% %     uniTTSync{ut} = [datestr(uniTimes(ut),'dd-mm-yyyy HH_MM_SS.FFF'),' ',num2str(reshape(uniTTSyncrows{:,:}.',1,[]))];
%     uniTTSync{ut,1} = uniTimesString(ut);
%     uniTTSync{ut,2} = uniTTSyncrows{:,1}';
%     uniTTSync{ut,3} = uniTTSyncrows{:,2}';
% end
% uniTTSyncT = cell2table(uniTTSync);
% cd('C:\Users\aalali\CameraIMUCorrelation\IMUFTM_Offline_Correlation_Result')
% writetable(uniTTSyncT,subFolders(k).name+".csv")%,'Delimiter',',')
% % writetable(uniTTSyncT,subFolders(k).name+".txt",'Delimiter',',','WriteRowNames',false);