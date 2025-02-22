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
    FrameMinDist(FrameMinDist==0|isnan(FrameMinDist))=inf;
    Frames_FrameMinDist{time_seg,5} = FrameMinDist;

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
    [assignment,cost] = munkres(Graph_FrameMinDist);
    [assrow,asscol] = find(assignment);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1 = cellfun(@(v)v(time_seg,1),Adindiv(:,asscol));
    x2 = cellfun(@(v)v(time_seg,2),Adindiv(:,asscol));
    x3 = cellfun(@(v)v(time_seg,3),Adindiv(:,asscol));

    Frames_FrameMinDist{time_seg,2} = x1;
    Frames_FrameMinDist{time_seg,3} = x2;
    Frames_FrameMinDist{time_seg,4} = x3;

    CorrelationResult(time_seg,1) = IsHere(time_seg,1);
    CorrelationResult(time_seg,asscol+1) = assrow'.*IsHere(time_seg,asscol+1);

    mulTimeSegment = mulTimeSegment+global_winShift;


end
TotalPeople = 0;
matchID = 0;
cd(sequences_path)
load(subFolders(k).name+"NoOfPohneHolsers")
if testSequences == 1
    if day == "20211004"
    AllphonesHolders = {'Subject1','Subject6'};
    subjects = {'Subject1','Subject6','Others'};
    else
        AllphonesHolders = {'Subject1','Subject6','Subject7'};
        subjects = {'Subject1','Subject6','Subject7'};
    end
end
if testSequences == 0
    AllphonesHolders = {'Subject1','Subject2','Subject3','Subject4','Subject5'};
    subjects = {'Subject1','Subject2','Subject3','Subject4','Subject5'};
end



phonesHolders = AllphonesHolders(1:size(PhoneIMUReadings,2));
% end

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
        TrackIDSub = string(ZedImgPosition{1,trackID(ti)}.SubName(1));
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
            ZedImgPosition{1,trackID(ti)}.AssPhoneID(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = CorrelationResult(ai,trackID(ti)+1);
            %ZedImgPosition{1,trackID(ti)}.Distance(ZedImgPosition{1,trackID(ti)}.frameNo==frno) = confidenceMatrix(ai,trackID(ti)+1);

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
            else
                FP = FP+1;
                sFP(CorrelationResult(ai,trackID(ti)+1))=sFP(CorrelationResult(ai,trackID(ti)+1))+1;
                phoneFP = [phoneFP,CorrelationResult(ai,trackID(ti)+1)];

            end
            
            ModDecision{end,3} = TrackIDSub;
            ModDecision{end,4} = phone;
        end
    end
    %             TN = TN + peopleInFrame-(FP+TP+FN)
    frameAcc(ai) = fTP/peopleInFrame;
    frameTP(ai) = fTP;
    %detected_people(ai)=peopleInFrame;
    detected_people(ai) = NoOfPohneHolsers(ai);
    frame_number(ai) = CorrelationResult(ai,1);
    seq_num(ai) = k;
    Frames_FrameMinDist{ai,1}=frameAcc(ai);

    Frames_FrameMinDist{ai,8}=detected_people(ai);
    Frames_FrameMinDist{ai,9}=phoneTP;
    Frames_FrameMinDist{ai,10}=phoneFP;
    Frames_FrameMinDist{ai,11}=phoneGND;
end
total_frameAcc = [total_frameAcc;seq_num',frame_number',frameAcc',detected_people',frameTP',repmat(trial,numel(seq_num),1)];
%%
overallFrameAccuracy = nansum(frameAcc(1,:))/sum(~isnan(frameAcc(1,:)));
%%
% Acc = matchID/sum(~isnan(subjTP));
%%

for bb=1:size(ZedImgPosition,2)
    if ~isempty(ZedImgPosition{1,bb}) %&& size(ZedImgPosition{1,bb},1)>=Kframes
        totalbb = totalbb+size(ZedImgPosition{1,bb},1);
        ZedImgPosition{1,bb}.SubName = string(ZedImgPosition{1,bb}.SubName);
        [matchFrame,~] = find(ZedImgPosition{1,bb}.SubName==ZedImgPosition{1,bb}.AssignedPhone);
        numberOfFramesToCorrectMatch = size(matchFrame,1);
        TrackIdLengthAcc = numberOfFramesToCorrectMatch/size(ZedImgPosition{1,bb},1);
        TrackIdLength = [TrackIdLength;size(ZedImgPosition{1,bb},1),TrackIdLengthAcc];
    end
end
%%
%ACC = (TP+TN) / (TP+TN+FP+FN);
ACC = TP/sum(detected_people)
sequenceAccuracy{k,1} = [subFolders(k).name,ACC];
PR = TP / (TP + FP);
REC = TP / (TP + FN);
NPV = TN / (TN + FP +FN + TP);
FS = 2 * (PR * REC) / (PR + REC);
%%
%     AllTPs = AllTPs+TP;
%     AllFPs = AllFPs+FP;
%     AllFNs = AllFNs+FN;
%     AllTNs = AllTNs+TN;
%     SeqAcc = SeqAcc+ACC;
AllSeqAcc = [AllSeqAcc,ACC]
%     AllSeqFS = [AllSeqFS,FS];
%     AllSeqRecall =  [AllSeqRecall,REC];
%     AllSeqPres = [AllSeqPres,PR];
%     AllSeqNPV = [AllSeqNPV,NPV];


cd(sequences_path+"/"+subFolders(k).name)
fname = subFolders(k).name+"CorrelationResult.mat";
save(fname,'CorrelationResult');
