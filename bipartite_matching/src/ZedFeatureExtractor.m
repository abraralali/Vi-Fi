
% ImgPosition{1, sub} = retime(ImgPosition{1, sub},'regular','nearest','TimeStep',dt);
[C,ia,ic] = unique(ZedImgPosition{1,sub}.timestamp);
ZedImgPosition{1, sub} = ZedImgPosition{1, sub}(ia,:);
% S = withtol(ImgPosition{1, sub}.timestamp,milliseconds(10));
sysTime = ZedframeTimestamp;
tlower = datetime(ValidTimestamps.VarName3(ValidTimestamps.VarName1==subFolders(k).name), 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tupper = datetime(ValidTimestamps.VarName5(ValidTimestamps.VarName1==subFolders(k).name), 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
tf = isbetween(sysTime,tlower,tupper);
ValidFrameNo = zedFrameNo(tf);
sysTime = sysTime(tf);
sysTime.TimeZone = 'America/New_York';
% sysTime.Day = 23;
% sysTime.Year = 2020;
% sysTime.Month = 12;
% ShowingframeTimestamp = sysTime(S,:);
ZedImgPosition{1, sub}.timestamp.TimeZone = 'America/New_York';
IsHere(:,1) = ValidFrameNo;
IsHere(:,sub+1) = ismember(sysTime,ZedImgPosition{1, sub}.timestamp);
% ImgPosition{1, sub} = ImgPosition{1, sub}(S,:);
IMU = imuReadings{1,phone};
IMU.timestamp.TimeZone = 'America/New_York';
% IMU.timestamp.Day = 23;
% IMU.timestamp.Year = 2020;
% IMU.timestamp.Month = 12;
IMU.PhoneID = [];

VID = ZedImgPosition{1,sub};
VID.timestamp.TimeZone = 'America/New_York';
% VID.timestamp.Day = 23;
% VID.timestamp.Year = 2020;
% VID.timestamp.Month = 12;
% VID.subject = [];
VID.SubName = [];
% imuvisual = synchronize(VID,IMU,'first','nearest');
%
phoneFTM{1,phone} = AllFTMData{1,phone};
FTMTable = phoneFTM{1,phone};
FTMTable.FName = [];
%

% visual = ImgPosition{1, sub}(isbetween(ImgPosition{1, sub}.timestamp,lowert,uppert),:);
% imu = PhoneIMUReadings{1,phone}(isbetween(PhoneIMUReadings{1,phone}.timestamp,lowert,uppert),:);

% MFTMTable = FTMCorrelatorData(FTMCorrelatorData.trackID==sub-1,[1,3+phone]);
% MFTMTable = table2timetable(MFTMTable);
% MFTMTable.timestamp.TimeZone = 'America/New_York';

%imuvisFTM{phone,sub} = synchronize(VID,IMU,FTMTable,MFTMTable,'first','nearest');
imuvisFTM{phone,sub} = synchronize(VID,IMU,FTMTable,'first','nearest');

phonesHolders = {'Hansi','Nicholas','Bo','None'};
overlap = global_overlap;
windSize = global_windSize;
winShift = global_winShift;
Checkpoint = 1;
d1 = [];
headCorr = [];
trajdtw = [];
ftmDist = [];
Mftmdist = [];
dtwvisimu = [];
cuurImuVis = [];
trans = [];
cT = [];
VisCoor = [];
i=1;

Localframe = 1;
headWeight = global_headWeight;
trajWeight = global_trajWeight;
FTMWeight = global_FTMWeight;
mftmWeight = global_MFTMWeight;
% visual.subject(2);
% imu.PhoneID(1)
% minSize = min([size(visual,1) size(imu,1)])

while  i <= size(sysTime,1)
    %     windSize
    frame = frame+1;
    nimuTraj =[NaN,NaN];
    nvidTraj =[NaN,NaN];
    nvidheading =NaN;
    nimuyaw =NaN;
    Z =[NaN,NaN];
    
    fi = zedFrameNo(i);
    
    fwind = zedFrameNo(windSize);
    if sum(IsHere(windSize,sub+1))>0
        cuurImuVis = imuvisFTM{phone,sub}(isbetween(imuvisFTM{phone,sub}.timestamp,sysTime(1),sysTime(windSize)),:);
        
        %         mftm = cuurImuVis(:,end);
        %         mftm = timetable2table(mftm);
        %         mftm = table2array(mftm(end,2));
        %         %         mftm = mean(mftm);
        %         Mftmdist(Localframe) = mftm;%*mftmWeight;
        %         %
        %cuurImuVis.depth(isnan(cuurImuVis.depth))=0;
        
        %         [maxVar,maxInd] = max([cuurImuVis.hvar(end),cuurImuVis.fvar(end)]);
        %         if maxInd==1
        %             FTMWeight = 0.1;
        %             headWeight = 0.8;
        %         else
        %             FTMWeight = 0.8;
        %             headWeight = 0.1;
        %         end
        cuurImuVis.depth(isnan(cuurImuVis.depth)|isinf(cuurImuVis.depth))=0;
        ftmDist(Localframe) = dtw(cuurImuVis.depth,cuurImuVis.FTM)/length(cuurImuVis.depth);
        
        
        nimuTraj = [cuurImuVis.xPos cuurImuVis.yPos];
        
        nvidTraj = cuurImuVis.Fmed;
        
        [~,Z,tr] = procrustes(nvidTraj,nimuTraj);
        cuurImuVis.TransformedxPos = Z(:,1);
        cuurImuVis.TransformedyPos = Z(:,2);
        %         dx = [diff(Z(:,1));0];
        %         dy = [diff(Z(:,2));0];
        %         ZHeadingx = rad2deg(atan2(dy,dx));
        %         ZHeadingy = rad2deg(atan2(dx,dy));
        
        d = (dtw(nvidTraj(:,1),Z(:,1))/length(cuurImuVis.depth))+(dtw(nvidTraj(:,2),Z(:,2))/length(cuurImuVis.depth));
        
        d1(Localframe) = d;%*trajWeight;
        
        ntrajdtw = dtw(nvidTraj,nimuTraj);
        trajdtw(Localframe) = ntrajdtw;
        
        % standerdize heading data from both sources before using DTW
        nimuyaw = cuurImuVis.yaw;
        nimuyaw(nimuyaw<0) = nimuyaw(nimuyaw<0)+360;
        cuurImuVis.nimuyaw = nimuyaw;
        %                 if std(nimuyaw)~=0
        %                     nimuyaw = normalize(nimuyaw,'scale');
        %                 end
        %
        nvidheading = cuurImuVis.visHeadingy;
        nvidheading(nvidheading<0) = nvidheading(nvidheading<0)+360;
        cuurImuVis.nvidheading = nvidheading;% to save it after modified
        %         if length(nvidheading)>=2
        %             nvidheading = smooth((1:length(nvidheading))',nvidheading,20,'rloess');
        %         end
        cdt = dtw(nimuyaw,nvidheading)/length(cuurImuVis.depth);
        %         if std(cdt)~=0
        %             cdt = normalize(cdt,'range');
        %         end
        %         zcdt = dtw(ZHeadingy,nvidheading);
        headCorr(Localframe) = cdt;%*headWeight;
        
        
        
        
        %         AllSequencesCorrelationData(frameCounter,:) = {k,phone,sub,IsHere(i,1),datestr(sysTime(windSize)),headCorr(Localframe),d1(Localframe),ftmDist(Localframe),...
        %             cuurImuVis.hvar(end),cuurImuVis.pophvar(end),...
        %             cuurImuVis.fvar(end),cuurImuVis.popfvar(end)};
        %AllSequencesCorrelationData(frameCounter,:) = {k,phone,sub,IsHere(i,1),datestr(sysTime(windSize)),headCorr(Localframe),d1(Localframe),ftmDist(Localframe),Mftmdist(Localframe)};
%         if ~ismember(cuurImuVis.trackID(Localframe),visitedtrackIDs)
%             trackFTMphone = [trackFTMphone;k,IsHere(i,1),cuurImuVis.trackID(Localframe),ftmDist(Localframe)];
%         else
%             [trackrow,~] = find(trackFTMphone==sub-1&&trackFTMphone==IsHere(i,1));
%             trackFTMphone(trackrow,end) = ftmDist(Localframe);
%         end
       
    else
        headCorr(Localframe) = NaN;
        trajdtw(Localframe) = NaN;
        d1(Localframe) = NaN;
        ftmDist(Localframe) = NaN;
        Mftmdist(Localframe) = NaN;
    end
    
    i = i + winShift;
    windSize = windSize + winShift ;
    if windSize >= size(sysTime,1)
        windSize = size(sysTime,1);
    end
    
    frameCounter = frameCounter+1;
    Localframe = Localframe+1;
    
end
if  ~isempty(cuurImuVis)
           %idx = ismember(cuurImuVis.frameNo,3727:3799);
           %if sum(idx)>0%
           frameInfo{phone,sub} = cuurImuVis;%{phone,sub,cuurImuVis.frameNo(idx),cuurImuVis.timestamp(idx),cuurImuVis.depth(idx),...
           fff=fff+1;
           %end
end
% if sum(ismember(cuurImuVis.ftm(localframe),[3727,3751,3775,3799]))>0%plot corresponding frame information
%     figure
%     plot(seconds(cuurImuVis.timestamp-cuurImuVis.timestamp(1)),cuurImuVis.depth,'-','LineWidth',2)
%     hold on
%     plot(seconds(cuurImuVis.timestamp-cuurImuVis.timestamp(1)),cuurImuVis.ftm,'LineWidth',2)
%     legend('Depth','FTM')
%     xlabel('Time [s]')
%     %ylabel('[]')
% 
%     figure
%     plot(nvidTraj(:,1),nvidTraj(:,2),'-','LineWidth',2)
%     hold on
%     plot(Z(:,1),Z(:,2),'LineWidth',2)
%     xlabel('X')
%     ylabel('Y')
%     axis ij
%     legend({'Video Trajectory','Phone Trajectory'})
% 
%     figure
%     plot(seconds(cuurImuVis.timestamp-cuurImuVis.timestamp(1)),nvidheading,'-','LineWidth',2)
%     hold on
%     plot(seconds(cuurImuVis.timestamp-cuurImuVis.timestamp(1)),nimuyaw,'LineWidth',2)
%     legend({'Video Heading^{\circ}','Phone Heading^{\circ}'})
%     xlabel('Timestamp')
% end
d1 = d1';
headCorr = headCorr';
trajdtw = trajdtw';
ftmDist = ftmDist';
% if size(Mftmdist,2)>1
%     Mftmdist = Mftmdist';
% end
% headCorr = headCorr./nanmax(headCorr);
% trajdtw = trajdtw./nanmax(trajdtw);
% Mftmdist = Mftmdist./nanmax(Mftmdist);
% idx1 = isnan(headCorr) ;  % get the positions of inf in A
% headCorr1 = headCorr ;   % B as A
% headCorr1(idx1) = 0 ;  % repalce inf's in B with 0
%
% idx2 = isnan(trajdtw) ;  % get the positions of inf in A
% trajdtw1 = trajdtw ;   % B as A
% trajdtw1(idx2) = 0 ;  % repalce inf's in B with 0
%
% idx3 = isnan(Mftmdist) ;  % get the positions of inf in A
% Mftmdist1 = Mftmdist ;   % B as A
% Mftmdist1(idx3) = 0 ;  % repalce inf's in B with 0
%
% normheadCorr = normalize(headCorr1,'range');
% normheadCorr(idx1) = NaN;
%
% normtrajdtw = normalize(trajdtw1,'range');
% normtrajdtw(idx2) = NaN;
%
% normMftmdist = normalize(Mftmdist1,'range');
% normMftmdist(idx3) = NaN;
%
% headCorr = normheadCorr;
% trajdtw = normMftmdist;
% Mftmdist = normMftmdist

% if std(headCorr)~=0
%     headCorr = normalize(headCorr);
% end
%
% if std(trajdtw)~=0
%     trajdtw = normalize(trajdtw);
% end
%
% if std(ftmDist)~=0
%     ftmDist = normalize(ftmDist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AdNoWeights{phone,sub} = [d1,headCorr];
% if (isnan(headCorr) && ~isnan(d1)) || (headCorr == 0 && d1 ~=0)
%     headWeight = 0;
%     trajWeight = 1;
% else

% end
Adindiv{phone,sub} = [d1, headCorr, ftmDist, headCorr.*headWeight +  d1.*trajWeight, headCorr.*headWeight +  d1.*trajWeight + Mftmdist.*mftmWeight];
%Ad{phone,sub} = headCorr.*headWeight +  d1.*trajWeight + Mftmdist.*mftmWeight;% + ;% similarity of all phones with all subjects
switch trial
    case 1
        Ad{phone,sub} = d1.*trajWeight;
    case 2
        Ad{phone,sub} = headCorr.*headWeight;
    case 3
        Ad{phone,sub} = headCorr.*headWeight +  d1.*trajWeight;
    case 4
        Ad{phone,sub} = ftmDist.*FTMWeight;
    case 5
        %Ad{phone,sub} = (headCorr +  d1 + ftmDist)/3;
        Ad{phone,sub} = headCorr.*headWeight +  d1.*trajWeight + ftmDist.*FTMWeight;
        %Ad{phone,sub} = mean([headCorr,d1,ftmDist],2);
end
% Ad{phone,sub} = headCorr.*headWeight ;
alltrans{phone,sub} = trans;

%         norm_dtwvisimu = (dtwvisimu - min(dtwvisimu)) / ( max(dtwvisimu) - min(dtwvisimu) );
alld{sub} = d1(~isempty(d1)).*trajWeight + headCorr(~isempty(headCorr)).*headWeight + trajdtw(~isempty(trajdtw)).*trajdtWeight; % dissimilarity of all subjects with current phone
% alld{sub} = headCorr.*headWeight;
% allcT{phone,sub} = cT;
%         MATLABdtw =  AllSequencesCorrelationData;
%         ggg2(all(cellfun('isempty',ggg2),2),:) = [];
%         ggg2 = cell2table(ggg2);
%         writetable(ggg2,subFolders(k).name+"MATLABdtw.csv",'Delimiter',',');