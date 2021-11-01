ClearExeptions
close all
% sequences_path = "C:/RANProject/WINLABData";
% src_path = "C:/RANProject/src";

AllIMU = {};
AllPositions = {};
PhoneIMUReadings = {};
fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
cd(sequences_path+"/"+subFolders(k).name)


IMUfname = subFolders(k).name+"IMU.mat";
AllIMU = load(IMUfname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure('name',subFolders(k).name)
alltable = [];
for i=1:size(AllIMU.PhoneIMUReadings,2)
    cd(src_path)
 
    [stepDetector,walkingFrequency] = DetectSteps(AllIMU.PhoneIMUReadings{1,i}.timestamp,AllIMU.PhoneIMUReadings{1,i}.gm);
    
    cd(sequences_path+"/"+subFolders(k).name)
   
    AllIMU.PhoneIMUReadings{1,i}.stepDetector = stepDetector;
    AllIMU.PhoneIMUReadings{1,i}.walkingFrequency = walkingFrequency;
    
    PhoneIMUReadings{1,i} = AllIMU.PhoneIMUReadings{1,i};
    
%     subplot(size(AllIMU.PhoneIMUReadings,2),1,i)
    % %         t2 = seconds(AllIMU.PhoneIMUReadings{1,i}.timestamp-AllIMU.PhoneIMUReadings{1,i}.timestamp(1));
%     plot(AllIMU.PhoneIMUReadings{1,i}.timestamp, AllIMU.PhoneIMUReadings{1,i}.gm,'LineWidth',2);
%     hold on
    % %         t1 = seconds(AllIMU.PhoneIMUReadings{1,i}.timestamp(AllIMU.PhoneIMUReadings{1,i}.stepDetector==1)-AllIMU.PhoneIMUReadings{1,i}.timestamp(1));
    walkingS = find(AllIMU.PhoneIMUReadings{1,i}.stepDetector==1);
%     plot(AllIMU.PhoneIMUReadings{1,i}.timestamp(walkingS),AllIMU.PhoneIMUReadings{1,i}.gm(walkingS),'.');
    %
    % %         plot(AllIMU.PhoneIMUReadings{1,i}.timestamp, AllIMU.PhoneIMUReadings{1,i}.mg,'DisplayName','Gyro mg');
    %axis normal;
    %         hold on
    %         plot(AllIMU.PhoneIMUReadings{1,i}.timestamp, AllIMU.PhoneIMUReadings{1,i}.yawTurns*200,'DisplayName','YawTurns');
    %         hold on
    %         plot(AllIMU.PhoneIMUReadings{1,i}.timestamp, AllIMU.PhoneIMUReadings{1,i}.yaw,'DisplayName','Yaw');
    %         title(AllIMU.PhoneIMUReadings{1,i}.PhoneID(end))
    %         hold on
    AllIMU.PhoneIMUReadings{1,i}.PhoneID = replace(AllIMU.PhoneIMUReadings{1,i}.PhoneID," ","");
    %fprintf('IMUPhoneID #%s = \n',AllIMU.PhoneIMUReadings{1,i}.PhoneID(end))
    %fprintf('Start Walking #%s = \n',AllIMU.PhoneIMUReadings{1,i}.timestamp(AllIMU.PhoneIMUReadings{1,i}.stepDetector==1));
    %         if  subFolders(k).name=="20200805_185546"
    %             for n=1:size(SubjTurns,2)
    %                 fprintf('IMUPhoneID #%s = \n',AllIMU.PhoneIMUReadings{1,i}.PhoneID(end))
    % %                 AllIMU.PhoneIMUReadings{1,i}.PhoneID = lower(AllIMU.PhoneIMUReadings{1,i}.PhoneID);
    % %                 SubjTurns{1,n}.PhoneId = lower(SubjTurns{1,n}.PhoneId);
    %                 if string(AllIMU.PhoneIMUReadings{1,i}.PhoneID(end))==string(SubjTurns{1,n}.PhoneId(end))
    %                     cd('C:\Users\aalali\CameraIMUCorrelation')
    %
    %                     % Get Accuracy of the Detector:
    %                     DetTurTimes = AllIMU.PhoneIMUReadings{1,i}.timestamp(find(AllIMU.PhoneIMUReadings{1,i}.gyroTurns==1));
    %                     tupper = PhoneIMUReadings{1,i}.timestamp(end);
    %                     tlower = PhoneIMUReadings{1,i}.timestamp(1);
    %                     GNDTurTimes = SubjTurns{1,n}.timestamp(isbetween(SubjTurns{1,n}.timestamp,tlower,tupper));
    %                     [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes);
    %                     GyroTurnAccPerSub{1,i} = {TP,FP,FN,F1,Turnsdelay};
    %                     vline(GNDTurTimes);
    % %                     DetTurTimes = AllIMU.PhoneIMUReadings{1,i}.timestamp(AllIMU.PhoneIMUReadings{1,i}.yawTurns==1);
    % %                     GNDTurTimes = SubjTurns{1,n}.timestamp;
    % %                     [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes);
    % %                     YawTurnAccPerSub{1,i} = {TP,FP,FN,F1,Turnsdelay};
    % %                     cd("./Trials/"+subFolders(k).name)
    %                     break
    %                 end
    %             end
    %         end
    %         legend('Z-axis Angular Velocity [deg/sec]','Turns Detected by Gyro')
    % axis normal;
    % hold off
end
%ax = findobj(f1,'Type','Axes');
%linkaxes(ax, 'x');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('name',subFolders(k).name)
%     for j=1:size(AllPositions.ImgPosition,2)
% %         AllPositions.ImgPosition{1,j} = AllPositions.ImgPosition{1,j}((AllPositions.ImgPosition{1,j}.timestamp >= min(AllIMU.PhoneIMUReadings{1,i}.timestamp))&...
% %                                 (AllPositions.ImgPosition{1,j}.timestamp <= max(AllIMU.PhoneIMUReadings{1,1}.timestamp)),:);
%         cd('C:\Users\aalali\CameraIMUCorrelation')
%
%         dt = milliseconds(100);
%         AllPositions.ImgPosition{1,j} = retime(AllPositions.ImgPosition{1,j},'regular','previous','TimeStep',dt);
%
%         XYChangeTurns = DetectXYChangeTurns(AllPositions.ImgPosition{1,j}.Fmed);
%         XYSlopeTurns  = DetectXYSlopTurns(AllPositions.ImgPosition{1,j}.Fmed);
%         [FX,FY] = XYGradient(AllPositions.ImgPosition{1,j}.Fmed);
%         [XYGradientTurns,XYGradientHead] = DetectXYGradientTurns(AllPositions.ImgPosition{1,j}.Fmed,AllPositions.ImgPosition{1,j}.timestamp);
% %         CamstepDetector = DetectStepsCam(AllPositions.ImgPosition{1,j}.Velocity);
% %         XYResdTurns = DetectXYTurnsResiduals(AllPositions.ImgPosition{1,j}.Fmed);
%
%         cd("./Trials/NewTrials/"+subFolders(k).name)
%
% %         if size(CamstepDetector,1)==1
% %             AllPositions.ImgPosition{1,j}.CamstepDetector = CamstepDetector';
% %         else
% %             AllPositions.ImgPosition{1,j}.CamstepDetector = CamstepDetector;
% %         end
% %
%         if size(XYGradientTurns,1)<10
%             AllPositions.ImgPosition{1,j}.XYGradientTurns = XYGradientTurns';
%             AllPositions.ImgPosition{1,j}.XYGradientHead = XYGradientHead';
%         else
%             AllPositions.ImgPosition{1,j}.XYGradientTurns = XYGradientTurns;
%             AllPositions.ImgPosition{1,j}.XYGradientHead = XYGradientHead;
%         end
%
% %         if size(XYSlopeTurns,1)==1
% %             AllPositions.ImgPosition{1,j}.XYSlopeTurns = XYSlopeTurns';
% %         else
% %             AllPositions.ImgPosition{1,j}.XYSlopeTurns = XYSlopeTurns;
% %         end
%
% %         if size(XYResdTurns,1)==1
% %             AllPositions.ImgPosition{1,j}.XYResdTurns = XYResdTurns';
% %         else
% %             AllPositions.ImgPosition{1,j}.XYResdTurns = XYResdTurns;
% %         end
%
% %         AllPositions.ImgPosition{1,j}.FX = FX;
% %         AllPositions.ImgPosition{1,j}.FY = FY;
%
%         ImgPosition{1,j} = AllPositions.ImgPosition{1,j};
%
% %         subplot(size(AllPositions.ImgPosition,2),1,j)
% %         plot(AllPositions.ImgPosition{1,j}.timestamp, AllPositions.ImgPosition{1,j}.XYGradientTurns*max(AllPositions.ImgPosition{1,j}.Dist));
% %         hold on
% %         plot(AllPositions.ImgPosition{1,j}.timestamp,AllPositions.ImgPosition{1,j}.Dist);
% %         hold on
% %         title(AllPositions.ImgPosition{1,j}.subject(end))
% %         hold on
% %         if  subFolders(k).name=="20200805_185546"
% %             cd('C:\Users\aalali\CameraIMUCorrelation')
% %             GetSubjectTurns;
% % %             AllIMU.PhoneIMUReadings{1,i}.timestamp.Hour = AllIMU.PhoneIMUReadings{1,i}.timestamp.Hour-4;
% %             for m=1:size(SubjTurns,2)
% %                 if string(AllPositions.ImgPosition{1,j}.subject(end))==string(SubjTurns{1,m}.subject(end))
% %                         vline(SubjTurns{1,m}.timestamp);
% %                         % Get Accuracy of the Detector:
%                         DetTurTimes = AllIMU.PhoneIMUReadings{1,i}.timestamp(AllPositions.ImgPosition{1,j}.XYResdTurns==1);
%                         GNDTurTimes = SubjTurns{1,m}.timestamp;
%                         [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes);
%                         XYResTurnAccPerSub{1,j} = {TP,FP,FN,F1,Turnsdelay};
% %
%                         DetTurTimes = AllPositions.ImgPosition{1,j}.timestamp(AllPositions.ImgPosition{1,j}.XYSlopeTurns==1);
%                         GNDTurTimes = SubjTurns{1,m}.timestamp;
%                         [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes);
%                         XYSlopeTurnsAccPerSub{1,j} = {TP,FP,FN,F1,Turnsdelay};
% %
%                         DetTurTimes =
% AllPositions.ImgPosition{1,j}.timestamp(find(AllPositions.ImgPosition{1,j}.XYGradientTurns==1));
%                         tupper =
% AllPositions.ImgPosition{1,j}.timestamp(end); % % tlower =
% AllPositions.ImgPosition{1,j}.timestamp(1); % % GNDTurTimes =
% SubjTurns{1,n}.timestamp(isbetween(SubjTurns{1,n}.timestamp,tlower,tupper));
%                         GNDTurTimes = SubjTurns{1,m}.timestamp; % %
% [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes); % %
% XYGradientTurnsAccPerSub{1,j} = {TP,FP,FN,F1,Turnsdelay};
% %
%                         DetTurTimes = AllIMU.PhoneIMUReadings{1,i}.timestamp(AllPositions.ImgPosition{1,j}.XYChangeTurns==1);
%                         GNDTurTimes = SubjTurns{1,m}.timestamp;
%                         [TP,FP,FN,F1,Turnsdelay] = calcAccuracy(DetTurTimes,GNDTurTimes);
%                         XYChangeTurnAccPerSub{1,j} = {TP,FP,FN,F1,Turnsdelay};
%                 end
%             end
%             legend('x','y','Detected Turns by Gradient')
%             cd("./Trials/"+subFolders(k).name)
%             fname = subFolders(k).name+"XYSlopeTurnsAcc";
%             save(fname,'XYSlopeTurnsAccPerSub')
%             fname = subFolders(k).name+"GyroTurnAcc";
%             save(fname,'GyroTurnAccPerSub')
%             fname = subFolders(k).name+"GradTurnAcc";
%             save(fname,'XYGradientTurnsAccPerSub')
%             fname = subFolders(k).name+"ChangeTurnAcc";
%             save(fname,'XYChangeTurnAccPerSub')
%             fname = subFolders(k).name+"ResTurnAcc";
%             save(fname,'XYResTurnAccPerSub')
% %
% %         end
% %         hold off
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure('name',subFolders(k).name)
%         subplot(size(AllPositions.ImgPosition,2),1,j)
%         plot(AllPositions.ImgPosition{1,j}.timestamp, AllPositions.ImgPosition{1,j}.XYSlopeTurns*200,'DisplayName','XYSlopeTurns');
%         hold on
%         plot(AllPositions.ImgPosition{1,j}.timestamp,AllPositions.ImgPosition{1,j}.Fmed);
%         hold on
%         title(unique(AllPositions.ImgPosition{1,j}.subject))
%         hold on
%         if  strcmp(subFolders(k).name,'20200805_185546')
%             cd('C:\Users\aalali\CameraIMUCorrelation')
%             GetSubjectTurns;
%             for m=1:size(SubjTurns,2)
%                 if strcmp(cellstr(unique(AllPositions.ImgPosition{1,j}.subject)),cellstr(unique(SubjTurns{1,m}.subject)))
%                     fprintf('Subject # = %s\n',unique(SubjTurns{1,m}.subject))
%                     vline(SubjTurns{1,m}.timestamp);
%                 end
%             end
%             cd("./Trials/"+subFolders(k).name)
%         end
%         hold off
%     end
%     fname = subFolders(k).name+"SubPosition";
%     save(fname,'ImgPosition');
fname = subFolders(k).name+"IMU";
save(fname,'PhoneIMUReadings')
cd(src_path)
% end