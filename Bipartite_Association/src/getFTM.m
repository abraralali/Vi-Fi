ClearExeptions
% close all
% cd('C:/RANProject')
% files = dir("C:/RANProject/WINLABData/*");
% files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
% % Get a logical vector that tells which is a directory.
% dirFlags = [files.isdir];
% % Extract only those that are directories.
% subFolders = files(dirFlags);
% % Print folder names to command window.
%%
% 
% for k = 1 : size(subFolders,1)
%     
    %IMPORT FTM DATA
    FTMfiles = dir(sequences_path+"/"+subFolders(k).name+"/WiFi/*.csv");
    num_depthfiles = length(FTMfiles);
    subplot(3,5,k)
    for kf = 1:num_depthfiles
        FTMData = importFTMfile(sequences_path+"/"+subFolders(k).name+"/WiFi/"+FTMfiles(kf).name);
        FTMData.timestamp = FTMData.timestamp/10^3;
        FTMData.timestamp = datetime(FTMData.timestamp,'ConvertFrom','posixtime','timezone','America/New_York','Format','HH:mm:ss.SSSSSS');
        FTMData = table2timetable(FTMData);
        FTMData.FTM = FTMData.FTM./1000;
        FTMData.FName = repmat(FTMfiles(kf).name,size(FTMData,1),1);
        AllFTMData{1,kf}=FTMData;
        plot(AllFTMData{1,kf}.timestamp,AllFTMData{1,kf}.FTM+kf*2)
        title(subFolders(k).name)
        phname{kf} = AllFTMData{1,kf}.FName(1);
        ylabel('FTM [m]')
        hold on
    end
    legend(phname)
    hold off
    
%%
%     %SPLIT DEPTH DATA BASED ON TRACK ID
%     track2start = find(isnan(FTMData.FTM));
%     FTMP = FTMData(1:track2start(1)-1,:);
%     PhoneIMUReadings{1, 1} = synchronize(PhoneIMUReadings{1, 1},FTMP,'first','nearest');
%     FTMP3b = FTMData(track2start(1)+1:track2start(2)-1,:);
%     PhoneIMUReadings{1, 2} = synchronize(PhoneIMUReadings{1, 2},FTMP3b,'first','nearest');
%     FTMP3w = FTMData(track2start(2)+1:end-1,:);
%     PhoneIMUReadings{1, 3} = synchronize(PhoneIMUReadings{1, 3},FTMP3w,'first','nearest');
%     t12 = synchronize(FTMP,FTMP3b,FTMP3w,'union','nearest');
    fname = subFolders(k).name+"FTM";
    cd(sequences_path+"/"+subFolders(k).name)
    save(fname,'AllFTMData')
    cd(sequences_path)
% end