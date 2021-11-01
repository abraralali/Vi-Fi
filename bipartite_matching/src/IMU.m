ClearExeptions
% close all
% % Get a list of all files and folders in this folder.
% files = dir("C:/RANProject/WINLABData/*");
% files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and ..
% % Get a logical vector that tells which is a directory.
% dirFlags = [files.isdir];
% % Extract only those that are directories.
% subFolders = files(dirFlags);
% % Print folder names to command window.
% figure
% for k = 1 : size(subFolders,1)%for each trial
%     k
    
%     fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
    disp(sequences_path+"/"+subFolders(k).name+"/IMU/*.csv")
    seqname = subFolders(k).name;
    [~, PhoneIMUReadings] = ImportPhoneIMU(subFolders(k).name,dir(sequences_path+"/"+subFolders(k).name+"/IMU/*.csv"),src_path,sequences_path,seqname);
   
%     subplot(3,5,k)
%     for i=1:size(PhoneIMUReadings,2)
%         
%         plot(PhoneIMUReadings{1,i}.timestamp,PhoneIMUReadings{1,i}.gz+i*2)%filtered accelerometer x-axis
%         hold on
% %         plot(PhoneIMUReadings{1,i}.gy)%filtered accelerometer y-axis
% %         hold on
% %         plot(PhoneIMUReadings{1,i}.gz)%filtered accelerometer z-axis
% %         legend({'x','y','z'})
%         subs{i} = PhoneIMUReadings{1,i}.PhoneID(2);
%     end
%     xlabel('Timestamp')
%     ylabel('Acc_z')
%     title(subFolders(k).name)
%     legend(subs)
%   fname="C:\Users\aalali\CameraIMUCorrelation\Yawfigures\"+"Sequence_"+ subFolders(k).name+".png";
%   saveas(gcf,fname)  
%     figure('name',subFolders(k).name + " Gyroscope Readings")
%     for i=1:size(PhoneIMUReadings,2)
%         subplot(size(PhoneIMUReadings,2),1,i)
%         plot(PhoneIMUReadings{1,i}.x_Phonegyro)%filtered Gyroscope x-axis
%         hold on
%         plot(PhoneIMUReadings{1,i}.y_Phonegyro)%filtered Gyroscope y-axis
%         hold on
%         plot(PhoneIMUReadings{1,i}.z_Phonegyro)%filtered Gyroscope z-axis
%         legend({'x','y','z'})
%         title(unique(PhoneIMUReadings{1,i}.PhoneID))
%         hold off
%     end
    
%     figure('name',subFolders(k).name + " Orientation ")% Quaternions resulted from ecompass function
%     for i=1:size(PhoneIMUReadings,2)
%         subplot(size(PhoneIMUReadings,2),1,i)
%         plot(PhoneIMUReadings{1,i}.roll)
%         hold on
%         plot(PhoneIMUReadings{1,i}.pitch)
%         hold on
%         plot(PhoneIMUReadings{1,i}.yaw)
%         legend({'roll','pitch','yaw'})
%         title(unique(PhoneIMUReadings{1,i}.PhoneID))
%         hold off
%     end
    fname = subFolders(k).name+"IMU";
    cd(sequences_path+"/"+subFolders(k).name)
    save(fname,'PhoneIMUReadings')
    cd(src_path)
% end
%%