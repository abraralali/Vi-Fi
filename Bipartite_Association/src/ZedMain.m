close all;
clear;
%DataPreparation

frameCounter = 1;
frameLength = [];
AllSequencesCorrelationData = {};
ALLcount = 1;
AllcountImg = 1;
fff = 1;

% samplingRate = 10;
samplingRate = 3;
specificDay = 0;
testSequences = 1;
% testSequences = 1;

%trial = 1;
for trial = 5:5
    cd('D:\RANProject\src2\src')
    if ismember(trial,[1,2])
        continue
    end
    tic
    total_frameAcc = [];
    BayesSeqAccuracy = [];
    trial
    num_of_phones = 3;
    AllTPs = 0;
    AllFPs = 0;
    AllTNs = 0;
    AllFNs = 0;
    AllSeqFS = [];
    AllSeqRecall =  [];
    AllSeqPres = [];
    AllSeqNPV = [];
    totalbb=0;
    SeqAcc = 0;
    AllSeqAcc = [];
    ModDecision = {};
    TrackIdLength = [];
    MatchedFrameDistSTD = [];
    CorrelationVariation = {};

    ZedSequencesCorrelationOnline;
%     ZedSequencesCorrelationOffline;
    %     calculateZedAccuracy2021
    %All_sequences_total_frameAcc = [All_sequences_total_frameAcc;total_frameAcc];
    All_ModalityDecision{trial} = ModDecision;
    AvgSeqAcc = SeqAcc/5;
    AllSeqAcc
    AllModalityAccuracy(trial,:) = AllSeqAcc;
    ROCmat(trial,:) = [AllTPs,AllFPs,AllTNs,AllFNs,AllTPs+AllFPs+AllTNs+AllFNs,totalbb,mean(AllSeqPres),mean(AllSeqRecall),mean(AllSeqFS),mean(AllSeqAcc)];
    AccuracyMatrix(trial,:) = [AvgSeqAcc,mean(AllSeqFS),mean(AllSeqRecall),mean(AllSeqPres)];
    BoxPlotMat{trial} = [AllSeqPres',AllSeqRecall',AllSeqNPV',AllSeqFS',AllSeqAcc'];
    %trial = trial+1;
    toc
    % save("Frames_FrameMinDistTrial"+num2str(trial),'Frames_FrameMinDist')
    cd('D:\RANProject\src2\src')
end
%%
% figure;
% als = All_sequences_total_frameAcc;
% x = linspace(0,5);
%  hold on;
%  plot(als{1,1})
%  plot(als{1,2})
%  plot(als{1,3})
%  plot(als{1,4})
%  plot(als{1,5})
%  legend({'Trajectory','Heading','Trajectory+Heading','FTM','All'})
%  xlabel('Sequence Number')
%  ylabel('Accuracy')
%  ylim([0 1])
%  grid on
%%
% filename="ModalitiesAccuracy2021"+datestr(datetime('now')); 
% save(filename,'All_sequences_total_frameAcc');
% save('ModalityDecisionCov','All_ModalityDecision')
%sum(All_sequences_total_frameAcc{1, 1}(:,5))/sum(All_sequences_total_frameAcc{1, 1}(:,4))
%%
% ROCmat
%  TrackIdLength = sortrows(TrackIdLength)
%%

% figure
% avgTrackAcc = [];
% maxTrackAcc = [];
% minTrackAcc = [];
% uniTrackIdLength = unique(trackslength(:,1));
% for tt=1:numel(uniTrackIdLength)
%     uniTrackIdLength(tt)
%     numel(TrackIdLength(TrackIdLength(:,1)==uniTrackIdLength(tt)));
%     avgTrackAcc = [avgTrackAcc;std(TrackIdLength(TrackIdLength(:,1)==uniTrackIdLength(tt),2))];
%     maxTrackAcc = [maxTrackAcc;max(TrackIdLength(TrackIdLength(:,1)==uniTrackIdLength(tt),2))-avgTrackAcc];
%     minTrackAcc = [minTrackAcc;avgTrackAcc-min(TrackIdLength(TrackIdLength(:,1)==uniTrackIdLength(tt),2))];
% end
% plot(uniTrackIdLength,avgTrackAcc,'LineWidth',2)
% % hold on
% % er = errorbar(uniTrackIdLength,avgTrackAcc,minTrackAcc,maxTrackAcc);
% % er.Color = [1 0 0];
% % er.LineStyle = 'none';
% title('Average Correlation Accuracy')
% xlabel('TrackId Length')
% ylabel('Avg. Accuracy')
% %%
% summat = [];
% for bbx=1:size(trackIDHisLength,1)
%     for bbxi=1:size(trackIDHisLength,2)
%         if ~isempty(trackIDHisLength(bbx,bbxi))
%             summat = [summat;trackIDHisLength(bbx,bbxi)];
%         end
%     end
% end
% summat = rmmissing(cell2mat(summat));
% figure
% ecdf(summat)
% %%
% figure
% histogram(summat,10)
% grid on
% xlabel('Track Length')
% ylabel('Number of Tracks')
% %%
% % data = mean(nonZeroROCmat);
% % errhigh = max(nonZeroROCmat)-mean(nonZeroROCmat);
% % errlow  = mean(nonZeroROCmat)-min(nonZeroROCmat);
% % x = [1];
% % bar(x,data)
% % hold on
% % er = errorbar(x,data,errlow,errhigh);
% % er.Color = [0 0 0];
% % er.LineStyle = 'none';
% % xlabel('K [Frames]')
% % ylabel('Avg. Accuracy')
% % g = 5*1:5:5*size(nonZeroROCmat,2);
% for bbx=1:size(BoxPlotMat,2)
%     summat(:,bbx) = BoxPlotMat{1,bbx}(:,5);
% end
% %%
% boxplot(summat)
% % hold off
% % nonZeroROCmat = sort(nonZeroROCmat,1,'descend')
% % %%
% % % The x-axis showing 1 – specificity (= false positive fraction = FP/(FP+TN))
% % % The y-axis showing sensitivity (= true positive fraction = TP/(TP+FN))
% %
% % % FPR = nonZeroROCmat(:,2)./(nonZeroROCmat(:,2)+nonZeroROCmat(:,4));
% % % TPR = nonZeroROCmat(:,1)./(nonZeroROCmat(:,1)+nonZeroROCmat(:,3));
% % plot(rIdcs,nonZeroROCmat,'LineWidth',2)
% % % xticks(Idcs(1):Idcs(end))
% xlabel('Length of Execluded TrackID')
% ylabel('Accuracy Sequences')
% %%
% %Track Lengths Histogram
% %%
% avgPR = [];
% avgRec = [];
% avgNPV = [];
% avgFS = [];
% avgAcc = [];
%%
% for bb=1:size(BoxPlotMat,2)
%     avgPR = [avgPR,mean(BoxPlotMat{1,bb}(:,1))];
%     avgRec = [avgRec,mean(BoxPlotMat{1,bb}(:,2))];
%     avgNPV = [avgNPV,mean(BoxPlotMat{1,bb}(:,3))];
%     avgFS = [avgFS,mean(BoxPlotMat{1,bb}(:,4))];
%     avgAcc = [avgAcc,mean(BoxPlotMat{1,bb}(:,5))];
%     nlessmavgAcc = [avgAcc,mean(BoxPlotMat{1,bb}(:,5))];
% end
% figure
% bar([avgPR',avgRec',avgNPV',avgFS',avgAcc'])
% legend({'Precision','Recall','TNR','F1-Score','Accuracy'})
% xlabel('No. Of Connected Phones')
%%
% %%
% figure;
% plot(AllSeqAcc);
% plot(AllSeqAcc,'-o','LineWidth',2)
% yyaxis left
% ylim([0 1])
% xlim([1 numel(AllSeqAcc)])
% ylabel('Sequence Accuracy')
% hold on;
% hline(mean(AllSeqAcc))
% hold on
% yyaxis right
% plot(frameLength,'-*','LineWidth',2)
% ylabel('Video Length [frames]')
% xlabel('Sequence Number')
%%
% figure
% mtlen = [];
% for tl=1:size(trackslength,2)
% %    subplot(3,5,tl)
% %    boxplot(trackslength{1,tl}(trackslength{1,tl}>0))
% %    xxx = trackslength{1,tl};
% %    title(num2str(tl))
% %    ylim([0 100])
% %    ylabel('Track Length')
% %    grid on
%     mtlen = [mtlen,mean(trackslength{1,tl}{1,:})];
% end
% plot(mtlen)
% %%
% % Zed Accuracy:
% figure;
% plot(zedAccuracy,'r');
% plot(zedAccuracy,'-or','LineWidth',2)
% ylim([0 1])
% xlim([1 numel(zedAccuracy)])
% hold on;
% hline(mean(zedAccuracy))
% hold on
% grid on
% ylabel('Tracktor++ Detection Accuracy')
% xlabel('Sequence Number')
% %%
%%
% ggg2 =  AllSequencesCorrelationData;
% ggg2(all(cellfun('isempty',ggg2),2),:) = [];
% %%
% ggg2 = cell2table(ggg2);
% writetable(ggg2,"ZedAllSequencesCorrelationData.csv",'Delimiter',',');