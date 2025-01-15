%% RPS HEP routine
close all hidden
close all force
RPS_t_config_script
doImport=0;
doConcat=0;
doProcessing=0;
doDesConcat=0;
doEpoch=0;
plotECG=0;
doERP=1;
reClean=0;
doGROUP=0;
doSingle=0;
getStats=0;
doERPstat=0;
doPlots=0;
genHEP=0;
epochSize=[-1.5 1.5];
eeglab;close all;
%% Convert to SET
if doImport
    root='/Volumes/EEG_UDD/RPS/RPS_t';
    root=pwd;
    chanlocsFolder=[root,filesep,'EGI_128.ced'];
    webRoot='/Users/fjparada/Library/CloudStorage/GoogleDrive-francisco.parada@mail.udp.cl/My Drive/DATA/EEG/RPS/2023';
    webRoot='/Volumes/EEG_UDD/RPS';
    webRoot='/Users/fjparada/Library/CloudStorage/GoogleDrive-francisco.parada@mail.udp.cl/My Drive/DATA/EEG/RPS/2023';
    dataFolder=[webRoot,filesep,'1_RAW_EGI'];
    matFolder=[webRoot,filesep,'1_MAT_EGI'];
    targetFolder=[root,filesep,'2_raw-EEGLAB'];
    files2load=dir(fullfile(dataFolder,'*.raw'));
    MATfiles2load=dir(fullfile(matFolder,'*.mat'));
    ECGfileEnd='.mat';
    thisLog={};
    for fIdx=1:size(files2load,1)
        thisFile=files2load(fIdx).name;
        if ~exist([targetFolder,filesep,thisFile(1:7),'.set'])
            try
                EEG=pop_readegi([dataFolder,filesep,thisFile],[],[],'auto');
                EEG=pop_chanedit(EEG,'load',chanlocsFolder);
                thisMATfile=MATfiles2load(fIdx).name;
                thisMAT=load([matFolder,filesep,thisMATfile]);
                thisMATnames=fieldnames(thisMAT);
                ECGfile=[thisMATfile(1:7),'_',thisMATfile(9:16),'_',thisMATfile(18:end-4),'ECG'];
                % ECGfile=[thisMATnames{4}];
                try
                    thisECG=thisMAT.(ECGfile);
                    EEG.ECG=thisMAT.(thisECG);
                    EEG.ECGfile=ECGfile;
                    EEG.setname=thisFile(1:7);
                    EEG.filename=[thisFile(1:7),'.set'];
                    EEG=pop_saveset(EEG,'filename',EEG.setname,'filepath',targetFolder);
                catch ME
                    tmp=fields(thisMAT);
                    thisECG=[tmp{1},'ECG'];
                    EEG.ECG=thisMAT.(thisECG);
                    EEG.ECGfile=thisECG;
                    EEG.setname=thisFile(1:7);
                    EEG.filename=[thisFile(1:7),'.set'];
                    EEG=pop_saveset(EEG,'filename',EEG.setname,'filepath',targetFolder);
                    try
                        tmp=fields(thisMAT);
                        thisECG=[tmp{1},'ECG'];
                        EEG.ECG=thisMAT.(thisECG);
                        % ECGfile=['a',thisFile(1:6),'20230124_1231',ECGfileEnd];                        
                        EEG.ECGfile=thisECG;
                        EEG.setname=thisFile(1:7);
                        EEG.filename=[thisFile(1:7),'.set'];
                        EEG=pop_saveset(EEG,'filename',EEG.setname,'filepath',targetFolder);
                    catch DEMOLEDORA
                        thisLog{fIdx}=sprintf('File %d (%s) is corrupted',fIdx,thisFile);
                        continue
                    end
                end
            catch MEE
                thisLog{fIdx}=sprintf('File %d (%s) is corrupted',fIdx,thisFile);
                continue
            end
        else
            fprintf('File %s already imported. Skipping...\n\n\n',[thisFile(1:7),'.set'])
            continue
        end
    end
    if isempty(thisLog);thisLog='No errors while importing data';end
    save([root,filesep,'ImportLog.mat'],'thisLog')
end
%% CONCAT data for better, more stable AMICA
if doConcat
    root='/Volumes/EEG_UDD/RPS/RPS_t';
    dataFolder=[root,filesep,'2_raw-EEGLAB'];
    targetFolder=[root,filesep,'2_concat-EEGLAB'];
    tmp=dir(fullfile(dataFolder));
    participants=[];
    inx=1;
    for pId=1:size(tmp,1)
        if tmp(pId).name(1)=='.'
            continue
        else
            participants(inx).name=tmp(pId).name;
            participants(inx).folder=tmp(pId).folder;
            participants(inx).date=tmp(pId).date;
            inx=inx+1;
        end
    end
    clear tmp
    files2load=dir(fullfile(dataFolder,'*.set'));
    numSub=27; %this is hardcoded number of participants
    for thisSub=1:numSub
        ALLEEG=[];
        [tmp{1}{thisSub},theseFiles{1}{thisSub}]=LoadSetFiles(dataFolder,...
            files2load(thisSub).name);
        [tmp{2}{thisSub},theseFiles{2}{thisSub}]=LoadSetFiles(dataFolder,...
            files2load(thisSub+(numSub)).name);
        [tmp{3}{thisSub},theseFiles{3}{thisSub}]=LoadSetFiles(dataFolder,...
            files2load(thisSub+(numSub*2)).name);
        [concatData,limits,dimensions]=concatdata({tmp{1}{thisSub}.data ...
            tmp{2}{thisSub}.data tmp{3}{thisSub}.data});
        ALLEEG=tmp{1}{thisSub};
        ALLEEG.data=concatData;
        ALLEEG.concatlimits=limits;
        ALLEEG.concatdimensions=dimensions;
        ALLEEG.files=[theseFiles{1}{thisSub} theseFiles{2}{thisSub} theseFiles{3}{thisSub}];
        ALLEEG.subject=thisSub;
        ALLEEG.pnts=size(ALLEEG.data,2);
        ALLEEG.xmax=size(ALLEEG.data,2);
        for i=1:3
            ALLEEG.xmaxall{i}=tmp{i}{thisSub}.xmax;
            ALLEEG.timesall{i}=tmp{i}{thisSub}.times;
            ALLEEG.eventall{i}=tmp{i}{thisSub}.event;
            ALLEEG.ureventall{i}=tmp{i}{thisSub}.urevent;
            ALLEEG.ECGall{i}=tmp{i}{thisSub}.ECG;
            ALLEEG.ECGallfiles{i}=tmp{i}{thisSub}.ECGfile;
        end
        ALLEEG.datfile=[];
        saveName=sprintf('sub%1.2d_RPS.set',thisSub);
        saveFolder='2_concat-EEGLAB';
        if ~exist([root,filesep,saveFolder]);mkdir([root,filesep,saveFolder]);end
        ALLEEG=pop_saveset(ALLEEG,'filename',saveName,'filepath',[root,filesep,saveFolder]);
    end
end
%% Run cleaning procedures on concatenated data
if doProcessing
    RPS_EEG_processing;
end
%%
if doDesConcat
    root='/Volumes/EEG_UDD/RPS/RPS_t';
    dataFolder=[root,filesep,'5_single-subject-EEG-analysis'];
    tmp=dir(fullfile(dataFolder));
    participants=[];
    inx=1;
    for pId=1:size(tmp,1)
        if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
            continue
        else
            participants(inx).name=tmp(pId).name;
            participants(inx).folder=tmp(pId).folder;
            participants(inx).date=tmp(pId).date;
            inx=inx+1;
        end
    end
    clear tmp
    nF=3;
    for tFix=1:size(participants,2)
        files2load=dir(fullfile([participants(tFix).folder,filesep,participants(tFix).name],...
            '*RPS_cleaned_with_ICA.set'));
        EEG=pop_loadset([files2load.folder,filesep,files2load.name]);
        targetFolder=[root,filesep,'5.1_desconcat'];
        if ~exist(targetFolder);mkdir(targetFolder);end
        for nFix=1:nF
            if nFix==1
                bgn=1;
                ending=EEG.xmaxall{nFix};
            else
                bgn=EEG.xmaxall{nFix-1};
                ending=bgn+(EEG.xmaxall{nFix-1}+EEG.xmaxall{nFix});...
            end
            TMP=pop_select(EEG,'time',[bgn ending]);
            saveName=[EEG.files{nFix}];
            TMP.filename=saveName;
            TMP.filepath=targetFolder;
            TMP.condition=saveName(1:end-4);
            TMP.event=EEG.eventall{nFix};
            pop_saveset(TMP,'filename',saveName,'filepath',targetFolder);
        end
    end
    numSub=27; %this is hardcoded number of participants
end
%% Integrate ECG channels as 129 and identify cardiac events
if doEpoch
    ERPtypenames={'ERP' 'BP' 'HEP'};
    for ERPtype=3
        root='/Volumes/EEG_UDD/RPS/RPS_t';
        dataFolder=[root,filesep,'5_single-subject-EEG-analysis'];
        chanlocsName='EGI128_ECG.sfp';
        chanlocsFolderECG=[root,filesep,chanlocsName];
        tmp=dir(fullfile(dataFolder));
        participants=[];
        inx=1;
        for pId=1:size(tmp,1)
            if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                continue
            else
                participants(inx).name=tmp(pId).name;
                participants(inx).folder=tmp(pId).folder;
                participants(inx).date=tmp(pId).date;
                inx=inx+1;
            end
        end
        rootECG='/Volumes/EEG_UDD/RPS';
        dataFolderECG=[rootECG,filesep,'RAWDATA'];
        tmp=dir(fullfile([dataFolderECG],'*.mat'));
        participantsECG=[];
        inx=1;
        for pId=1:size(tmp,1)
            if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                continue
            else
                participantsECG(inx).name=tmp(pId).name;
                participantsECG(inx).folder=tmp(pId).folder;
                participantsECG(inx).date=tmp(pId).date;
                inx=inx+1;
            end
        end
        clear tmp
        for sIdx=1:size(participants,2)
            file2load=dir(fullfile([dataFolder,filesep,participants(sIdx).name],'*_cleaned_with_ICA.set'));
            thisFolder=[participants(sIdx).folder,filesep,participants(sIdx).name];
            try
                EEG=pop_loadset([thisFolder,filesep,file2load.name]);
                if floor(mean(participantsECG(sIdx).name(1:7)==EEG.setname(1:7)))
                    ECGfile2load=dir(fullfile([dataFolderECG,filesep,participantsECG(sIdx).name]));
                    ECG=load([ECGfile2load.folder,filesep,ECGfile2load.name]);
                    EEG.ECG=ECG.([EEG.ECGfile,'ECG']);
                else
                    return
                end
                targetFolder=[root filesep '6_Epoched'];if ~exist(targetFolder);mkdir(targetFolder);end
                if ERPtype==1
                    EEG_trials=pop_epoch(EEG,...
                        {'PP1L' 'PP1R' 'RK2L' 'RK2R' 'SC3L' 'SC3R'},epochSize,'epochinfo','yes');
                    finalFolder=[targetFolder filesep '6.1_ERP'];if ~exist(finalFolder);mkdir(finalFolder);end
                elseif ERPtype==2
                    try
                        EEG_trials=pop_epoch(EEG,...
                            {'btPP' 'btRK' 'btSC'},epochSize,'epochinfo','yes');
                    catch ME
                        EEG_trials=pop_epoch(EEG,...
                            {'bt'},epochSize,'epochinfo','yes');
                    end
                    finalFolder=[targetFolder filesep '6.2_BP'];if ~exist(finalFolder);mkdir(finalFolder);end
                elseif ERPtype==3
                    % downsample ECG to 500 Hz
                    EEG.data(EEG.nbchan+1,:)=EEG.ECG(1,1:2:end);
                    EEG.nbchan=size(EEG.data,1);
                    EEG=pop_chanedit(EEG,'load',{chanlocsFolderECG,'filetype','autodetect'});
                    % get ECG metrics and store in EEG
                    [data Rpeak Qpeak,Speak,Twave,Pwave IBIntervals IBItimes]=...
                        getECG(EEG,EEG.nbchan,[3 30],'slow',9,0);
                    EEG=data;
                    EEG.data=EEG.data(1:EEG.nbchan-1,:);
                    EEG.nbchan=size(EEG.data,1);
                    chanlocsName='EGI_128.ced';
                    chanlocsFolderECG=[root,filesep,chanlocsName];
                    EEG=pop_chanedit(EEG,'load',{chanlocsFolderECG,'filetype','autodetect'});
                    EEG.Rpeak=Rpeak;
                    EEG.Qpeak=Qpeak;
                    EEG.Speak=Speak;
                    EEG.Twave=Twave;
                    EEG.Pwave=Pwave;
                    EEG.IBIntervals=IBIntervals;
                    EEG.IBItimes=IBItimes;
                    if plotECG
                        figure
                        heplab_ecgplot(EEG.data(65,:),EEG.srate,Rpeak,EEG.xmin,[],9);
                        [maxY]=max(EEG.data(65,:));
                        [minY]=min(EEG.data(65,:));
                        set(gca,'ylim',[ceil(minY)-1 ceil(maxY)+10]);
                    end
                    [centralEvent,centralEventIdx,centralEventTime,centralEventOriginal,...
                        previousEvent,previousEventIdx,previousEventTime,previousEventOriginal...
                        nextEvent,nextEventIdx,nextEventTime,nextEventOriginal]=...
                        findSpecificEvents(EEG,'bt','ECG-R');

                    for sR=1:size(previousEvent,2)
                        thisEvId=previousEvent{sR}.eventIdx;
                        if size(EEG.event(thisEvId).type,2)==5
                            EEG.event(thisEvId).type=[EEG.event(thisEvId).type,'_prev'];
                        end
                    end
                    tmp=EEG;
                    tmp.event=EEG.event;
                    EEG_trials=pop_epoch(tmp,...
                        {'ECG-R_prev'},[-1.5 1.5],'epochinfo','yes');                    

                    % Eliminate trials with no HEP
                    timeThresh=190;
                    sortingEvent='bt ';
                    EEG_trials.sortVar=createSortVar(EEG_trials,timeThresh,sortingEvent);
                    EEG_trials=pop_select(EEG_trials,'trial',EEG_trials.sortVar(2,:));
                    finalFolder=[targetFolder filesep '6.3_HEP'];if ~exist(finalFolder);mkdir(finalFolder);end
                end
                pop_saveset(EEG_trials,'filename',[EEG_trials.subject,'_',...
                    ERPtypenames{ERPtype},],'filepath',finalFolder);
            catch ME
                continue
            end
        end
    end
end
%% Do extra cleaning
if reClean
    ERPtypenames={'6.1_ERP' '6.2_BP' '6.3_HEP'};
    condNames={'c01_NoTask' 'c02_RPSm' 'c03_RPSh'};
    cNam={'No Task' 'RPS AI' 'RPS human'};
    root='/Volumes/EEG_UDD/RPS/RPS_t';
    for ERPtype=2
        dataFolder=[root,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype}];
        tmp=dir(fullfile([dataFolder],'*.set'));
        participants=[];
        inx=1;
        for pId=1:size(tmp,1)
            if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                continue
            else
                participants(inx).name=tmp(pId).name;
                participants(inx).folder=tmp(pId).folder;
                participants(inx).date=tmp(pId).date;
                inx=inx+1;
            end
        end
        if isempty(participants)
            dataFolder=[pwd,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype}];
            tmp=dir(fullfile([dataFolder],'*.set'));
            participants=[];
            inx=1;
            for pId=1:size(tmp,1)
                if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                    continue
                else
                    participants(inx).name=tmp(pId).name;
                    participants(inx).folder=tmp(pId).folder;
                    participants(inx).date=tmp(pId).date;
                    inx=inx+1;
                end
            end
        end
        for cIdx=1:size(condNames,2)
            if cIdx==1;cnt=1;
            elseif cIdx==2;cnt=26;
            elseif cIdx==3;cnt=51;
            end
            for sIdx=1:size(participants,2)
                thisOne=[root,filesep,'6_Epoched',filesep,...
                    ERPtypenames{ERPtype},filesep,...
                    participants(sIdx).name];
                EEG=pop_loadset(thisOne);
                EEG=pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:128] ,...
                    'computepower',1,'linefreqs',50,'newversion',0,...
                    'normSpectrum',1,'p',0.01,'pad',2,'plotfigures',0,...
                    'scanforlines',0,'sigtype','Channels','taperbandwidth',2,...
                    'tau',100,'verb',1,'winsize',3,'winstep',3);
            end
        end

    end
end

%% Compute measures for single subjects
if doERP
    reCompute=1;
    baseTime=[-250 -50];
    muTime=[-50 1000];
    filtLims=[0.5 40];
    numSubs=25;
    nPerm=10;
    alphaVal=0.05;    
    lowFx=2;
    higFx=55;
    numFx=60;
    fx2Compute=[lowFx higFx numFx];
    ERPtypenames={'6.1_ERP' '6.2_BP' '6.3_HEP'};
    condNames={'c01_NoTask' 'c02_RPSm' 'c03_RPSh'};
    cNam={'No Task' 'RPS AI' 'RPS human'};
    root='/Volumes/EEG_UDD/RPS/RPS_t';
    root=pwd;
    for ERPtype=2
        dataFolder=[root,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype}];
        tmp=dir(fullfile(dataFolder,'*.set'));
        participants=[];
        inx=1;
        for pId=1:size(tmp,1)
            if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                continue
            else
                participants(inx).name=tmp(pId).name;
                participants(inx).folder=tmp(pId).folder;
                participants(inx).date=tmp(pId).date;
                inx=inx+1;
            end
        end
        if isempty(participants)
            dataFolder=[pwd,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype}];
            tmp=dir(fullfile([dataFolder],'*.set'));
            participants=[];
            inx=1;
            for pId=1:size(tmp,1)
                if tmp(pId).name(1)=='.' || tmp(pId).name(1)=='I'
                    continue
                else
                    participants(inx).name=tmp(pId).name;
                    participants(inx).folder=tmp(pId).folder;
                    participants(inx).date=tmp(pId).date;
                    inx=inx+1;
                end
            end
        end
        for cIdx=1:size(condNames,2)
            if cIdx==1;cnt=1;
            elseif cIdx==2;cnt=26;
            elseif cIdx==3;cnt=51;
            end
            for sIdx=1:size(participants,2)
                thisOne=[root,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype},filesep,...
                    participants(sIdx).name(1:7),'_',...
                    ERPtypenames{ERPtype}(5:end),'.mat'];                
                savePath=[root,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype}];
                if reCompute
                    file2load=dir(fullfile([dataFolder,filesep,participants(sIdx).name]));
                    fprintf('Working %s... subject %d | %s\n',ERPtypenames{ERPtype}(5:end),sIdx,file2load.name);
                    EEG=pop_loadset([file2load.folder,filesep,file2load.name]);
                    btIdx(1)=dsearchn(EEG.times',baseTime(1));
                    btIdx(2)=dsearchn(EEG.times',baseTime(2));
                    muIdx(1)=dsearchn(EEG.times',muTime(1));
                    muIdx(2)=dsearchn(EEG.times',muTime(2));
                    TFpower=[];rawTF=[];channelID=[];frequencies=[];timeVector=[];measureInfo=[];
                    fprintf('Computing TF measures... Subject %d | %s\n',sIdx,file2load.name);
                    if ERPtype~=2
                        % TF function does baseline when needed
                        [TFpower,~,~,channelID,...
                            frequencies,timeVector,measureInfoPower]=...
                            Compute_TFpower(EEG,1:EEG.nbchan,...
                            fx2Compute,0,baseTime,0,savePath);                        
                        tpow=squeeze(TFpower(1,:,:,:));
                        bPow=squeeze(mean(TFpower(1,:,:,btIdx(1):btIdx(2)),4));
                        ERD_ERS=tpow-bPow;                        
                    else
                        [TFpower,~,~,channelID,...
                            frequencies,timeVector,measureInfoPower]=...
                            Compute_TFpower(EEG,1:EEG.nbchan,...
                            fx2Compute,0,baseTime,0,savePath);
                        tpow=squeeze(TFpower(1,:,:,:));
                        bPow=squeeze(mean(TFpower(1,:,:,btIdx(1):btIdx(2)),4));
                        ERD_ERS=tpow-bPow;                        
                    end
                    % for ERP we do baseline
                    EEG=pop_rmbase(EEG,[baseTime(1) baseTime(2)],[]);
                    fprintf('Computing ERP... Subject %d | %s\n',sIdx,file2load.name);
                    [ERP(:,:,:),...
                        CI(:,:,:),...
                        PS(:,:,:)]=...
                        limo_pbci(EEG.data,nPerm,0.05);
                    fprintf('Saving files... %s with ERP, TF, rawTF computations\n',[file2load.name(1:7),...
                        '_',ERPtypenames{ERPtype}(5:end),'.mat']);                                       
                    save([root,filesep,'6_Epoched',filesep,ERPtypenames{ERPtype},filesep,...
                        file2load.name(1:7),'_',ERPtypenames{ERPtype}(5:end),'.mat'],...
                        'TFpower','ERD_ERS','channelID','frequencies','timeVector','measureInfoPower',...
                        'EEG','ERP','CI','PS','nPerm','alphaVal','-v7.3');
                    ERP=[];CI=[];PS=[];TFpower=[];rawTF=[];
                    fprintf('DONE!\n\n\n');
                else
                    fprintf('File \n\n%s\n\n already exists... Skipping\n\n',thisOne);
                    return
                end

            end
        end
    end
end