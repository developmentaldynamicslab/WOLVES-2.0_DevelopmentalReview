nObjects = 6; %%total number of objects
maxTR = 24; % total # of training trials
t_max = historyDuration;%floor((4000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration

sprintf('running production task file')

colorPropertiesIndices=1:6;
colorLabelIndices=[3 4 5 6 7 8];%
shapePropertiesIndices=9:14;
shapeLabelIndices=[11 12 13 14 15 16];
nLabels= 12;

sim.init();

%% Initialize

if (mode == 1 ) % GUI initilaizations and looking history updation
    gui.init();
    delete(gui.visualizations{1}.axesHandle);
    gui.visualizations{1}.axesProperties{2} = [-t_max, 10];
    gui.visualizations{1}.plotProperties{1}{end} = 0:-1:-t_max+1;
    gui.visualizations{1}.plotProperties{2}{end} = 0:-1:-t_max+1;
    gui.visualizations{1}.plotProperties{3}{end} = 0:-1:-t_max+1;
    gui.visualizations{1}.plotProperties{4}{end} = 0:-1:-t_max+1;
    gui.visualizations{1}.plotProperties{5}{end} = 0:-1:-t_max+1;
    gui.visualizations{1}.init(gui.figureHandle);
end
hStore=sim.getElement('wf1').h;
hStoreW=sim.getElement('word').h;

%% parfor loop -- use parfor for mode 2, use for for mode 1
parfor subject = 1:numSubjects 

    f1=repmat(1:nObjects,1,maxTR/nObjects);
    f2=repmat(1:nObjects,1,maxTR/nObjects);

    objectList{subject}=zeros(numel(f1),2);

    objectList{subject} = [f1(randperm(numel(f1)))' f2(randperm(numel(f2)))'];

    sim2 = 5;
    if (mode == 1 )
        sim2 = sim;
    elseif (mode == 0)
        sim2 = sim;
    elseif (mode == 2)
        sim2 = sim.copy();
    end

    %Note: should probably add 'history motor' and 'history inpn' for completeness
    handle_historyL=sim2.getElement('history lookL'); %reset lookingduration history variable
    handle_historyL.timeSlots = t_max;
    handle_historyR=sim2.getElement('history lookR');
    handle_historyR.timeSlots = t_max;
    handle_historyLWB=sim2.getElement('history lookLWB'); %reset lookingduration history variable
    handle_historyLWB.timeSlots = t_max;
    handle_historyRWB=sim2.getElement('history lookRWB');
    handle_historyRWB.timeSlots = t_max;
    handle_historyC=sim2.getElement('history lookC');
    handle_historyC.timeSlots = t_max;


    %% Make_WordFeatureMapping
    gaus = @(x,mu,sig,amp)amp*exp(-(((x-mu).^2)/(2*sig.^2)));
    memTraceRange=15; %10;
    for i=1:nLabels/2 %label-feature mappings
        x = 1:fieldSize_ftr;
        sig = memTraceRange;  % unknown width
        amp = memTraceStrength(subject);  % unknown amplitude

        mu = Property(1,colorPropertiesIndices(i));  % unknown mean
        y = gaus(x,mu,sig,amp);
        sim2.getElement('hwf1').output(:,colorLabelIndices(i))=y;

        mu = Property(2,shapePropertiesIndices(i)-8);  % unknown mean
        y = gaus(x,mu,sig,amp);
        sim2.getElement('hwf2').output(:,shapeLabelIndices(i))=y;

        % %             sim2.getElement('hwf1').output(Property(1,colorPropertiesIndices(i)):Property(1,colorPropertiesIndices(i))+memTraceRange,colorLabelIndices(i))=memTraceStrength;
        % %             sim2.getElement('hwf2').output(Property(2,shapePropertiesIndices(i)-8):Property(2,shapePropertiesIndices(i)-8)+memTraceRange,shapeLabelIndices(i))=memTraceStrength;
    end

    memTraceStrengthJ=memTraceStrengthJbase(randperm(length(memTraceStrengthJbase)));
    createjitter = 1;
    while createjitter
        tempcolorPropertiesIndices = colorPropertiesIndices(randperm(length(colorPropertiesIndices)));
        tempshapePropertiesIndices = shapePropertiesIndices(randperm(length(shapePropertiesIndices)));
        diffcolor = tempcolorPropertiesIndices - colorPropertiesIndices;
        diffshape = tempshapePropertiesIndices - shapePropertiesIndices;
        if size(find(diffcolor == 0),2) == 0 && size(find(diffshape == 0),2) == 0
            createjitter=0;
        end
    end

    for i=1:nLabels/2 %label-feature mappings
        x = 1:fieldSize_ftr;
        sig = memTraceRange;  % unknown width
        amp = memTraceStrengthJ(1,i);  % unknown amplitude

        mu = Property(1,tempcolorPropertiesIndices(i));  % unknown mean
        y = gaus(x,mu,sig,amp);
        z = sim2.getElement('hwf1').output(:,colorLabelIndices(i))';
        sim2.getElement('hwf1').output(:,colorLabelIndices(i))=z+y;

        mu = Property(2,tempshapePropertiesIndices(i)-8);  % unknown mean
        y = gaus(x,mu,sig,amp);
        z = sim2.getElement('hwf2').output(:,shapeLabelIndices(i))';
        sim2.getElement('hwf2').output(:,shapeLabelIndices(i))=z+y;
    end


    %% Production
    visuals_On = floor(500/scale_factor);    visuals_Off = floor(3750/scale_factor); % visual display on off timing
    word1_On  = floor(50/scale_factor);     word1_Off = floor(3750/scale_factor); % word 1 on off timing

    % DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    for i = 1 : nFeatures
        n=num2str(i);
        handle_hcon_f{subject}(i)=sim2.getElement(['hcon_f' n]);
        handle_hwm_f{subject}(i)=sim2.getElement(['hwm_f' n]);
        handle_hwf{subject}(i)=sim2.getElement(['hwf' n]);
        handle_hwm_c{subject}(i)=sim2.getElement(['hwm_c' n]);
    end
    handle_hcon_s{subject}=sim2.getElement('hcon_s');
    handle_hwm_s{subject}=sim2.getElement('hwm_s');
    handle_hword{subject}=sim2.getElement('hword');

    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_wf{subject}=zeros(fieldSize_wd);
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    respLocTmp{subject}=zeros(maxTR,1);
    training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%

    %% RUN PRODUCTION TRIALS
    for tr=1:maxTR %n trials loop

        %[subject tr]

        sim2.setElementParameters({'wf1' 'wf2' 'word'},{'h' 'h' 'h'},[hStore hStore hStoreW]);
        sim2.t =sim2.tZero;

        if tr > 1
            for i = 1 : nFeatures
                n=num2str(i);
                handle_hcon_f{subject}(i).output = squeeze(savestate_hcon_f{subject}(i,:));
                handle_hwm_f{subject}(i).output = squeeze(savestate_hwm_f{subject}(i,:));
                handle_hwf{subject}(i).output = squeeze(savestate_hwf{subject}(i,:,:));
                handle_hwm_c{subject}(i).output = squeeze(savestate_hwm_c{subject}(i,:,:));
            end
            handle_hcon_s{subject}.output = savestate_hcon_s{subject};
            handle_hwm_s{subject}.output = savestate_hwm_s{subject};
            handle_hword{subject}.output = savestate_hword{subject};
        end
        respMade=0;
        respTime=0;
        respThresh=10;
        while sim2.t <= t_max % run the ime loop of a trial

            t = sim2.t;

            %%%%This is now the same as what is happening in the DCCS
            if t>visuals_On
                if max(sim2.getElement('word').output)<.9
                    hprev=sim2.getElement('wf1').h;
                    hnew=hprev+.006;
                    sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hnew hnew]);
                    hwordprev=sim2.getElement('word').h;
                    hwordnew=hwordprev+.0007;
                    sim2.setElementParameters('word','h',hwordnew);

                else
                    hprev=sim2.getElement('wf1').h;
                    hnew=hprev-.03;
                    if hnew > hStore
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hnew hnew]);
                    else
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hStore hStore]);
                    end
                    sim2.setElementParameters('word','h',hStoreW);
                end
            end

            if isalmost(t, visuals_On,tolerance) %present visual stimuli
                for f=1:nFeatures

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set feature dimension value
                        },{'positionY'}, ...
                        {Property(f,colorPropertiesIndices((objectList{subject}(tr,f))))}); %Make list of features to present

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set spatial dimension value
                        }, ...
                        {'positionX'}, {spt_C});   %Make value for center of field

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set left right stimuli up
                        }, ...
                        {'amplitude'}, vstrength * ones(1, 1));
                end
            end
            if isalmost(t, visuals_Off,tolerance) %present visual stimuli
                for f=1:nFeatures

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set feature dimension value
                        },{'positionY'}, ...
                        {Property(f,colorPropertiesIndices((objectList{subject}(tr,f))))}); %Make list of features to present

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set spatial dimension value
                        }, ...
                        {'positionX'}, {spt_C});   %Make value for center of field

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Centre'), ... % set left right stimuli up
                        }, ...
                        {'amplitude'}, 0);
                end
            end

            if isalmost(t, word1_On,tolerance) %present audio stimuli
                for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                    sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                        {(colorLabelIndices(wi)),wstrength(subject)});
                end
            end
            if isalmost(t, word1_Off,tolerance) %end audio stimuli
                for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                    sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                        {(colorLabelIndices(wi)),0});
                end
            end

            if isalmost(t, t_max,tolerance) % save loking data and memory activity states
                for i = 1 : nFeatures
                    n=num2str(i);
                    savestate_hcon_f{subject}(i,:) = sim2.getComponent(['hcon_f' n], 'output');
                    savestate_hwm_f{subject}(i,:) = sim2.getComponent(['hwm_f' n], 'output');
                    savestate_hwf{subject}(i,:,:) = sim2.getComponent(['hwf' n], 'output');
                    savestate_hwm_c{subject}(i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
                end
                savestate_hcon_s{subject} = sim2.getComponent('hcon_s', 'output');
                savestate_hwm_s{subject} = sim2.getComponent('hwm_s', 'output');
                savestate_hword{subject} = sim2.getComponent('hword', 'output');
                savestate_historyL{subject}(tr,:) = flipud(sim2.getComponent('history lookL', 'output'));
                savestate_historyR{subject}(tr,:) =flipud(sim2.getComponent('history lookR', 'output'));
            end

            if (mode == 1)
                if gui.quitSimulation
                    gui.close();
                    break;
                end
                if (mod(t,gui_speed)==0)
                    gui.step();
                else
                    sim2.step();
                end
            else
                sim2.step();
            end

            if respMade==0 %index response locations
                if max(sim2.getElement('word').output)>.5
                    respTime=respTime+1
                else
                    respTime=0;
                end
                if respTime==10
                    savestate_wf{subject}(tr,:)=sim2.getComponent('word','output');
                    respPeak=find(sim2.getElement('word').output>.5);
                    respLocTmp{subject}(tr)=mean(respPeak);
                    respRT(subject,tr)=sim2.t;
                    respMade=1;
                end
            end
        end % Time sim2.t
    end % Trial maxTR loop


    %% Wrap Up
    try %save the data
        OutName = [simName,num2str(subject),'_test.mat'];
        matTrain=matfile(OutName,'writable',true);

        matTrain.hcon_f = savestate_hcon_f{subject};
        matTrain.hwm_f = savestate_hwm_f{subject};
        matTrain.hwm_c = savestate_hwm_c{subject};
        matTrain.hwf = savestate_hwf{subject};
        matTrain.hcon_s = savestate_hcon_s{subject};
        matTrain.hwm_s = savestate_hwm_s{subject};
        matTrain.hword = savestate_hword{subject};
        matTrain.respLoc=respLocTmp{subject};
        matTrain.corrLoc=objectList{subject};
        matTrain.wf = savestate_wf{subject};


        matTrain.historyL = savestate_historyL{subject};
        matTrain.historyR = savestate_historyR{subject};
    catch
        disp('Error saving a test file');
    end

end
