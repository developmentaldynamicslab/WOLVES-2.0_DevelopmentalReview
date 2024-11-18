sprintf('running standard DCCS version')
wmc_boost=0;

nObjects = 6; %%total number of objects
t_max = historyDuration; %specify simulation time  % scale with Experiment training trial Duration
wstrengthInstruct=15;
sptStrength=15;

colorPropertiesIndices=1:6;
colorLabelIndices=[3 4 5 6 7 8];
shapePropertiesIndices=9:14;
shapeLabelIndices=[11 12 13 14 15 16];

nLabels= 12;

sim.init();

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

targetList=[80 80;
    230 230];
testList=[80 230;
    230 80];

sr_sa=sim.getElement('atn_sr -> atn_sa').amplitude;
vis1_wms=sim.getElement('vis_f1 -> wm_s').amplitude;
vis2_wms=sim.getElement('vis_f2 -> wm_s').amplitude;
vis1_cons=sim.getElement('vis_f1 -> con_s').amplitude;
vis2_cons=sim.getElement('vis_f2 -> con_s').amplitude;
ior_atn=sim.getElement('ior_s -> atn_sa').amplitude;

maxTRall = 15; %total trials over phases

%% parfor loop -- use parfor for mode 2, use for for mode 1
for subject = 1:numSubjects 

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

    %subject
    trnum_cont = 0; %counter of trials across phases...

    %% Reset stored data
    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyL{subject}=zeros(maxTRall,t_max);
    savestate_historyR{subject}=zeros(maxTRall,t_max);
    savestate_atnsa{subject}=zeros(maxTRall,fieldSize_spt);
    corrRespTmp{subject}=zeros(maxTRall,1);

    %%Reset Hebbian Fields
    sim2.getElement('hwm_c1').output=zeros(fieldSize_ftr,fieldSize_spt);
    sim2.getElement('hwm_c2').output=zeros(fieldSize_ftr,fieldSize_spt);
    sim2.getElement('hwm_s').output=zeros(1,fieldSize_spt);
    sim2.getElement('hword').output=zeros(1,fieldSize_wd);
    sim2.getElement('hwf1').output=zeros(fieldSize_ftr,fieldSize_wd);
    sim2.getElement('hwf2').output=zeros(fieldSize_ftr,fieldSize_wd);
    sim2.getElement('hcon_s').output=zeros(1,fieldSize_spt);
    sim2.getElement('hwm_f1').output=zeros(1,fieldSize_ftr);
    sim2.getElement('hwm_f2').output=zeros(1,fieldSize_ftr);
    sim2.getElement('hcon_f1').output=zeros(1,fieldSize_ftr);
    sim2.getElement('hcon_f2').output=zeros(1,fieldSize_ftr);

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

    %% Re-initialize toggled parameters
    sim2.setElementParameters('atn_sr -> atn_sa','amplitude',sr_sa);
    sim2.setElementParameters('vis_f1 -> wm_s','amplitude',vis1_wms);
    sim2.setElementParameters('vis_f2 -> wm_s','amplitude',vis2_wms);
    sim2.setElementParameters('vis_f1 -> con_s','amplitude',vis1_cons);
    sim2.setElementParameters('vis_f2 -> con_s','amplitude',vis2_cons);
    sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

    %% DCCS_Instructions

    visuals_On = floor(250/scale_factor);    visuals_Off = t_max - floor(500/scale_factor); % visual display on off timing

    instruct1_On=floor(500/scale_factor);    instruct1_Off = floor(3250/scale_factor);
    instruct2_On=floor(3750/scale_factor);    instruct2_Off = floor(6500/scale_factor);

    % DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[];
    Words=[];

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
    handle_atn_sa{subject}=sim2.getElement('atn_sa');

    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyL{subject}=zeros(maxTRall,t_max);
    savestate_historyR{subject}=zeros(maxTRall,t_max);
    training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    savestate_atnsa{subject}=zeros(maxTRall,fieldSize_spt);

    trnum_cont = trnum_cont+1;
    %[subject trnum_cont]

    sim2.setElementParameters({'wf1' 'wf2' 'word'},{'h' 'h' 'h'},[hStore hStore hStoreW]);
    sim2.t =sim2.tZero;

    for wi=1:6
        sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
            {(colorLabelIndices(wi)),wstrength(subject)});
    end

    while sim2.t <= t_max % run the ime loop of a trial

        t = sim2.t;

        if isalmost(t, visuals_On,tolerance) %present visual stimuli
            for f=1:nFeatures
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(1,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Lm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));



                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(2,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Rm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));

            end
        end

        if isalmost(t, instruct1_On,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacel', ... % set feature dimension value
                },{'position'}, {spt_Lm}); %Make list of features to present
            sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

        elseif isalmost(t, instruct1_On+100,tolerance) %present visual stimuli
            sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},...
                {(colorLabelIndices(2)),wstrengthInstruct});

        elseif isalmost(t, instruct1_Off,tolerance) %present visual stimuli
            for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),0});
            end
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field
            sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);
        end

        if isalmost(t, instruct2_On,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacel', ... % set feature dimension value
                },{'position'}, {spt_Rm}); %Make list of features to present
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

        elseif isalmost(t, instruct2_On+100,tolerance) %present visual stimuli
            for wi=1:6%nLabels \
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),wstrength(subject)});
            end

            sim2.setElementParameters({'Word_Label_5'}, {'position','amplitude'},...
                {(colorLabelIndices(5)),wstrengthInstruct});

        elseif isalmost(t, instruct2_Off,tolerance) %present visual stimuli
            for wi=1:6%nLabels
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),0});
            end
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);
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
            savestate_historyL{subject}(trnum_cont,:) = flipud(sim2.getComponent('history lookL', 'output'));
            savestate_historyR{subject}(trnum_cont,:) =flipud(sim2.getComponent('history lookR', 'output'));
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
    end % Time sim.t

    %% Initialize

    for f=1:nFeatures
        sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
            }, {'amplitude'}, 0 * ones(1, 1));
        sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
            }, {'amplitude'}, 0 * ones(1, 1));
    end

    %% DCCS_Demo

    demo1_On=floor(250/scale_factor);    demo1_Off = floor(1250/scale_factor);
    demo2_On=floor(1500/scale_factor);    demo2_Off = floor(2500/scale_factor);

    targets_On=floor(2750/scale_factor);    targets_Off=floor(6750/scale_factor);

    trnum_cont = trnum_cont+1;

    sim2.t =sim2.tZero;

    while sim2.t <= t_max % run the ime loop of a trial

        t = sim2.t;

        if isalmost(t, demo1_On,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field
            sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                }, {'position'}, {spt_Lm});   %Make value for center of field

            for f=1:nFeatures
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                    },{'positionY'}, {testList(1,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Lm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));
            end

            for wi=1:6%nLabels
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),wstrength(subject)});
            end
        elseif isalmost(t, demo1_Off,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field

            for f=1:nFeatures
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, 0 * ones(1, 1));
                sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                    }, {'amplitude'}, {0});   %Make value for center of field

            end

            for wi=1:6%nLabels
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),0});
            end
        end


        if isalmost(t, demo2_On,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacer', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field
            sim2.setElementParameters({'inputSpacer', ... % set spatial dimension value
                }, {'position'}, {spt_Rm});   %Make value for center of field
            for f=1:nFeatures
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                    },{'positionY'}, {testList(2,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Rm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));
            end
            for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),wstrength(subject)});
            end
        elseif isalmost(t, demo2_Off,tolerance) %present visual stimuli
            sim2.setElementParameters({'inputSpacer', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field

            for f=1:nFeatures
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, 0 * ones(1, 1));
                sim2.setElementParameters({'inputSpacer', ... % set spatial dimension value
                    }, {'amplitude'}, {0});   %Make value for center of field

            end
            for wi=1:6%nLabels
                sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                    {(colorLabelIndices(wi)),0});
            end
        end

        if isalmost(t, targets_On,tolerance)
            for f=1:2
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(1,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Lm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));


                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(2,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Rm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));
            end
            sim2.setElementParameters('atn_sr -> atn_sa','amplitude',sr_sa);
            sim2.setElementParameters('vis_f1 -> wm_s','amplitude',vis1_wms);
            sim2.setElementParameters('vis_f2 -> wm_s','amplitude',vis2_wms);
            sim2.setElementParameters('vis_f1 -> con_s','amplitude',vis1_cons);
            sim2.setElementParameters('vis_f2 -> con_s','amplitude',vis2_cons);
            sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

        elseif isalmost(t, targets_Off,tolerance)
            for f=1:2
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, 0);
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, 0);
            end
            sim2.setElementParameters('atn_sr -> atn_sa','amplitude',0);
            sim2.setElementParameters('vis_f1 -> wm_s','amplitude',0);
            sim2.setElementParameters('vis_f2 -> wm_s','amplitude',0);
            sim2.setElementParameters('vis_f1 -> con_s','amplitude',0);
            sim2.setElementParameters('vis_f2 -> con_s','amplitude',0);
            sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);
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
            savestate_historyL{subject}(trnum_cont,:) = flipud(sim2.getComponent('history lookL', 'output'));
            savestate_historyR{subject}(trnum_cont,:) =flipud(sim2.getComponent('history lookR', 'output'));
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

    end


    %% DCCS_Preswitch

    test_On=floor(500/scale_factor);    test_Off = floor(5500/scale_factor);
    targets_On=floor(5750/scale_factor);    targets_Off=floor(6650/scale_factor);

    maxTR=6;
    trList=[1 1 1 2 2 2];
    trList=trList(randperm(6));
    hStoreS=sim2.getElement('wm_s').h;

    sim2.setElementParameters('atn_sr -> atn_sa','amplitude',0);
    sim2.setElementParameters('vis_f1 -> wm_s','amplitude',0);
    sim2.setElementParameters('vis_f2 -> wm_s','amplitude',0);
    sim2.setElementParameters('vis_f1 -> con_s','amplitude',0);
    sim2.setElementParameters('vis_f2 -> con_s','amplitude',0);
    sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

    wf_h_init=sim2.getElement('wf1').h;
    word_h_init=sim2.getElement('word').h;
    wmc1_h_init=sim2.getElement('wm_c1').h;
    wmc2_h_init=sim2.getElement('wm_c2').h;

    for tr=1:maxTR %n trials loop

        trnum_cont = trnum_cont+1;
        %[subject trnum_cont]
        corrRespTmp{subject}(trnum_cont)=trList(tr);

        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hStore hStore]);
        sim2.setElementParameters('wm_s','h',hStoreS);
        sim2.t =sim2.tZero;

        test_Off_sorted = 0;

        while sim2.t <= t_max % run the time loop of a trial

            t = sim2.t;

            if t>test_On && t<test_Off
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
                    if hnew > wf_h_init
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hnew hnew]);
                    else
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[wf_h_init wf_h_init]);
                    end
                    sim2.setElementParameters('word','h',word_h_init);
                end
                if max(sim2.getElement('wm_s').output)<.9
                    %
                    sim2.setElementParameters('wm_c1','h',wmc1_h_init+wmc_boost);
                    sim2.setElementParameters('wm_c2','h',wmc2_h_init+wmc_boost);

                elseif max(sim2.getElement('atn_sa').activation)>8
                    savestate_atnsa{subject}(trnum_cont,:)  = sim2.getComponent('atn_sa', 'output');
                    test_Off_sorted = 1;
                    sim2.setElementParameters('wm_c1','h',wmc1_h_init);
                    sim2.setElementParameters('wm_c2','h',wmc2_h_init);
                end
            end


            if isalmost(t, test_On,tolerance) %present visual stimuli

                for f=1:nFeatures
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                        },{'positionY'}, {testList(trList(tr),f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_C});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 1.25*vstrength * ones(1, 1));
                end
                sim2.setElementParameters('boost wm_c','amplitude',2.1);


                for wi=1:6%nLabels
                    sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                        {(colorLabelIndices(wi)),wstrength(subject)});
                end


            elseif isalmost(t, test_Off,tolerance) || test_Off_sorted %present visual stimuli

                for f=1:nFeatures
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0 * ones(1, 1));
                    sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                        }, {'amplitude'}, {0});   %Make value for center of field

                end
                sim2.setElementParameters('boost wm_c','amplitude',0);

                for wi=1:6%nLabels
                    sim2.setElementParameters({['Word_Label_' num2str(wi)]}, {'position','amplitude'},...
                        {(colorLabelIndices(wi)),0});
                end

                sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[wf_h_init wf_h_init]);
                test_Off_sorted = 0;

            end
            if isalmost(t, targets_On,tolerance)
                for f=1:2
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                        },{'positionY'}, {targetList(1,f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_Lm});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, vstrength * ones(1, 1));


                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                        },{'positionY'}, {targetList(2,f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_Rm});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                        }, {'amplitude'}, vstrength * ones(1, 1));
                end
                sim2.setElementParameters('atn_sr -> atn_sa','amplitude',sr_sa);
                sim2.setElementParameters('vis_f1 -> wm_s','amplitude',vis1_wms);
                sim2.setElementParameters('vis_f2 -> wm_s','amplitude',vis2_wms);
                sim2.setElementParameters('vis_f1 -> con_s','amplitude',vis1_cons);
                sim2.setElementParameters('vis_f2 -> con_s','amplitude',vis2_cons);
                sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

            elseif isalmost(t, targets_Off,tolerance)
                for f=1:2
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0);
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0);
                end
                sim2.setElementParameters('atn_sr -> atn_sa','amplitude',0);
                sim2.setElementParameters('vis_f1 -> wm_s','amplitude',0);
                sim2.setElementParameters('vis_f2 -> wm_s','amplitude',0);
                sim2.setElementParameters('vis_f1 -> con_s','amplitude',0);
                sim2.setElementParameters('vis_f2 -> con_s','amplitude',0);
                sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);
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
                savestate_historyL{subject}(trnum_cont,:) = flipud(sim2.getComponent('history lookL', 'output'));
                savestate_historyR{subject}(trnum_cont,:) =flipud(sim2.getComponent('history lookR', 'output'));
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

        end
    end

    %% Reset
    sim2.setElementParameters('atn_sr -> atn_sa','amplitude',sr_sa);
    sim2.setElementParameters('vis_f1 -> wm_s','amplitude',vis1_wms);
    sim2.setElementParameters('vis_f2 -> wm_s','amplitude',vis2_wms);
    sim2.setElementParameters('vis_f1 -> con_s','amplitude',vis1_cons);
    sim2.setElementParameters('vis_f2 -> con_s','amplitude',vis2_cons);
    sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

    %% DCCS_Instructions_Post

    visuals_On = floor(250/scale_factor);    visuals_Off = t_max - floor(500/scale_factor); % visual display on off timing
    %word1_On  = floor(50/scale_factor);     word1_Off = floor(3000/scale_factor); % word 1 on off timing

    instruct1_On=floor(500/scale_factor);    instruct1_Off = floor(3250/scale_factor);
    instruct2_On=floor(3750/scale_factor);    instruct2_Off = floor(6500/scale_factor);

    % DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[];
    Words=[];

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
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%

    % RUN TRIAL

    trnum_cont = trnum_cont+1;
    %[subject trnum_cont]

    sim2.setElementParameters({'wf1' 'wf2' 'word'},{'h' 'h' 'h'},[hStore hStore hStoreW]);
    sim2.t =sim2.tZero;

    for wi=1:6%nLabels
        sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
            {(shapeLabelIndices(wi)),wstrength(subject)});
    end
    while sim2.t <= t_max % run the ime loop of a trial

        t = sim2.t;

        if isalmost(t, visuals_On,tolerance) %present visual stimuli
            for f=1:nFeatures

                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(1,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Lm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));

                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                    },{'positionY'}, {targetList(2,f)}); %Make list of features to present
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                    }, {'positionX'}, {spt_Rm});   %Make value for center of field
                sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                    }, {'amplitude'}, vstrength * ones(1, 1));

            end
        end

        if isalmost(t, instruct1_On,tolerance) %present visual stimuli

            sim2.setElementParameters({'inputSpacel', ... % set feature dimension value
                },{'position'}, {spt_Lm}); %Make list of features to present
            sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

        elseif isalmost(t, instruct1_On+100,tolerance) %present visual stimuli
            sim2.setElementParameters({'Word_Label_8'}, {'position','amplitude'},...
                {(shapeLabelIndices(2)),wstrengthInstruct});

        elseif isalmost(t, instruct1_Off,tolerance) %present visual stimuli
            for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
                    {(shapeLabelIndices(wi)),0});
            end
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

        end


        if isalmost(t, instruct2_On,tolerance) %present visual stimuli

            sim2.setElementParameters({'inputSpacel', ... % set feature dimension value
                },{'position'}, {spt_Rm}); %Make list of features to present
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {sptStrength});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

        elseif isalmost(t, instruct2_On+100,tolerance) %present visual stimuli
            for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
                    {(shapeLabelIndices(wi)),wstrength(subject)});
            end

            sim2.setElementParameters({'Word_Label_11'}, {'position','amplitude'},...
                {(shapeLabelIndices(5)),wstrengthInstruct});

        elseif isalmost(t, instruct2_Off,tolerance) %present visual stimuli
            for wi=1:6%nLabels %%NEED TO MAKE nLabels PARAMETERS
                sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
                    {(shapeLabelIndices(wi)),0});
            end
            sim2.setElementParameters({ 'inputSpacel', ... % set spatial dimension value
                }, {'amplitude'}, {0});   %Make value for center of field

            sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

        end

        if isalmost(t, t_max,tolerance) % save looking data and memory activity states
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
            savestate_historyL{subject}(trnum_cont,:) = flipud(sim2.getComponent('history lookL', 'output'));
            savestate_historyR{subject}(trnum_cont,:) =flipud(sim2.getComponent('history lookR', 'output'));
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

    end % Time sim.t


    %% Initialize
    for f=1:nFeatures
        sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
            }, {'amplitude'}, 0 * ones(1, 1));
        sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
            }, {'amplitude'}, 0 * ones(1, 1));
    end

    %% DCCS_Postswitch

    test_On=floor(500/scale_factor);    test_Off = floor(5500/scale_factor);
    targets_On=floor(5750/scale_factor);    targets_Off=floor(6650/scale_factor);

    maxTR=6;
    trList=[1 1 1 2 2 2];
    trList=trList(randperm(6));
    hStoreS=sim2.getElement('wm_s').h;

    sim2.setElementParameters('atn_sr -> atn_sa','amplitude',0);
    sim2.setElementParameters('vis_f1 -> wm_s','amplitude',0);
    sim2.setElementParameters('vis_f2 -> wm_s','amplitude',0);
    sim2.setElementParameters('vis_f1 -> con_s','amplitude',0);
    sim2.setElementParameters('vis_f2 -> con_s','amplitude',0);
    sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);

    wf_h_init=sim2.getElement('wf1').h;
    word_h_init=sim2.getElement('word').h;
    wmc1_h_init=sim2.getElement('wm_c1').h;
    wmc2_h_init=sim2.getElement('wm_c2').h;

    for tr=1:maxTR %n trials loop

        trnum_cont = trnum_cont+1;
        %[subject trnum_cont]
        corrRespTmp{subject}(trnum_cont)=abs(trList(tr)-3); %inverting so 1=2 and 2=1

        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hStore hStore]);
        sim2.setElementParameters('wm_s','h',hStoreS);
        sim2.t =sim2.tZero;

        test_Off_sorted = 0;

        while sim2.t <= t_max % run the time loop of a trial

            t = sim2.t;

            if t>test_On && t<test_Off
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
                    if hnew > wf_h_init
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[hnew hnew]);
                    else
                        sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[wf_h_init wf_h_init]);
                    end
                    sim2.setElementParameters('word','h',word_h_init);
                end
                if max(sim2.getElement('wm_s').output)<.9
                    %
                    sim2.setElementParameters('wm_c1','h',wmc1_h_init+wmc_boost);
                    sim2.setElementParameters('wm_c2','h',wmc2_h_init+wmc_boost);
                elseif max(sim2.getElement('atn_sa').activation)>8
                    savestate_atnsa{subject}(trnum_cont,:)  = sim2.getComponent('atn_sa', 'output');
                    test_Off_sorted = 1;
                    sim2.setElementParameters('wm_c1','h',wmc1_h_init);
                    sim2.setElementParameters('wm_c2','h',wmc2_h_init);
                end
            end

            if isalmost(t, test_On,tolerance) %present visual stimuli

                for f=1:nFeatures
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                        },{'positionY'}, {testList(trList(tr),f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_C});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 1.25*vstrength * ones(1, 1));
                end
                sim2.setElementParameters('boost wm_c','amplitude',2.1);


                for wi=1:6%nLabels
                    sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
                        {(shapeLabelIndices(wi)),wstrength(subject)});
                end

            elseif isalmost(t, test_Off,tolerance) || test_Off_sorted %present visual stimuli

                for f=1:nFeatures
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0 * ones(1, 1));
                    sim2.setElementParameters({'inputSpacel', ... % set spatial dimension value
                        }, {'amplitude'}, {0});   %Make value for center of field

                end
                sim2.setElementParameters('boost wm_c','amplitude',0);

                for wi=1:6%nLabels
                    sim2.setElementParameters({['Word_Label_' num2str(wi+6)]}, {'position','amplitude'},...
                        {(shapeLabelIndices(wi)),0});
                end

                sim2.setElementParameters({'wf1' 'wf2'},{'h' 'h'},[wf_h_init wf_h_init]);
                test_Off_sorted = 0;

            end

            if isalmost(t, targets_On,tolerance)
                for f=1:2
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set feature dimension value
                        },{'positionY'}, {targetList(1,f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_Lm});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, vstrength * ones(1, 1));


                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set feature dimension value
                        },{'positionY'}, {targetList(2,f)}); %Make list of features to present
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set spatial dimension value
                        }, {'positionX'}, {spt_Rm});   %Make value for center of field
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                        }, {'amplitude'}, vstrength * ones(1, 1));
                end
                sim2.setElementParameters('atn_sr -> atn_sa','amplitude',sr_sa);
                sim2.setElementParameters('vis_f1 -> wm_s','amplitude',vis1_wms);
                sim2.setElementParameters('vis_f2 -> wm_s','amplitude',vis2_wms);
                sim2.setElementParameters('vis_f1 -> con_s','amplitude',vis1_cons);
                sim2.setElementParameters('vis_f2 -> con_s','amplitude',vis2_cons);
                sim2.setElementParameters('ior_s -> atn_sa','amplitude',ior_atn);

            elseif isalmost(t, targets_Off,tolerance)
                for f=1:2
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_LeftMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0);
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_RightMid'), ... % set left right stimuli up
                        }, {'amplitude'}, 0);
                end
                sim2.setElementParameters('atn_sr -> atn_sa','amplitude',0);
                sim2.setElementParameters('vis_f1 -> wm_s','amplitude',0);
                sim2.setElementParameters('vis_f2 -> wm_s','amplitude',0);
                sim2.setElementParameters('vis_f1 -> con_s','amplitude',0);
                sim2.setElementParameters('vis_f2 -> con_s','amplitude',0);
                sim2.setElementParameters('ior_s -> atn_sa','amplitude',0);
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
                savestate_historyL{subject}(trnum_cont,:) = flipud(sim2.getComponent('history lookL', 'output'));
                savestate_historyR{subject}(trnum_cont,:) =flipud(sim2.getComponent('history lookR', 'output'));
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

        end
    end

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

        matTrain.atnsa = savestate_atnsa{subject};
        matTrain.corrResp = corrRespTmp{subject};

        matTrain.historyL = savestate_historyL{subject};
        matTrain.historyR = savestate_historyR{subject};
    catch
        disp('Error saving a test file');
    end

end


