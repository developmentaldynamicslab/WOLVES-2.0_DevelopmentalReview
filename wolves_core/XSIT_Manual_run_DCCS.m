close all; clear all;
tic

if 0 %choose for mac...
    run('..\COSIVINA\setpath.m') % add cosivina and jsoblab files to the matlab path
end
if 0 % choose for hpc at Tennessee...
    run('cosivina/setpath.m') % add cosivina and jsoblab files to the matlab path
    addpath('jsonlab');
    % addpath('../experiments_code');
    addpath('experiments_code_2.0');
    addpath('support_code');
    addpath('wolves_core');
    parpool('Processes',48);

end
if 0 % choose for hpc at UEA...
    run('../../../cosivina/setpath.m') % add cosivina and jsoblab files to the matlab path
    addpath('../../../jsonlab');
    addpath('../experiments_code');
    addpath('../experiments_code_2.0');
    addpath('../support_code');
    addpath('../wolves_core');
    parpool('SlurmProfile1',96)
end

mode=1; %1 = auto/gui, 0 = single batch; 2 for multicore batch mode (also need to switch to 'parfor' in the experiment file)
gui_speed=10; %10; %update gui after every n iterations: ideal values from 1 to 20.
notes = ['wolvesPaperPR.json'];% notes for experimenting simulations


%% loop over ages for DCCS tasks
for agect = 1:2
%for agect = 2 %uncomment this line (and comment line above) to run orig XSit tasks below

    if agect == 1
        youngDCCS=1;
        age='young';
        tasklist = [1 2 3 4 5 6];
    else
        youngDCCS=0;
        age='old';
        tasklist = [1 2 3];
    end

    for DCCS = tasklist
    %for DCCS = 0 %uncomment this line (and comment line above) to run orig XSit tasks below

        if DCCS==1 %original DCCS
            task=15;
        elseif DCCS==2 %Production
            task=16;
        elseif DCCS==3 %Comprehension
            task=17;
        elseif DCCS==4 %Munakata Uninformative
            task=18;
        elseif DCCS==5 %Munakata Novel Stimuli
            task=19;
        elseif DCCS==6 %Munakata PPC (Baseline sims)
            task=20;
        end

        createComboSim_DCCS; % create the model

        if sim.loadSettings('wolvesPaperPR.json','changeable') == 0; disp('json file load ERROR!!'); end % load the parameters file
        sim.init(); %initialize the model

        createComboGUI_paper; %GUI use in DevReview paper
        %createComboGUI; %original WOLVES GUI
        %createComboControls; %create and initialize GUI controls

        %% Update Memory Build and Decay Parameters
        parBuild=1200;
        parDecay=5500;
        setMemoryTraceTimescales(sim, parBuild, parDecay);

        %% Changes to parameters -- see table 1 DevReview paper
        sim.setElementParameters('wf1 -> word','sigma',0);
        sim.setElementParameters('wf2 -> word','sigma',0);
        sim.setElementParameters('wm_c1 -> wm_s','amplitude',0.6); %0.3 in PR
        sim.setElementParameters('wm_c2 -> wm_s','amplitude',0.6); %0.3 in PR
        sim.setElementParameters('atn_f1 -> wf1','amplitude',1.4); %1 in PR
        sim.setElementParameters('atn_f2 -> wf2','amplitude',1.4); %1 in PR
        sim.setElementParameters('hword -> word','amplitude',0.2); %0.1 in PR
        sim.setElementParameters('hwm_c1 -> wm_c1','amplitude',0.4); %1 in PR
        sim.setElementParameters('hwm_c2 -> wm_c2','amplitude',0.4); %1 in PR        
        sim.setElementParameters('wm_f1 -> wm_f1','amplitudeExc',22); %24 in PR
        sim.setElementParameters('wm_f2 -> wm_f2','amplitudeExc',22); %24 in PR
        sim.setElementParameters('wm_c1 -> wm_c1 (global)','amplitude',-.01); %0 in PR
        sim.setElementParameters('wm_c2 -> wm_c2 (global)','amplitude',-.01); %0 in PR
        sim.setElementParameters('wm_c1 -> atn_c1','amplitude',9.5); %5.5 in PR
        sim.setElementParameters('wm_c2 -> atn_c2','amplitude',9.5); %5.5 in PR

        %% new to WOLVES 2.0
        sim.setElementParameters('wf1 -> wm_c1','amplitude',1.1); %not in PR
        sim.setElementParameters('wf2 -> wm_c2','amplitude',1.1); %not in PR

        %% Parameters uniquely changed for DCCS, Production, Comprehension
        if DCCS > 0
            sim.setElementParameters('hwf1 -> wf1','amplitude',1); %4 in PR; 1 in DCCS
            sim.setElementParameters('hwf2 -> wf2','amplitude',1); %4 in PR; 1 in DCCS
            sim.setElementParameters('wf1 -> word','amplitude',0.2); %0.05 in PR; 0.2 in DCCS
            sim.setElementParameters('wf2 -> word','amplitude',0.2); %0.05 in PR; 0.2 in DCCS

            %model needs to attend more to the words in comprehension...like XSit
            if DCCS == 3
                sim.setElementParameters('word -> wf','amplitude',4.25); %4.25 in PR; 2.2 in DCCS
            else
                sim.setElementParameters('word -> wf','amplitude',2.2); %4.25 in PR; 2.2 in DCCS
            end

        else
            %% parameters changed back to original values for XSit
            sim.setElementParameters('hwf1 -> wf1','amplitude',4); %4 in PR; 1 in DCCS
            sim.setElementParameters('hwf2 -> wf2','amplitude',4); %4 in PR; 1 in DCCS
            sim.setElementParameters('wf1 -> word','amplitude',0.05); %0.05 in PR; 0.2 in DCCS
            sim.setElementParameters('wf2 -> word','amplitude',0.05); %0.05 in PR; 0.2 in DCCS
            sim.setElementParameters('word -> wf','amplitude',4.25); %4.25 in PR; 2.2 in DCCS
        end

        %% To run a batch of simulations on a HPC/ multicore pc, change variable mode = 2 (above).
        % running 96 on HPC to match cores available 
        if (mode == 2), numSubjects = 96; tolerance = 0; % specify the number of simulations/subjects to run. 300 default
        else,   numSubjects = 1; tolerance = 3; end % tolerance is used for gui's only to ensure they don't skip equality(==)conditionals

        %% developmental parameters for DCCS
        if DCCS > 0
            if youngDCCS
                wstrength = sort(pearsrnd(5,.15,3.6,14,numSubjects,1));%5.4;%4.9;%8;% word stimulus strength
                memTraceStrength=sort(pearsrnd(.375,.12,3.6,14,numSubjects,1));%.35;%.5;
                memTraceStrengthJbase=[.2, .1, .3, .45, .4, .1]; %.1
            else
                wstrength =sort(pearsrnd(5.4,.12,-3.6,14,numSubjects,1));% 5.4*ones(numSubjects,1);%8;% word stimulus strength
                memTraceStrength=sort(pearsrnd(.8,.12,-3.6,14,numSubjects,1));%.9*ones(numSubjects,1);%.1;
                memTraceStrengthJbase=[.1, .1, .1, .1, .1, .1]; %.1;
            end
        end        

        for tasknum = task
        %for tasknum = [1 3 8]  %uncomment this line (and comment line above to run orig XSit tasks 

            taskvar = tasknum; %choose the task/experiment (taskvar value) to simulate: default Smith & Yu (2008,11)

            %% create a name for your sims
            simNamePar = strcat('FullOrig2_', age, '_',num2str(tasknum),'_'); % give a name to your simulation.

            if (taskvar==1)
                %% Task - Smith & Yu, Dev Sci, 2008 - Yu & Smith, Dev Sci, 2011 - Infants - standard cross-sit
                simName = [simNamePar '_Smith_Yu_2008_'];
                Smith_Yu_2008_2011;
            elseif (taskvar==2)
                %% Task - Smith & Yu, Lang Learn Dev, 2013 - Infants - novelty trap
                simName = [simNamePar '_Smith_Yu_2013_'];
                Smith_Yu_2013;
            elseif (taskvar==3)
                %% Task - Vlach & Johnson, Cognition, 2013 - Infants - massed vs interleaved
                simName = [simNamePar '_Vlach_Johnson_2013_'];
                %Vlach_Johnson_2013;
                Vlach_CSWL_2013_17_19;
            elseif (taskvar==401)
                %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
                simName = [simNamePar '_Vlach_DeBrockOR_2017_'];
                %Vlach_DeBrockOR_2017;
                Vlach_DeBrock_WOB_2017;
            elseif (taskvar==402)
                %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
                simName = [simNamePar '_Vlach_DeBrockWR_2017_'];
                %Vlach_DeBrockWR_2017;
                Vlach_DeBrock_WOB_2017;
            elseif (taskvar==403)
                %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
                simName = [simNamePar '_Vlach_DeBrockWOB_2017_'];
                %Vlach_DeBrockWOB_2017;
                Vlach_DeBrock_WOB_2017;
            elseif (taskvar==5)
                %% Task - Fitneva & Christiansen, Cognitive Science (2015) - Partial Knowledge, Initial Accuracy
                simName = [simNamePar '_Fitneva_Christiansen_2015_'];
                Fitneva_Christiansen_2015; % use 300x2 conditions numSubjects
            elseif (taskvar==6)
                %% Task - Kachergis_Yu_Shiffrin, Psychon Bull Rev (2012) - Adults- Prior Learning, ME
                simName = [simNamePar '_Kachergis_Yu_Shiffrin_2012_'];
                Kachergis_Yu_Shiffrin_2012_Modified; % use 300x12 conditions numSubjects
            elseif (taskvar==7)
                %% Task - Yurovsky_Yu_Smith, Cognitive Science (2013) - Adults - Competitive Processes
                simName = [simNamePar '_Yurovsky_Yu_Smith_2013_'];
                Yurovsky_Yu_Smith_2013; % use 300x3 conditions numSubjects
            elseif (taskvar==8)
                %% Task - Suanda_Mugwanya_Namy, Jrl Exp Child Psy (2014) - 6 year olds - Referential Ambiguity
                simName = [simNamePar '_Suanda_Mugwanya_Namy_2014_'];
                Suanda_Mugwanya_Namy_2014; % use 300x3 conditions numSubjects
            elseif (taskvar==902)
                %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 2x2
                simName = [simNamePar '_Yu_Smith_Two_2007_'];
                Yu_Smith_Two_2007;
            elseif (taskvar==903)
                %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 3x3
                simName = [simNamePar '_Yu_Smith_Three_2007_'];
                Yu_Smith_Three_2007;
            elseif (taskvar==904)
                %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 4x4
                simName = [simNamePar '_Yu_Smith_Four_2007_'];
                Yu_Smith_Four_2007;
            elseif (taskvar==11)
                %% Yu Zhong & Fricker, Frontiers in Psychology(2012) - Adults - Eyetracking
                simName = [simNamePar '_Yu_Zhong_Fricker_2012_'];
                Yu_Zhong_Fricker_2012;
            elseif (taskvar==12)
                %% Trueswell Medina Hafri & Gleitman Cognitive Psychology(2013) - Adults - Eyetracking
                simName = [simNamePar '_Trueswell_Medina_Hafri_Gleitman_2013_'];
                Trueswell_Medina_Hafri_Gleitman_2013;
            elseif (taskvar==13)
                %% Mather Schafer Houston-Price (2011) - Infants Silent Condition
                simName = [simNamePar '_Mather_Schafer_Houston_Price_2011_'];
                Labelling_condition_ON = 0;%%change to 0 for Silent condition
                Mather_Schafer_Houston_Price;
            elseif (taskvar==14)
                %% Mather Schafer Houston-Price (2011) - Infants Labelling Condition
                simName = [simNamePar '_Mather_Schafer_Houston_Price_2011_'];
                Labelling_condition_ON = 1;%%change to 0 for Silent condition
                Mather_Schafer_Houston_Price;
            elseif (taskvar==15)
                %% Task - DCCS
                simName = [simNamePar 'Standard_'];
                DCCS_Basic;
            elseif (taskvar==16)
                %% Task - Label Production
                simName = [simNamePar 'Production_'];
                Production;
            elseif (taskvar==17)
                %% Task - Comprehension
                simName = [simNamePar 'Comprehension_'];
                Comprehension;
            elseif (taskvar==18)
                %% Task - Uninformative DCCS (no labels during pre-switch)
                simName = [simNamePar 'Uninfo_'];
                DCCS_Uninform;
            elseif (taskvar==19)
                %% Task - Uninformative DCCS (no labels during pre-switch)
                simName = [simNamePar 'NovelStimuli_'];
                DCCS_NovelStimuli;
            elseif (taskvar==20)
                simName = [simNamePar 'PPC_'];
                DCCS_PPC;
            end

            %% Save Simulated Results
            train=[];test=[];
            for subject=1:numSubjects
                try
                    if (taskvar ~= 13) & (taskvar ~= 14) &  (taskvar ~= 15) &  (taskvar ~= 16) &  (taskvar ~= 17) &  (taskvar ~= 18) &  (taskvar ~= 19) &  (taskvar ~= 20)
                        OutName1 = [simName,num2str(subject),'_train.mat'];
                        OutName2 = [simName,num2str(subject),'_test.mat'];
                        tempTrn=load(OutName1);
                        tempTst=load(OutName2);
                        [tempTrn(:).subject] = subject;
                        train = [train; tempTrn];
                        [tempTst(:).subject] = subject;
                        test = [test; tempTst];
                        delete(OutName1);
                        delete(OutName2);
                    else
                        OutName2 = [simName,num2str(subject),'_test.mat'];
                        tempTst=load(OutName2);
                        [tempTst(:).subject] = subject;
                        test = [test; tempTst];
                        delete(OutName2);
                    end
                catch
                    disp('Error on concatenating subject number ');
                    disp(subject);
                    continue;
                end
            end
            OutName = [simName,'results.mat'];
            save(OutName,'train','test','sim','notes');

        end %tasks

        toc

        if mode == 2 & DCCS > 0

            if DCCS==1 || DCCS==4 || DCCS==5 || DCCS==6

                DCCS_ComputeAccuracy;

                passpre=accuracy(:,1)>.6;
                failpre=accuracy(:,1)<.4;

                passpass=accuracy(passpre,2)>.6;
                passfail=accuracy(passpre,2)<.4;
                passmix=1-mean(passpass)-mean(passfail);

                sprintf('Passed Preswitch: %0.02f \nFailed Preswitch: %0.02f \n\nPassed Postswitch: %0.02f \nFailed Postswitch: %0.02f \nMixed Postswitch: %0.02f',mean(passpre),mean(failpre), mean(passpass), mean(passfail),passmix)

                outname=[simName,'DCCSResults.txt'];
                outfile=fopen(outname,'a+');
                if DCCS==1
                    taskWrite='Standard';
                elseif DCCS==4
                    taskWrite='Uninfo';
                elseif DCCS==5
                    taskWrite='NovelFeat';
                elseif DCCS==6
                    taskWrite='PPC';
                end

                if youngDCCS
                    fprintf(outfile,sprintf('Young Model %s\n',taskWrite));
                else
                    fprintf(outfile,sprintf('Old Model\n'))
                end

                date=clock;
                fprintf(outfile,sprintf('Run completed on %d/%d/%d at %d:%d \n',date(1),date(2),date(3),date(4),date(5)));
                fprintf(outfile,sprintf('Passed Preswitch: %0.02f \nFailed Preswitch: %0.02f \n\nPassed Postswitch: %0.02f \nFailed Postswitch: %0.02f \nMixed Postswitch: %0.02f\n\n\n',mean(passpre),mean(failpre), mean(passpass), mean(passfail),passmix));

            elseif DCCS==2

                Production_ComputeAccuracy;
                sprintf('Production Accuracy: %0.02f',mean(accuracy));
                outname=[simName,'ProductionResults.txt'];
                outfile=fopen(outname,'a+');
                date=clock;
                fprintf(outfile,sprintf('Run completed on %d/%d/%d at %d:%d \n',date(1),date(2),date(3),date(4),date(5)));
                if youngDCCS
                    fprintf(outfile,sprintf('Young Model\n'))
                else
                    fprintf(outfile,sprintf('Old Model\n'))
                end
                fprintf(outfile,sprintf('Production Accuracy: %0.02f\n\n',mean(accuracy)));

            elseif DCCS==3

                Comprehension_ComputeAccuracy;
                sprintf('Comprehension Accuracy: %0.02f',mean(mean(accuracy)))
                outname=[simName,'ComprehensionResults.txt'];
                outfile=fopen(outname,'a+');
                date=clock;
                fprintf(outfile,sprintf('Run completed on %d/%d/%d at %d:%d \n',date(1),date(2),date(3),date(4),date(5)));
                if youngDCCS
                    fprintf(outfile,sprintf('Young Model\n'))
                else
                    fprintf(outfile,sprintf('Old Model\n'))
                end
                fprintf(outfile,sprintf('Comprehension Accuracy: %0.02f\n\n',mean(mean(accuracy))));

            end

        end
    end
end