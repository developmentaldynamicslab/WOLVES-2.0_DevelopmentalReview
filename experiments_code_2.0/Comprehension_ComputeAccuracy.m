respOptions=[8 25 42 59 75 92];

for sj =1:numSubjects
    for tr=1:maxTR
        respLocTmp=mean(find(test(sj).atnsa(tr,:)>.5));
        respLocDiff=abs(respOptions-respLocTmp);
        respLoc=find(respLocDiff==min(respLocDiff));
        corrLoc=find(test(sj).objList(tr,:,1)==test(sj).wordList(tr));
        clStore(sj,tr)=corrLoc;

        if isempty(respLoc)
            rlStore(sj,tr)=-1;
            accuracy(sj,tr)=0;
        else
            rlStore(sj,tr)=respLoc;
            acc=respLoc==corrLoc;
            accuracy(sj,tr)=acc;
        end
    end

end

name=[simName,'Comp_'];
save([name '_accuracy.mat'],'accuracy','-mat');