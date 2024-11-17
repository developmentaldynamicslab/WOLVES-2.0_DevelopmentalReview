accuracy=zeros(numSubjects,1);

for sj=1:numSubjects

    resp=test(sj).respLoc(:,1)-2;
    corrResp=test(sj).corrLoc(:,1);
    accuracy(sj,1)=mean(resp==corrResp);
end

name=[simName,'Prod_'];
save([name '_accuracy.mat'],'accuracy','-mat');


