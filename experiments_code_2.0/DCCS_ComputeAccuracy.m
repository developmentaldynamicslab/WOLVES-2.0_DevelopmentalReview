prepostTR=[3:8;
    10:15];
accuracy=zeros(size(test(:),1),2);

for sj=1:size(test(:),1)
    for phase=1:2
        accuracyTmp=[];
        for trI=prepostTR(phase,:)
            respLoc=(mean(find(test(sj).atnsa(trI,:)>.5))>50)+1;
            accuracyTmp=[accuracyTmp respLoc==test(sj).corrResp(trI)];
        end
        accuracy(sj,phase)=mean(accuracyTmp);
    end
end

mean(accuracy,1)

