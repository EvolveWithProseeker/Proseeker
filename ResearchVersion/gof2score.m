 function [letsgo]= gof2score(letsgo)

%Function:          Converts gofchecks to scores.

%inputs:            gofcheck files

%outputs:           Compacted scores for individual proteins based on
%                   functional motif detection

clear files files2
files=table2cell(readtable('1klistgof.txt','ReadVariableNames',false));
altseqs=transpose(importdata('altseqary.csv'));

for gofmodindx=1:1000
    gofrescore=importdata(files{gofmodindx,1});
    toDelete=(sum(gofrescore,1)==1275);
    gofrescore(:,toDelete)=[];
    tmpscrcnt=1;
    for scrmtchindx=1:12:size(gofrescore,2)
        tmpscrmat(tmpscrcnt,2)=sum(sum(gofrescore(:,scrmtchindx:scrmtchindx+3)));
        tmpscrmat(tmpscrcnt,3)=sum(sum(gofrescore(:,scrmtchindx+6:scrmtchindx+9)));
        tmpscrmat(tmpscrcnt,1)=sum(sum(tmpscrmat(tmpscrcnt,2:3)));
        tmpscrmat(tmpscrcnt,4:size(altseqs,2)+3)=altseqs(gofmodindx,1:size(altseqs,2));
        tmpscrcnt=tmpscrcnt+1;
    end
    if gofmodindx==1
        scrcnt=tmpscrmat;
    else
        scrcnt=[scrcnt;tmpscrmat];
    end
    simsearchres=sortrows(scrcnt,3);
end

dlmwrite('simsearchres',simsearchres);

