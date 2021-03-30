%Function:          Create and refine profiles.
%inputs:            manual import of aaindexvalstrim2.txt
%outputs:           Physical profile of each amino acid on a residue by
%                   residue basis.
%other              Nil.

clear

sophies=input('Is this operation for BRES generation or to make profiles to be assessed? (1 for BRES, 2 for assessment)     ');
aastoprof=input('What is the filename containing the unaligned sequences for which you want me to generate physical profiles?     ', 's');
aaset=importdata(aastoprof);
aaset2=reshape(aaset', 1, numel(aaset));
aastr=strjoin(aaset2(1,:));
aaset3=strsplit(aastr,'>');
aaset3=aaset3';
aaset3(1,:)=[];

for profcrtindx=1:size(aaset3,1)
    tmpstr=char(aaset3(profcrtindx,1));
    tmpstrary=strsplit(tmpstr,' ');
    aaset3(profcrtindx,1:2)=tmpstrary(1,1:2);
end

clear profcrtindx

aaset3tmp=aaset3(:,2);

for aaset3tmpindx=1:numel(aaset3tmp)
    val(aaset3tmpindx)=numel(aaset3tmp{aaset3tmpindx});
end

vpout=aaset3tmp(val==max(val));
vout=vpout{1,1};
vout=char(vout);

clear aaset3tmp

aaset3a=num2cell(zeros(size(aaset3,1),length(vout)));

for profcrtindx=1:size(aaset3,1)
    cellContents = aaset3{profcrtindx,2};
    aaset3{profcrtindx,2}=cellContents(1:end);
    s={aaset3(profcrtindx,2)};
    sout=s{1,1};
    sout=char(sout);
    sout(sout==' ')=[];
    d={aaset3(profcrtindx,2)};
    dout=d{1,1};
    dout=char(dout);
    dout=strrep(dout, ' ', ''); 
    
    for splitdx = 1:length(dout)
        dout2(splitdx,1)=cellstr(dout(splitdx));
    end
    
    aaset3{profcrtindx,2}=sout;
    dout2=cellstr(dout2);
    doutdim=size(dout2,1);
    dout2=transpose(dout2);
    aaset3a(profcrtindx,1:size(dout2,2))=dout2;
    clear dout2 s sout
end
clear aaset aaset2 aaset3tmp aaset3tmpindx val
aasetcombprot=[aaset3, aaset3a];
clear aaset aaset2 aaset3 aaset3a

aaindexident1=importdata('aaindex1ident.txt');
aaindex1valstrim2=num2cell(importdata('aaindex1valstrim2.txt'));
aaindexvalsub=num2cell(zeros((size(aaindex1valstrim2,1)/2),20));
for subsindex=1:2:size(aaindex1valstrim2,1)
    aaindexvalsub((subsindex/2)+0.5,1:10)=aaindex1valstrim2(subsindex,1:10);
end
clear subsindex
for subsindex=2:2:size(aaindex1valstrim2,1)
    aaindexvalsub((subsindex/2),11:20)=aaindex1valstrim2(subsindex,1:10);
end
clear subsindex

for k=1:size(aaindexident1,1)
    cellContents=aaindexident1{k};
    aaindexident1{k} = cellContents(3:end);
end

aaindexident1(473:482,:)=[];
aaindexident1(514:515,:)=[];
aaindexvalsub(473:482,:)=[];
aaindexvalsub(514:515,:)=[];
aaindexvalsub(146,:)=[];
aaindexident1(146,:)=[];
aaindexvalsub(379,:)=[];
aaindexident1(379,:)=[];
aaindexvalsub(509,:)=[];
aaindexident1(509,:)=[];
aaindexvalsub(26:31,:)=[];
aaindexident1(26:31,:)=[];
aaindexvalsub(76,:)=[];
aaindexident1(76,:)=[];

profmat=num2cell(zeros((size(aaindexident1,1)+1),size(aasetcombprot,2)));
profmat(2:size(profmat,1),1)=aaindexident1(:,1);

for subsindex=1:size(aasetcombprot,1);
    profmat(1,:)=aasetcombprot(subsindex,:);
    for subsubsindex=3:size(profmat,2)
        switchvar=profmat{1,subsubsindex};
        switch switchvar;
            case 'A'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,1);
            case 'R'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,2);
            case 'N'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,3);
            case 'D'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,4);
            case 'C'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,5);
            case 'Q'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,6);
            case 'E'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,7);
            case 'G'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,8);
            case 'H'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,9);
            case 'I'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,10);
            case 'L'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,11);
            case 'K'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,12);
            case 'M'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,13);
            case 'F'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,14);
            case 'P'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,15);
            case 'S'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,16);
            case 'T'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,17);
            case 'W'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,18);
            case 'Y'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,19);
            case 'V'
                profmat(2:size(profmat,1),subsubsindex)=aaindexvalsub(:,20);
        end
    end
    
    textfilename2 = ['p' num2str(subsindex) '.csvx'];
    textfilename3 = ['p' num2str(subsindex) 'bp.csvx'];
    profmatwr=cell2mat(profmat(2:end,3:end));
    csvwrite(textfilename2,profmatwr);
   if sophies==1
    profmatwr( :, ~any(profmatwr,1) ) = [];
    profmatcap1=profmatwr(:,end-9:end);
    profmatcap2=profmatwr(:,1:10);
    profmatwr=[profmatcap1 profmatwr profmatcap2];
    csvwrite(textfilename3,profmatwr);
   end
    indextracker{subsindex,1}=textfilename2;
    indextracker(subsindex,2)=profmat(1,1);
    indextracker(subsindex,3)=profmat(1,2);
        
end

tabletrack=cell2table(indextracker(:,1:3), 'VariableNames', {'Filename' 'SeqID' 'FullSeq'});
writetable(tabletrack, 'indexguide.txt');

