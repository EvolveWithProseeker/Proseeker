 function proseeker

%Version:           v 9.1.5 RESEARCH

%Env:               MIN: 10 core 4.0 gHz, 1.5 gb RAM per core

%Function:          Iteratively evolve proteins based upon cumulative
%                   probability of mutation with residue substitutions based on quantitative
%                   amino acid descriptors.

%inputs:            FASTA formatted protein sequence for start.
%                   Blocking file - Comma seperated text file with 1 for a
%                   residue to mutate and 0 for a residue not to mutate.
%                   Easy to make, requires no skill BUT it MUST be exactly
%                   the same length as the start sequence AND is oriented
%                   in precisely the same way.
%                   Profiles of both key sites (this was made with sequences
%                   having two binding sites in mind).

%outputs:           Sequences.txt contains a complete list of sequence
%                   iterations.

%other:             This script is very sensitive to mutation rate. It is extremely important so it is
%                   suggested that it be carefully considered prior to
%                   acceptance of sequence output. Some sequence output may
%                   well represent what can be considered in real terms to
%                   be an evolutionary dead end so functionality of all
%                   sequences is far from garaunteed. It is also impossible
%                   for this script to add amino acids, mutate to stop
%                   codons or delete AAs. Essentially it just substitutes
%                   natural amino acids (if you want non-natural amino acids
%                   that will mean a substantial modification so look out for 
%                   that in later versions).

mkdir evolvescrs
warning('off','MATLAB:DELETE:FileNotFound');
warning('off','all');
poolobj = gcp('nocreate');
delete(poolobj);
gencounter=1; 
delete *.csvx *.csvres simsearchres simsearchresstr runstats.csv altseqary.csv *alts.txt
mkdir evolvescrs
rmdir ('evolvescrs', 's'); 
detfile='jobdet.strt.csv';
zeropoint='DHR8.txt';
BLOCK='BLOCK.txt';

jobstart=importdata(detfile);

generations=20;
generations=generations-1;
mkdir evolvescrs
genwidth=100;

blokyn=1; 

STRTTB = importdata(zeropoint,'\n');
STRTTB([1,1],:) = [];
Tab2Conv = char(STRTTB);
STRTTB2=cellstr(reshape(Tab2Conv,1,[])');
STRTTB2a=char(STRTTB2);
STRTTB2b=uint8(STRTTB2a);
wrkngcolmod=zeros(size(STRTTB2b,1),8);  %   This is used later. It is not immediately apparent but it is VITAL to the implementation of selection pressure.
profrefmata=[65,82,78,68,67,81,69,71,72,73,76,75,77,70,80,83,84,87,89,86;(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100),(100)];
profrefmat=transpose(profrefmata);

STRTTB2b=255.*(im2double(STRTTB2b));

s=double(STRTTB2b);
s(s==65)=100;
s(s==82)=100;
s(s==78)=100;
s(s==68)=100;
s(s==67)=100;
s(s==81)=100;
s(s==69)=100;
s(s==71)=100;
s(s==72)=100;
s(s==73)=100;
s(s==76)=100;
s(s==75)=100;
s(s==77)=100; %M - 90 to fix abundance issue
s(s==70)=100;
s(s==80)=100;
s(s==83)=100;
s(s==84)=100;
s(s==87)=100; %W - 90 to fix abundance issue
s(s==89)=100;
s(s==86)=100;
clcary=zeros(size(STRTTB2b,1),2);
clcary(:,1)=STRTTB2b;
clcary(:,2)=s;

sumdifmata=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,];
sumdifmata=[sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata;sumdifmata];
sumdifmat=zeros(21,21);
sumdifmat(2:21,1)=[65;82;78;68;67;81;69;71;72;73;76;75;77;70;80;83;84;87;89;86]; %  This line and the one below it might look redundant but they are not. Do not change them.
sumdifmat(1,2:21)=[65,82,78,68,67,81,69,71,72,73,76,75,77,70,80,83,84,87,89,86];
sumdifmat(2:21,2:21)=sumdifmata;

freqchoice=1;

if blokyn==1
    blockmat=transpose(importdata(BLOCK));
end
prodcol=(generations*2)+2;
prodary=zeros(size(clcary,1),prodcol);
prodary(:,1:2)=clcary;
if blokyn==1
    for blkindx=1:size(prodary,1)
        if blockmat(blkindx,1)==0
            prodary(blkindx,2)=0;
        end
    end
end

if freqchoice==1
    freq=0.05; % Mutation rate.
    realfreqa=size(clcary,1)*freq;
    realfreq=floor(realfreqa);
else
    freq=0.025;
    realfreqa=size(clcary,1)*freq;
    realfreq=floor(realfreqa);
    
end
disp(realfreq);
clear realfreqa
f = zeros(size(prodary,1),1);

for indxf=1:size(prodary,1)
    f(indxf,1)=indxf;
end

g = zeros(size(sumdifmat,1),1);
xca=[65;82;78;68;67;81;69;71;72;73;76;75;77;70;80;83;84;87;89;86];

for indxe=3:2:prodcol-1
    wrkngcol=prodary(:,indxe-2);
    wrkbp=wrkngcol;

    clear C indxg X xcmod xc
    
    if gencounter==1
        for wdindx=1:genwidth
            if wdindx==1
                C = cumsum(prodary(:,indxe-1));
                prodary(:,indxe+1)=prodary(:,indxe-1);
                
                for indxg=1:realfreq
                    X=f(1+sum(C(end)*rand>C));
                    if wrkngcol(X,1)==65
                        g(:,1)=sumdifmat(:,2);
                    elseif wrkngcol(X,1)==82
                        g(:,1)=sumdifmat(:,3);
                    elseif wrkngcol(X,1)==78
                        g(:,1)=sumdifmat(:,4);
                    elseif wrkngcol(X,1)==68
                        g(:,1)=sumdifmat(:,5);
                    elseif wrkngcol(X,1)==67
                        g(:,1)=sumdifmat(:,6);
                    elseif wrkngcol(X,1)==81
                        g(:,1)=sumdifmat(:,7);
                    elseif wrkngcol(X,1)==69
                        g(:,1)=sumdifmat(:,8);
                    elseif wrkngcol(X,1)==71
                        g(:,1)=sumdifmat(:,9);
                    elseif wrkngcol(X,1)==72
                        g(:,1)=sumdifmat(:,10);
                    elseif wrkngcol(X,1)==73
                        g(:,1)=sumdifmat(:,11);
                    elseif wrkngcol(X,1)==76
                        g(:,1)=sumdifmat(:,12);
                    elseif wrkngcol(X,1)==75
                        g(:,1)=sumdifmat(:,13);
                    elseif wrkngcol(X,1)==77
                        g(:,1)=sumdifmat(:,14);
                    elseif wrkngcol(X,1)==70
                        g(:,1)=sumdifmat(:,15);
                    elseif wrkngcol(X,1)==80
                        g(:,1)=sumdifmat(:,16);
                    elseif wrkngcol(X,1)==83
                        g(:,1)=sumdifmat(:,17);
                    elseif wrkngcol(X,1)==84
                        g(:,1)=sumdifmat(:,18);
                    elseif wrkngcol(X,1)==87
                        g(:,1)=sumdifmat(:,19);
                    elseif wrkngcol(X,1)==89
                        g(:,1)=sumdifmat(:,20);
                    elseif wrkngcol(X,1)==86
                        g(:,1)=sumdifmat(:,21);
                    end

                    sdevg=std(g(2:21,1));
                    gprep=g(2:21);
                    ming=min(gprep(gprep>0));
                    g(g>(ming+(sdevg)))=0;
                    gc = cumsum(g(:,1));
                    xcmod=1+sum(gc(end)*rand>gc);
                    while xcmod>size(xca,1)
                        xcmod=1+sum(gc(end)*rand>gc);
                    end
                    
                    xc=xca(xcmod);

                    if X==size(prodary,1)
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)= profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    elseif X==1
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                    else
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    end
                    clear X
                end

                scrseqary(:,wdindx)=wrkngcol;
                scrprbary(:,wdindx)=wrkngcol;
            else
                wrkngcol=wrkbp;
                C = cumsum(prodary(:,indxe-1));
                prodary(:,indxe+1)=prodary(:,indxe-1);
                
                for indxg=1:realfreq
                    X=f(1+sum(C(end)*rand>C));
                    if wrkngcol(X,1)==65
                        g(:,1)=sumdifmat(:,2);
                    elseif wrkngcol(X,1)==82
                        g(:,1)=sumdifmat(:,3);
                    elseif wrkngcol(X,1)==78
                        g(:,1)=sumdifmat(:,4);
                    elseif wrkngcol(X,1)==68
                        g(:,1)=sumdifmat(:,5);
                    elseif wrkngcol(X,1)==67
                        g(:,1)=sumdifmat(:,6);
                    elseif wrkngcol(X,1)==81
                        g(:,1)=sumdifmat(:,7);
                    elseif wrkngcol(X,1)==69
                        g(:,1)=sumdifmat(:,8);
                    elseif wrkngcol(X,1)==71
                        g(:,1)=sumdifmat(:,9);
                    elseif wrkngcol(X,1)==72
                        g(:,1)=sumdifmat(:,10);
                    elseif wrkngcol(X,1)==73
                        g(:,1)=sumdifmat(:,11);
                    elseif wrkngcol(X,1)==76
                        g(:,1)=sumdifmat(:,12);
                    elseif wrkngcol(X,1)==75
                        g(:,1)=sumdifmat(:,13);
                    elseif wrkngcol(X,1)==77
                        g(:,1)=sumdifmat(:,14);
                    elseif wrkngcol(X,1)==70
                        g(:,1)=sumdifmat(:,15);
                    elseif wrkngcol(X,1)==80
                        g(:,1)=sumdifmat(:,16);
                    elseif wrkngcol(X,1)==83
                        g(:,1)=sumdifmat(:,17);
                    elseif wrkngcol(X,1)==84
                        g(:,1)=sumdifmat(:,18);
                    elseif wrkngcol(X,1)==87
                        g(:,1)=sumdifmat(:,19);
                    elseif wrkngcol(X,1)==89
                        g(:,1)=sumdifmat(:,20);
                    elseif wrkngcol(X,1)==86
                        g(:,1)=sumdifmat(:,21);
                    end

                    sdevg=std(g(2:21,1));
                    gprep=g(2:21);
                    ming=min(gprep(gprep>0));
                    g(g>(ming+(sdevg*1)))=0;
                    gc = cumsum(g(:,1));
                    xcmod=1+sum(gc(end)*rand>gc);
                    while xcmod>size(xca,1)
                        xcmod=1+sum(gc(end)*rand>gc);
                    end
                    
                    xc=xca(xcmod);
                    
                    if X==size(prodary,1)
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)= profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    elseif X==1
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                    else
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    end
                    clear X
                end

                scrseqary(:,wdindx)=wrkngcol;
                scrprbary(:,wdindx)=prodary(:,indxe+1);
            end
        end

    else
        for wdindx=1:genwidth
            if wdindx==1
                C = cumsum(prodary(:,indxe-1));
                prodary(:,indxe+1)=prodary(:,indxe-1);
                
                for indxg=1:realfreq
                    X=f(1+sum(C(end)*rand>C));
                    if wrkngcol(X,1)==65
                        g(:,1)=sumdifmat(:,2);
                    elseif wrkngcol(X,1)==82
                        g(:,1)=sumdifmat(:,3);
                    elseif wrkngcol(X,1)==78
                        g(:,1)=sumdifmat(:,4);
                    elseif wrkngcol(X,1)==68
                        g(:,1)=sumdifmat(:,5);
                    elseif wrkngcol(X,1)==67
                        g(:,1)=sumdifmat(:,6);
                    elseif wrkngcol(X,1)==81
                        g(:,1)=sumdifmat(:,7);
                    elseif wrkngcol(X,1)==69
                        g(:,1)=sumdifmat(:,8);
                    elseif wrkngcol(X,1)==71
                        g(:,1)=sumdifmat(:,9);
                    elseif wrkngcol(X,1)==72
                        g(:,1)=sumdifmat(:,10);
                    elseif wrkngcol(X,1)==73
                        g(:,1)=sumdifmat(:,11);
                    elseif wrkngcol(X,1)==76
                        g(:,1)=sumdifmat(:,12);
                    elseif wrkngcol(X,1)==75
                        g(:,1)=sumdifmat(:,13);
                    elseif wrkngcol(X,1)==77
                        g(:,1)=sumdifmat(:,14);
                    elseif wrkngcol(X,1)==70
                        g(:,1)=sumdifmat(:,15);
                    elseif wrkngcol(X,1)==80
                        g(:,1)=sumdifmat(:,16);
                    elseif wrkngcol(X,1)==83
                        g(:,1)=sumdifmat(:,17);
                    elseif wrkngcol(X,1)==84
                        g(:,1)=sumdifmat(:,18);
                    elseif wrkngcol(X,1)==87
                        g(:,1)=sumdifmat(:,19);
                    elseif wrkngcol(X,1)==89
                        g(:,1)=sumdifmat(:,20);
                    elseif wrkngcol(X,1)==86
                        g(:,1)=sumdifmat(:,21);
                    end

                    sdevg=std(g(2:21,1));
                    gprep=g(2:21);
                    ming=min(gprep(gprep>0));
                    g(g>(ming+(sdevg*1)))=0;
                    gc = cumsum(g(:,1));
                    xcmod=1+sum(gc(end)*rand>gc);
                    while xcmod>size(xca,1)
                        xcmod=1+sum(gc(end)*rand>gc);
                    end
                    
                    xc=xca(xcmod);
  
                    if X==size(prodary,1)
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)= profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    elseif X==1
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                    else
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    end
                    clear X
                end

                scrseqary(:,wdindx)=wrkngcol;
                scrprbary(:,wdindx)=wrkngcol;
            else
                wrkngcol=wrkbp;
                C = cumsum(prodary(:,indxe-1));
                prodary(:,indxe+1)=prodary(:,indxe-1);
                
                for indxg=1:realfreq
                    X=f(1+sum(C(end)*rand>C));
                    if wrkngcol(X,1)==65
                        g(:,1)=sumdifmat(:,2);
                    elseif wrkngcol(X,1)==82
                        g(:,1)=sumdifmat(:,3);
                    elseif wrkngcol(X,1)==78
                        g(:,1)=sumdifmat(:,4);
                    elseif wrkngcol(X,1)==68
                        g(:,1)=sumdifmat(:,5);
                    elseif wrkngcol(X,1)==67
                        g(:,1)=sumdifmat(:,6);
                    elseif wrkngcol(X,1)==81
                        g(:,1)=sumdifmat(:,7);
                    elseif wrkngcol(X,1)==69
                        g(:,1)=sumdifmat(:,8);
                    elseif wrkngcol(X,1)==71
                        g(:,1)=sumdifmat(:,9);
                    elseif wrkngcol(X,1)==72
                        g(:,1)=sumdifmat(:,10);
                    elseif wrkngcol(X,1)==73
                        g(:,1)=sumdifmat(:,11);
                    elseif wrkngcol(X,1)==76
                        g(:,1)=sumdifmat(:,12);
                    elseif wrkngcol(X,1)==75
                        g(:,1)=sumdifmat(:,13);
                    elseif wrkngcol(X,1)==77
                        g(:,1)=sumdifmat(:,14);
                    elseif wrkngcol(X,1)==70
                        g(:,1)=sumdifmat(:,15);
                    elseif wrkngcol(X,1)==80
                        g(:,1)=sumdifmat(:,16);
                    elseif wrkngcol(X,1)==83
                        g(:,1)=sumdifmat(:,17);
                    elseif wrkngcol(X,1)==84
                        g(:,1)=sumdifmat(:,18);
                    elseif wrkngcol(X,1)==87
                        g(:,1)=sumdifmat(:,19);
                    elseif wrkngcol(X,1)==89
                        g(:,1)=sumdifmat(:,20);
                    elseif wrkngcol(X,1)==86
                        g(:,1)=sumdifmat(:,21);
                    end
 
                    sdevg=std(g(2:21,1));
                    gprep=g(2:21);
                    ming=min(gprep(gprep>0));
                    g(g>(ming+(sdevg*1)))=0;
                    gc = cumsum(g(:,1));
                    xcmod=1+sum(gc(end)*rand>gc);
                    while xcmod>size(xca,1)
                        xcmod=1+sum(gc(end)*rand>gc);
                    end
                    
                    xc=xca(xcmod);
                    
                    if X==size(prodary,1)
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)= profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    elseif X==1
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                    else
                        wrkngcol(X,1)=xc;
                        prodary(X,indxe+1)=profrefmat(find(profrefmat(:,1)==xc),2);
                        prodary(X+1,indxe+1)=prodary(X+1,indxe+1)*1.075;
                        prodary(X-1,indxe+1)=prodary(X-1,indxe+1)*1.075;
                    end
                    clear X
                end
                scrseqary(:,wdindx)=wrkngcol;
                scrprbary(:,wdindx)=prodary(:,indxe+1);
            end
            
            simsearchres=sortrows(importdata('simsearchres'),-1);
            simsearchres=unique(simsearchres,'rows');
            simsearchres=sortrows(simsearchres,-1);
            simsearchres=transpose(simsearchres);
            wrkngcol=simsearchres(4:end,ceil(wdindx/100));
            wrkbp=wrkngcol;
            
        end
    end
    
    clear simsearchres

    csvwrite('altseqary.csv', scrseqary);
    scrseqarytrans=transpose(scrseqary);
    scrseqaryfinproc=zeros(size(scrseqarytrans,1)+1,size(scrseqarytrans,2));
    scrseqaryfinproc(2:size(scrseqaryfinproc,1),:)=scrseqarytrans;
    scrseqarytaba=char(scrseqaryfinproc);
    scrseqarytab=cellstr(scrseqarytaba);
    sILOG=1;
    scrseqarytabtmp=num2cell(zeros(2*genwidth+1,1));
    for indxh=1:2:size(scrseqarytabtmp,1)-1;
        scrseqarytabtmp{indxh,1}=sprintf('>Alt%d', sILOG);
        scrseqarytabtmp{indxh+1,1}=scrseqarytab{sILOG+1,1};
        sILOG=sILOG+1;
    end
    
    for indxi=1:size(scrseqarytabtmp,1)
        if scrseqarytabtmp{indxi,1}==0
            scrseqarytabtmp(indxi,:)=[];
        end
    end
    scrseqarytab=scrseqarytabtmp;
    clear scrseqarytabtmp indxh indxi
    GENALTS=cell2table(scrseqarytab);
    textfilenamegenset = ['gen' num2str(gencounter+1) 'alts.txt'];
    writetable(GENALTS,textfilenamegenset);
    
    letsgo=gencounter+1;

profilecreatorp(letsgo);
parpool(40);
parfor SHAZAAM=1:40
   if SHAZAAM==1; RCLSTSCP001(letsgo);end
   if SHAZAAM==2; RCLSTSCP002(letsgo);end
   if SHAZAAM==3; RCLSTSCP003(letsgo);end
   if SHAZAAM==4; RCLSTSCP004(letsgo);end
end
disp(num2str(gencounter+1));
disp('Calling Gof2Score');
gof2score(letsgo);

    simsearchres=importdata('simsearchres');
    simsearchres=unique(simsearchres,'rows');
 
    simsearchres=transpose(simsearchres);
    
    prodary(:,indxe)=simsearchres(4:end,1);
    prodaryaccomp(:,indxe)=simsearchres(1:3,1);
    simsrchchar=cellstr(transpose(char(simsearchres(4:end,:))));
    simsscores=cellstr(num2str(transpose(simsearchres(1:3,:))));
    simcomb=[simsscores simsrchchar];
    
    SIMALTS=cell2table(simcomb);
    textfilenamesim = ['gen' num2str(gencounter+1) 'scrs.txt'];
    writetable(SIMALTS,textfilenamesim);
    
    movefile *scrs.txt evolvescrs
    
    gencounter=gencounter+1;
end

prodarytrans=transpose(prodary);
prodaryfinproc=zeros(size(prodarytrans,1)+1,size(prodarytrans,2));
prodaryfinproc(2:size(prodaryfinproc,1),:)=prodarytrans;
prodarytaba=char(prodaryfinproc);
prodarytab=cellstr(prodarytaba);
ILOG=1;

for indxh=1:2:size(prodarytab,1)-1
    prodarytab{indxh,1}=sprintf('>Iteration%d', ILOG);
    ILOG=ILOG+1;
end

prodaryaccomptrans=transpose(prodaryaccomp);
prodaryaccompfinproc=zeros(size(prodaryaccomptrans,1)+1,size(prodaryaccomptrans,2));
prodaryaccompfinproc(2:size(prodaryaccompfinproc,1),:)=prodaryaccomptrans;
prodaryaccomptaba=char(prodaryaccompfinproc);
prodaryaccomptab=cellstr(prodaryaccomptaba);
xILOG=1;

for indxfil=1:2:size(prodaryaccomptab,1)-1
    prodaryaccomptab{indxh,1}=sprintf('>Iteration%d', ILOG);
    xILOG=xILOG+1;
end

PAT=cell2table(prodarytab);
writetable(PAT,'PrimarySequences.txt');
PATAC=cell2table(prodaryaccomptab);
writetable(PATAC,'PrimarySequenceScores.txt');

runstats=num2cell(zeros(6,2));
runstats{1,1}='StartSeq';
runstats{1,2}=STRTTB;
runstats{2,1}='StartSeqMod';
runstats{2,2}=transpose(char(STRTTB2b));
runstats{3,1}='SdevProp';
runstats{3,2}=1;
runstats{4,1}='Iterations';
runstats{4,2}=generations+1;
runstats{5,1}='Mutation Freq';
runstats{5,2}=freq;
runstats{6,1}='Absolute Mutation Freq';
runstats{6,2}=realfreq;

RST=cell2table(runstats);
writetable(RST,'runstats.csv');

mkdir results

movefile PrimarySequences.txt results
movefile PrimarySequenceScores.txt results
movefile runstats.csv results



