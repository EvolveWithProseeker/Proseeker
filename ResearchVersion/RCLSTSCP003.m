function [letsgo]= RCLSTSCP003(letsgo) 

%Function:          Search Profiles and then compare selected windows from
%                   the protein under assessment to them.

%                   This replaces the deprecated simsearch and reclustscorePAR functions.

%Inputs:            Motif files (previously generated via interative
%                   profilecreator and then turned into motif sets by
%                   searchderive)
%                   Profiles to search against (created by profilecreatorp)

%Outputs            Goodness of fit statistics for individual proteins based on
%                   functional motif detection and for collating by gof2score

coder.extrinsic('join');
coder.extrinsic('mex_WriteMatrix');

finscrmatval=1; %#ok<*NASGU> 
files={"p1.bres1.x.csv";"p1.bres2.x.csv";"p2.bres1.x.csv";"p2.bres2.x.csv";"p3.bres1.x.csv";"p3.bres2.x.csv";"p4.bres1.x.csv";"p4.bres2.x.csv";"p5.bres1.x.csv";"p5.bres2.x.csv";"p6.bres1.x.csv";"p6.bres2.x.csv";"p7.bres1.x.csv";"p7.bres2.x.csv";"p8.bres1.x.csv";"p8.bres2.x.csv";"p9.bres1.x.csv";"p9.bres2.x.csv";"p10.bres1.x.csv";"p10.bres2.x.csv";"p11.bres1.x.csv";"p11.bres2.x.csv";"p12.bres1.x.csv";"p12.bres2.x.csv";"p13.bres1.x.csv";"p13.bres2.x.csv";"p14.bres1.x.csv";"p14.bres2.x.csv";"p15.bres1.x.csv";"p15.bres2.x.csv";"p16.bres1.x.csv";"p16.bres2.x.csv";"p17.bres1.x.csv";"p17.bres2.x.csv";"p18.bres1.x.csv";"p18.bres2.x.csv";"p19.bres1.x.csv";"p19.bres2.x.csv";"p20.bres1.x.csv";"p20.bres2.x.csv";"p21.bres1.x.csv";"p21.bres2.x.csv";"p22.bres1.x.csv";"p22.bres2.x.csv";"p23.bres1.x.csv";"p23.bres2.x.csv";"p24.bres1.x.csv";"p24.bres2.x.csv";"p25.bres1.x.csv";"p25.bres2.x.csv";"p26.bres1.x.csv";"p26.bres2.x.csv";"p27.bres1.x.csv";"p27.bres2.x.csv";"p28.bres1.x.csv";"p28.bres2.x.csv";"p29.bres1.x.csv";"p29.bres2.x.csv";"p30.bres1.x.csv";"p30.bres2.x.csv"};
flcount=size(files,1);
goodnessmatfin=zeros(50,16245);
assdsdcountsfin=[0;0;0];
blockcpy=coder.load('BLOCK.txt');
goodnessmat1=zeros(50,361);
rngmat=coder.load('rankingschemeMOD.csvspc');
rngmatstr=zeros(size(rngmat,1),3);
for rngindx=1:size(rngmat,1)
    rngmatstr(rngindx,1)=min(rngmat(rngindx,1:20)); %#ok<PFBNS>
end
for rngindx=1:size(rngmat,1)
    rngmatstr(rngindx,2)=max(rngmat(rngindx,1:20)); %#ok<PFBNS>
end
rngmatstr(:,3)=(rngmatstr(:,1)-rngmatstr(:,2));
posinrang=zeros(544,15*60);
meanwavs=zeros(100,13*60);
fourierwavs=zeros(50,13*60);
cmeanmat=zeros(50,1*60);
submatGOF=zeros(50,6*60);
posinrangcnt=16;
meanwavscnt=14;
fourierwavscnt=14;
cmeanmatcnt=2;
submatGOFcnt=7;

tmpsrc=zeros(694,210);
for bresindex=1:60
    if bresindex==1
        targetload=coder.load('p1.bres1.x.csv');
        tmpsrc(:,:)=targetload(:,:);
        tmpsrc3=tmpsrc(51:150,2:14);
        tmpsrc7=tmpsrc(151:694,1:15);
        tmpsrc9=tmpsrc(202:251,17:29);
        tmpsrc10=tmpsrc(152:201,26);
        tmpsrc11=tmpsrc(152:201,[21,22,23,24,25,27]);
        posinrang(:,1:15)=tmpsrc7;
        meanwavs(:,1:13)=tmpsrc3;
        fourierwavs(:,1:13)=tmpsrc9;
        cmeanmat(:,1)=tmpsrc10;
        submatGOF(:,1:6)=tmpsrc11;
    else
        switch bresindex
            case 2
                targetload=coder.load('p1.bres2.x.csv');
            case 3
                targetload=coder.load('p2.bres1.x.csv');
            case 4
                targetload=coder.load('p2.bres2.x.csv');
            case 5
                targetload=coder.load('p3.bres1.x.csv');
            case 6
                targetload=coder.load('p3.bres2.x.csv');
            case 7
                targetload=coder.load('p4.bres1.x.csv');
            case 8
                targetload=coder.load('p4.bres2.x.csv');
            case 9
                targetload=coder.load('p5.bres1.x.csv');
            case 10
                targetload=coder.load('p5.bres2.x.csv');
            case 11
                targetload=coder.load('p6.bres1.x.csv');
            case 12
                targetload=coder.load('p6.bres2.x.csv');
            case 13
                targetload=coder.load('p7.bres1.x.csv');
            case 14
                targetload=coder.load('p7.bres2.x.csv');
            case 15
                targetload=coder.load('p8.bres1.x.csv');
            case 16
                targetload=coder.load('p8.bres2.x.csv');
            case 17
                targetload=coder.load('p9.bres1.x.csv');
            case 18
                targetload=coder.load('p9.bres2.x.csv');
            case 19
                targetload=coder.load('p10.bres1.x.csv');
            case 20
                targetload=coder.load('p10.bres2.x.csv');
            case 21
                targetload=coder.load('p11.bres1.x.csv');
            case 22
                targetload=coder.load('p11.bres2.x.csv');
            case 23
                targetload=coder.load('p12.bres1.x.csv');
            case 24
                targetload=coder.load('p12.bres2.x.csv');
            case 25
                targetload=coder.load('p13.bres1.x.csv');
            case 26
                targetload=coder.load('p13.bres2.x.csv');
            case 27
                targetload=coder.load('p14.bres1.x.csv');
            case 28
                targetload=coder.load('p14.bres2.x.csv');
            case 29
                targetload=coder.load('p15.bres1.x.csv');
            case 30
                targetload=coder.load('p15.bres2.x.csv');
            case 31
                targetload=coder.load('p16.bres1.x.csv');
            case 32
                targetload=coder.load('p16.bres2.x.csv');
            case 33
                targetload=coder.load('p17.bres1.x.csv');
            case 34
                targetload=coder.load('p17.bres2.x.csv');
            case 35
                targetload=coder.load('p18.bres1.x.csv');
            case 36
                targetload=coder.load('p18.bres2.x.csv');
            case 37
                targetload=coder.load('p19.bres1.x.csv');
            case 38
                targetload=coder.load('p19.bres2.x.csv');
            case 39
                targetload=coder.load('p20.bres1.x.csv');
            case 40
                targetload=coder.load('p20.bres2.x.csv');
            case 41
                targetload=coder.load('p21.bres1.x.csv');
            case 42
                targetload=coder.load('p21.bres2.x.csv');
            case 43
                targetload=coder.load('p22.bres1.x.csv');
            case 44
                targetload=coder.load('p22.bres2.x.csv');
            case 45
                targetload=coder.load('p23.bres1.x.csv');
            case 46
                targetload=coder.load('p23.bres2.x.csv');
            case 47
                targetload=coder.load('p24.bres1.x.csv');
            case 48
                targetload=coder.load('p24.bres2.x.csv');
            case 49
                targetload=coder.load('p25.bres1.x.csv');
            case 50
                targetload=coder.load('p25.bres2.x.csv');
            case 51
                targetload=coder.load('p26.bres1.x.csv');
            case 52
                targetload=coder.load('p26.bres2.x.csv');
            case 53
                targetload=coder.load('p27.bres1.x.csv');
            case 54
                targetload=coder.load('p27.bres2.x.csv');
            case 55
                targetload=coder.load('p28.bres1.x.csv');
            case 56
                targetload=coder.load('p28.bres2.x.csv');
            case 57
                targetload=coder.load('p29.bres1.x.csv');
            case 58
                targetload=coder.load('p29.bres2.x.csv');
            case 59
                targetload=coder.load('p30.bres1.x.csv');
            case 60
                targetload=coder.load('p30.bres2.x.csv');
        end
    tmpsrc(:,:)=targetload(:,:);
    tmpsrc3=tmpsrc(51:150,2:14);
    tmpsrc7=tmpsrc(151:694,1:15);
    tmpsrc9=tmpsrc(202:251,17:29);
    tmpsrc10=tmpsrc(152:201,26);
    tmpsrc11=tmpsrc(152:201,[21,22,23,24,25,27]);
    posinrang(:,posinrangcnt:posinrangcnt+14)=tmpsrc7;
    posinrangcnt=posinrangcnt+15;
    meanwavs(:,meanwavscnt:meanwavscnt+12)=tmpsrc3;
    meanwavscnt=meanwavscnt+13;
    fourierwavs(:,fourierwavscnt:fourierwavscnt+12)=tmpsrc9;
    fourierwavscnt=fourierwavscnt+13;
    cmeanmat(:,cmeanmatcnt)=tmpsrc10;
    cmeanmatcnt=cmeanmatcnt+1;
    submatGOF(:,submatGOFcnt:submatGOFcnt+5)=tmpsrc11;
    submatGOFcnt=submatGOFcnt+6;
    end
end

tmpcountmant=zeros(50,1);
submatGOF=[tmpcountmant submatGOF];

for switchindex=1:25
    switch switchindex
        case 1
            subject=coder.load('p51.csvx');
        case 2
            subject=coder.load('p52.csvx');
        case 3
            subject=coder.load('p53.csvx');
        case 4
            subject=coder.load('p54.csvx');
        case 5
            subject=coder.load('p55.csvx');
        case 6
            subject=coder.load('p56.csvx');
        case 7
            subject=coder.load('p57.csvx');
        case 8
            subject=coder.load('p58.csvx');
        case 9
            subject=coder.load('p59.csvx');
        case 10
            subject=coder.load('p60.csvx');
        case 11
            subject=coder.load('p61.csvx');
        case 12
            subject=coder.load('p62.csvx');
        case 13
            subject=coder.load('p63.csvx');
        case 14
            subject=coder.load('p64.csvx');
        case 15
            subject=coder.load('p65.csvx');
        case 16
            subject=coder.load('p66.csvx');
        case 17
            subject=coder.load('p67.csvx');
        case 18
            subject=coder.load('p68.csvx');
        case 19
            subject=coder.load('p69.csvx');
        case 20
            subject=coder.load('p70.csvx');
        case 21
            subject=coder.load('p71.csvx');
        case 22
            subject=coder.load('p72.csvx');
        case 23
            subject=coder.load('p73.csvx');
        case 24
            subject=coder.load('p74.csvx');
        case 25
            subject=coder.load('p75.csvx');
    end
    
    for subsrchindx1=1:size(subject,1)
        for subsubsrchindx1=1:size(subject,2)
            subject(subsrchindx1,subsubsrchindx1)=(subject(subsrchindx1,subsubsrchindx1)-rngmatstr(subsrchindx1,1))/rngmatstr(subsrchindx1,3);
        end
    end
    
    subject=[subject(:,end-9:end) subject subject(:,1:10)];
    
    circblock=[blockcpy(:,end-9:end) blockcpy blockcpy(:,1:10)];
    
    for intindx1=11:size(subject,2)-10
        if circblock(1,intindx1)==1
            
            tblock=subject(:,intindx1-6:intindx1+6);
            
            for clustass=15:15:size(posinrang,2)
                subjasstblock=[posinrang(:,1) tblock posinrang(:,clustass)];
                subjasstblock=sortrows(subjasstblock,size(tblock,2));
                
                c1=subjasstblock;
                c2=subjasstblock;
                c3=subjasstblock;
                c4=subjasstblock;
                c5=subjasstblock;
                c6=subjasstblock;
                c7=subjasstblock;
                c8=subjasstblock;
                c9=subjasstblock;
                c10=subjasstblock;
                c11=subjasstblock;
                c12=subjasstblock;
                c13=subjasstblock;
                c14=subjasstblock;
                c15=subjasstblock;
                c16=subjasstblock;
                c17=subjasstblock;
                c18=subjasstblock;
                c19=subjasstblock;
                c20=subjasstblock;
                c21=subjasstblock;
                c22=subjasstblock;
                c23=subjasstblock;
                c24=subjasstblock;
                c25=subjasstblock;
                c26=subjasstblock;
                c27=subjasstblock;
                c28=subjasstblock;
                c29=subjasstblock;
                c30=subjasstblock;
                c31=subjasstblock;
                c32=subjasstblock;
                c33=subjasstblock;
                c34=subjasstblock;
                c35=subjasstblock;
                c36=subjasstblock;
                c37=subjasstblock;
                c38=subjasstblock;
                c39=subjasstblock;
                c40=subjasstblock;
                c41=subjasstblock;
                c42=subjasstblock;
                c43=subjasstblock;
                c44=subjasstblock;
                c45=subjasstblock;
                c46=subjasstblock;
                c47=subjasstblock;
                c48=subjasstblock;
                c49=subjasstblock;
                c50=subjasstblock;
                
                coder.varsize('c1');
                coder.varsize('c2');
                coder.varsize('c3');
                coder.varsize('c4');
                coder.varsize('c5');
                coder.varsize('c6');
                coder.varsize('c7');
                coder.varsize('c8');
                coder.varsize('c9');
                coder.varsize('c10');
                coder.varsize('c11');
                coder.varsize('c12');
                coder.varsize('c13');
                coder.varsize('c14');
                coder.varsize('c15');
                coder.varsize('c16');
                coder.varsize('c17');
                coder.varsize('c18');
                coder.varsize('c19');
                coder.varsize('c20');
                coder.varsize('c21');
                coder.varsize('c22');
                coder.varsize('c23');
                coder.varsize('c24');
                coder.varsize('c25');
                coder.varsize('c26');
                coder.varsize('c27');
                coder.varsize('c28');
                coder.varsize('c29');
                coder.varsize('c30');
                coder.varsize('c31');
                coder.varsize('c32');
                coder.varsize('c33');
                coder.varsize('c34');
                coder.varsize('c35');
                coder.varsize('c36');
                coder.varsize('c37');
                coder.varsize('c38');
                coder.varsize('c39');
                coder.varsize('c40');
                coder.varsize('c41');
                coder.varsize('c42');
                coder.varsize('c43');
                coder.varsize('c44');
                coder.varsize('c45');
                coder.varsize('c46');
                coder.varsize('c47');
                coder.varsize('c48');
                coder.varsize('c49');
                coder.varsize('c50');
                
                for cutindx=1:size(c1,1)
                    if c1(cutindx,size(c1,2))~=1
                        c1(cutindx,:)=0;
                    end
                    if c2(cutindx,size(c2,2))~=2
                        c2(cutindx,:)=0;
                    end
                    if c3(cutindx,size(c3,2))~=3
                        c3(cutindx,:)=0;
                    end
                    if c4(cutindx,size(c4,2))~=4
                        c4(cutindx,:)=0;
                    end
                    if c5(cutindx,size(c5,2))~=5
                        c5(cutindx,:)=0;
                    end
                    if c6(cutindx,size(c6,2))~=6
                        c6(cutindx,:)=0;
                    end
                    if c7(cutindx,size(c7,2))~=7
                        c7(cutindx,:)=0;
                    end
                    if c8(cutindx,size(c8,2))~=8
                        c8(cutindx,:)=0;
                    end
                    if c9(cutindx,size(c9,2))~=9
                        c9(cutindx,:)=0;
                    end
                    if c10(cutindx,size(c10,2))~=10
                        c10(cutindx,:)=0;
                    end
                    if c11(cutindx,size(c11,2))~=11
                        c11(cutindx,:)=0;
                    end
                    if c12(cutindx,size(c12,2))~=12
                        c12(cutindx,:)=0;
                    end
                    if c13(cutindx,size(c13,2))~=13
                        c13(cutindx,:)=0;
                    end
                    if c14(cutindx,size(c14,2))~=14
                        c14(cutindx,:)=0;
                    end
                    if c15(cutindx,size(c15,2))~=15
                        c15(cutindx,:)=0;
                    end
                    if c16(cutindx,size(c16,2))~=16
                        c16(cutindx,:)=0;
                    end
                    if c17(cutindx,size(c17,2))~=17
                        c17(cutindx,:)=0;
                    end
                    if c18(cutindx,size(c18,2))~=18
                        c18(cutindx,:)=0;
                    end
                    if c19(cutindx,size(c19,2))~=19
                        c19(cutindx,:)=0;
                    end
                    if c20(cutindx,size(c20,2))~=20
                        c20(cutindx,:)=0;
                    end
                    if c21(cutindx,size(c21,2))~=21
                        c21(cutindx,:)=0;
                    end
                    if c22(cutindx,size(c22,2))~=22
                        c22(cutindx,:)=0;
                    end
                    if c23(cutindx,size(c23,2))~=23
                        c23(cutindx,:)=0;
                    end
                    if c24(cutindx,size(c24,2))~=24
                        c24(cutindx,:)=0;
                    end
                    if c25(cutindx,size(c25,2))~=25
                        c25(cutindx,:)=0;
                    end
                    if c26(cutindx,size(c26,2))~=26
                        c26(cutindx,:)=0;
                    end
                    if c27(cutindx,size(c27,2))~=27
                        c27(cutindx,:)=0;
                    end
                    if c28(cutindx,size(c28,2))~=28
                        c28(cutindx,:)=0;
                    end
                    if c29(cutindx,size(c29,2))~=29
                        c29(cutindx,:)=0;
                    end
                    if c30(cutindx,size(c30,2))~=30
                        c30(cutindx,:)=0;
                    end
                    if c31(cutindx,size(c31,2))~=31
                        c31(cutindx,:)=0;
                    end
                    if c32(cutindx,size(c32,2))~=32
                        c32(cutindx,:)=0;
                    end
                    if c33(cutindx,size(c33,2))~=33
                        c33(cutindx,:)=0;
                    end
                    if c34(cutindx,size(c34,2))~=34
                        c34(cutindx,:)=0;
                    end
                    if c35(cutindx,size(c35,2))~=35
                        c35(cutindx,:)=0;
                    end
                    if c36(cutindx,size(c36,2))~=36
                        c36(cutindx,:)=0;
                    end
                    if c37(cutindx,size(c37,2))~=37
                        c37(cutindx,:)=0;
                    end
                    if c38(cutindx,size(c38,2))~=38
                        c38(cutindx,:)=0;
                    end
                    if c39(cutindx,size(c39,2))~=39
                        c39(cutindx,:)=0;
                    end
                    if c40(cutindx,size(c40,2))~=40
                        c40(cutindx,:)=0;
                    end
                    if c41(cutindx,size(c41,2))~=41
                        c41(cutindx,:)=0;
                    end
                    if c42(cutindx,size(c42,2))~=42
                        c42(cutindx,:)=0;
                    end
                    if c43(cutindx,size(c43,2))~=43
                        c43(cutindx,:)=0;
                    end
                    if c44(cutindx,size(c44,2))~=44
                        c44(cutindx,:)=0;
                    end
                    if c45(cutindx,size(c45,2))~=45
                        c45(cutindx,:)=0;
                    end
                    if c46(cutindx,size(c46,2))~=46
                        c46(cutindx,:)=0;
                    end
                    if c47(cutindx,size(c47,2))~=47
                        c47(cutindx,:)=0;
                    end
                    if c48(cutindx,size(c48,2))~=48
                        c48(cutindx,:)=0;
                    end
                    if c49(cutindx,size(c49,2))~=49
                        c49(cutindx,:)=0;
                    end
                    if c50(cutindx,size(c50,2))~=50
                        c50(cutindx,:)=0;
                    end
                end
                
                c1( ~any(c1,2), : ) = [];
                c2( ~any(c2,2), : ) = [];
                c3( ~any(c3,2), : ) = [];
                c4( ~any(c4,2), : ) = [];
                c5( ~any(c5,2), : ) = [];
                c6( ~any(c6,2), : ) = [];
                c7( ~any(c7,2), : ) = [];
                c8( ~any(c8,2), : ) = [];
                c9( ~any(c9,2), : ) = [];
                c10( ~any(c10,2), : ) = [];
                c11( ~any(c11,2), : ) = [];
                c12( ~any(c12,2), : ) = [];
                c13( ~any(c13,2), : ) = [];
                c14( ~any(c14,2), : ) = [];
                c15( ~any(c15,2), : ) = [];
                c16( ~any(c16,2), : ) = [];
                c17( ~any(c17,2), : ) = [];
                c18( ~any(c18,2), : ) = [];
                c19( ~any(c19,2), : ) = [];
                c20( ~any(c20,2), : ) = [];
                c21( ~any(c21,2), : ) = [];
                c22( ~any(c22,2), : ) = [];
                c23( ~any(c23,2), : ) = [];
                c24( ~any(c24,2), : ) = [];
                c25( ~any(c25,2), : ) = [];
                c26( ~any(c26,2), : ) = [];
                c27( ~any(c27,2), : ) = [];
                c28( ~any(c28,2), : ) = [];
                c29( ~any(c29,2), : ) = [];
                c30( ~any(c30,2), : ) = [];
                c31( ~any(c31,2), : ) = [];
                c32( ~any(c32,2), : ) = [];
                c33( ~any(c33,2), : ) = [];
                c34( ~any(c34,2), : ) = [];
                c35( ~any(c35,2), : ) = [];
                c36( ~any(c36,2), : ) = [];
                c37( ~any(c37,2), : ) = [];
                c38( ~any(c38,2), : ) = [];
                c39( ~any(c39,2), : ) = [];
                c40( ~any(c40,2), : ) = [];
                c41( ~any(c41,2), : ) = [];
                c42( ~any(c42,2), : ) = [];
                c43( ~any(c43,2), : ) = [];
                c44( ~any(c44,2), : ) = [];
                c45( ~any(c45,2), : ) = [];
                c46( ~any(c46,2), : ) = [];
                c47( ~any(c47,2), : ) = [];
                c48( ~any(c48,2), : ) = [];
                c49( ~any(c49,2), : ) = [];
                c50( ~any(c50,2), : ) = [];
                
                cextra=zeros(2,15);
                c1=[c1; cextra];
                c2=[c2; cextra];
                c3=[c3; cextra];
                c4=[c4; cextra];
                c5=[c5; cextra];
                c6=[c6; cextra];
                c7=[c7; cextra];
                c8=[c8; cextra];
                c9=[c9; cextra];
                c10=[c10; cextra];
                c11=[c11; cextra];
                c12=[c12; cextra];
                c13=[c13; cextra];
                c14=[c14; cextra];
                c15=[c15; cextra];
                c16=[c16; cextra];
                c17=[c17; cextra];
                c18=[c18; cextra];
                c19=[c19; cextra];
                c20=[c20; cextra];
                c21=[c21; cextra];
                c22=[c22; cextra];
                c23=[c23; cextra];
                c24=[c24; cextra];
                c25=[c25; cextra];
                c26=[c26; cextra];
                c27=[c27; cextra];
                c28=[c28; cextra];
                c29=[c29; cextra];
                c30=[c30; cextra];
                c31=[c31; cextra];
                c32=[c32; cextra];
                c33=[c33; cextra];
                c34=[c34; cextra];
                c35=[c35; cextra];
                c36=[c36; cextra];
                c37=[c37; cextra];
                c38=[c38; cextra];
                c39=[c39; cextra];
                c40=[c40; cextra];
                c41=[c41; cextra];
                c42=[c42; cextra];
                c43=[c43; cextra];
                c44=[c44; cextra];
                c45=[c45; cextra];
                c46=[c46; cextra];
                c47=[c47; cextra];
                c48=[c48; cextra];
                c49=[c49; cextra];
                c50=[c50; cextra];
                
                c1sz=size(c1,1);
                c1(c1sz-1:c1sz,15)=c1(c1sz-2,15)-(2*c1(c1sz-2,15));
                c1(c1sz,15)=c1(c1sz,15)-0.5;
                c2sz=size(c2,1);
                c2(c2sz-1:c2sz,15)=c2(c2sz-2,15)-(2*c2(c2sz-2,15));
                c2(c2sz,15)=c2(c2sz,15)-0.5;
                c3sz=size(c3,1);
                c3(c3sz-1:c3sz,15)=c3(c3sz-2,15)-(2*c3(c3sz-2,15));
                c3(c3sz,15)=c3(c3sz,15)-0.5;
                c4sz=size(c4,1);
                c4(c4sz-1:c4sz,15)=c4(c4sz-2,15)-(2*c4(c4sz-2,15));
                c4(c4sz,15)=c4(c4sz,15)-0.5;
                c5sz=size(c5,1);
                c5(c5sz-1:c5sz,15)=c5(c5sz-2,15)-(2*c5(c5sz-2,15));
                c5(c5sz,15)=c5(c5sz,15)-0.5;
                c6sz=size(c6,1);
                c6(c6sz-1:c6sz,15)=c6(c6sz-2,15)-(2*c6(c6sz-2,15));
                c6(c6sz,15)=c6(c6sz,15)-0.5;
                c7sz=size(c7,1);
                c7(c7sz-1:c7sz,15)=c7(c7sz-2,15)-(2*c7(c7sz-2,15));
                c7(c7sz,15)=c7(c7sz,15)-0.5;
                c8sz=size(c8,1);
                c8(c8sz-1:c8sz,15)=c8(c8sz-2,15)-(2*c8(c8sz-2,15));
                c8(c8sz,15)=c8(c8sz,15)-0.5;
                c9sz=size(c9,1);
                c9(c9sz-1:c9sz,15)=c9(c9sz-2,15)-(2*c9(c9sz-2,15));
                c9(c9sz,15)=c9(c9sz,15)-0.5;
                c10sz=size(c10,1);
                c10(c10sz-1:c10sz,15)=c10(c10sz-2,15)-(2*c10(c10sz-2,15));
                c10(c10sz,15)=c10(c10sz,15)-0.5;
                c11sz=size(c11,1);
                c11(c11sz-1:c11sz,15)=c11(c11sz-2,15)-(2*c11(c11sz-2,15));
                c11(c11sz,15)=c11(c11sz,15)-0.5;
                c12sz=size(c12,1);
                c12(c12sz-1:c12sz,15)=c12(c12sz-2,15)-(2*c12(c12sz-2,15));
                c12(c12sz,15)=c12(c12sz,15)-0.5;
                c13sz=size(c13,1);
                c13(c13sz-1:c13sz,15)=c13(c13sz-2,15)-(2*c13(c13sz-2,15));
                c13(c13sz,15)=c13(c13sz,15)-0.5;
                c14sz=size(c14,1);
                c14(c14sz-1:c14sz,15)=c14(c14sz-2,15)-(2*c14(c14sz-2,15));
                c14(c14sz,15)=c14(c14sz,15)-0.5;
                c15sz=size(c15,1);
                c15(c15sz-1:c15sz,15)=c15(c15sz-2,15)-(2*c15(c15sz-2,15));
                c15(c15sz,15)=c15(c15sz,15)-0.5;
                c16sz=size(c16,1);
                c16(c16sz-1:c16sz,15)=c16(c16sz-2,15)-(2*c16(c16sz-2,15));
                c16(c16sz,15)=c16(c16sz,15)-0.5;
                c17sz=size(c17,1);
                c17(c17sz-1:c17sz,15)=c17(c17sz-2,15)-(2*c17(c17sz-2,15));
                c17(c17sz,15)=c17(c17sz,15)-0.5;
                c18sz=size(c18,1);
                c18(c18sz-1:c18sz,15)=c18(c18sz-2,15)-(2*c18(c18sz-2,15));
                c18(c18sz,15)=c18(c18sz,15)-0.5;
                c19sz=size(c19,1);
                c19(c19sz-1:c19sz,15)=c19(c19sz-2,15)-(2*c19(c19sz-2,15));
                c19(c19sz,15)=c19(c19sz,15)-0.5;
                c20sz=size(c20,1);
                c20(c20sz-1:c20sz,15)=c20(c20sz-2,15)-(2*c20(c20sz-2,15));
                c20(c20sz,15)=c20(c20sz,15)-0.5;
                c21sz=size(c21,1);
                c21(c21sz-1:c21sz,15)=c21(c21sz-2,15)-(2*c21(c21sz-2,15));
                c21(c21sz,15)=c21(c21sz,15)-0.5;
                c22sz=size(c22,1);
                c22(c22sz-1:c22sz,15)=c22(c22sz-2,15)-(2*c22(c22sz-2,15));
                c22(c22sz,15)=c22(c22sz,15)-0.5;
                c23sz=size(c23,1);
                c23(c23sz-1:c23sz,15)=c23(c23sz-2,15)-(2*c23(c23sz-2,15));
                c23(c23sz,15)=c23(c23sz,15)-0.5;
                c24sz=size(c24,1);
                c24(c24sz-1:c24sz,15)=c24(c24sz-2,15)-(2*c24(c24sz-2,15));
                c24(c24sz,15)=c24(c24sz,15)-0.5;
                c25sz=size(c25,1);
                c25(c25sz-1:c25sz,15)=c25(c25sz-2,15)-(2*c25(c25sz-2,15));
                c25(c25sz,15)=c25(c25sz,15)-0.5;
                c26sz=size(c26,1);
                c26(c26sz-1:c26sz,15)=c26(c26sz-2,15)-(2*c26(c26sz-2,15));
                c26(c26sz,15)=c26(c26sz,15)-0.5;
                c27sz=size(c27,1);
                c27(c27sz-1:c27sz,15)=c27(c27sz-2,15)-(2*c27(c27sz-2,15));
                c27(c27sz,15)=c27(c27sz,15)-0.5;
                c28sz=size(c28,1);
                c28(c28sz-1:c28sz,15)=c28(c28sz-2,15)-(2*c28(c28sz-2,15));
                c28(c28sz,15)=c28(c28sz,15)-0.5;
                c29sz=size(c29,1);
                c29(c29sz-1:c29sz,15)=c29(c29sz-2,15)-(2*c29(c29sz-2,15));
                c29(c29sz,15)=c29(c29sz,15)-0.5;
                c30sz=size(c30,1);
                c30(c30sz-1:c30sz,15)=c30(c30sz-2,15)-(2*c30(c30sz-2,15));
                c30(c30sz,15)=c30(c30sz,15)-0.5;
                c31sz=size(c31,1);
                c31(c31sz-1:c31sz,15)=c31(c31sz-2,15)-(2*c31(c31sz-2,15));
                c31(c31sz,15)=c31(c31sz,15)-0.5;
                c32sz=size(c32,1);
                c32(c32sz-1:c32sz,15)=c32(c32sz-2,15)-(2*c32(c32sz-2,15));
                c32(c32sz,15)=c32(c32sz,15)-0.5;
                c33sz=size(c33,1);
                c33(c33sz-1:c33sz,15)=c33(c33sz-2,15)-(2*c33(c33sz-2,15));
                c33(c33sz,15)=c33(c33sz,15)-0.5;
                c34sz=size(c34,1);
                c34(c34sz-1:c34sz,15)=c34(c34sz-2,15)-(2*c34(c34sz-2,15));
                c34(c34sz,15)=c34(c34sz,15)-0.5;
                c35sz=size(c35,1);
                c35(c35sz-1:c35sz,15)=c35(c35sz-2,15)-(2*c35(c35sz-2,15));
                c35(c35sz,15)=c35(c35sz,15)-0.5;
                c36sz=size(c36,1);
                c36(c36sz-1:c36sz,15)=c36(c36sz-2,15)-(2*c36(c36sz-2,15));
                c36(c36sz,15)=c36(c36sz,15)-0.5;
                c37sz=size(c37,1);
                c37(c37sz-1:c37sz,15)=c37(c37sz-2,15)-(2*c37(c37sz-2,15));
                c37(c37sz,15)=c37(c37sz,15)-0.5;
                c38sz=size(c38,1);
                c38(c38sz-1:c38sz,15)=c38(c38sz-2,15)-(2*c38(c38sz-2,15));
                c38(c38sz,15)=c38(c38sz,15)-0.5;
                c39sz=size(c39,1);
                c39(c39sz-1:c39sz,15)=c39(c39sz-2,15)-(2*c39(c39sz-2,15));
                c39(c39sz,15)=c39(c39sz,15)-0.5;
                c40sz=size(c40,1);
                c40(c40sz-1:c40sz,15)=c40(c40sz-2,15)-(2*c40(c40sz-2,15));
                c40(c40sz,15)=c40(c40sz,15)-0.5;
                c41sz=size(c41,1);
                c41(c41sz-1:c41sz,15)=c41(c41sz-2,15)-(2*c41(c41sz-2,15));
                c41(c41sz,15)=c41(c41sz,15)-0.5;
                c42sz=size(c42,1);
                c42(c42sz-1:c42sz,15)=c42(c42sz-2,15)-(2*c42(c42sz-2,15));
                c42(c42sz,15)=c42(c42sz,15)-0.5;
                c43sz=size(c43,1);
                c43(c43sz-1:c43sz,15)=c43(c43sz-2,15)-(2*c43(c43sz-2,15));
                c43(c43sz,15)=c43(c43sz,15)-0.5;
                c44sz=size(c44,1);
                c44(c44sz-1:c44sz,15)=c44(c44sz-2,15)-(2*c44(c44sz-2,15));
                c44(c44sz,15)=c44(c44sz,15)-0.5;
                c45sz=size(c45,1);
                c45(c45sz-1:c45sz,15)=c45(c45sz-2,15)-(2*c45(c45sz-2,15));
                c45(c45sz,15)=c45(c45sz,15)-0.5;
                c46sz=size(c46,1);
                c46(c46sz-1:c46sz,15)=c46(c46sz-2,15)-(2*c46(c46sz-2,15));
                c46(c46sz,15)=c46(c46sz,15)-0.5;
                c47sz=size(c47,1);
                c47(c47sz-1:c47sz,15)=c47(c47sz-2,15)-(2*c47(c47sz-2,15));
                c47(c47sz,15)=c47(c47sz,15)-0.5;
                c48sz=size(c48,1);
                c48(c48sz-1:c48sz,15)=c48(c48sz-2,15)-(2*c48(c48sz-2,15));
                c48(c48sz,15)=c48(c48sz,15)-0.5;
                c49sz=size(c49,1);
                c49(c49sz-1:c49sz,15)=c49(c49sz-2,15)-(2*c49(c49sz-2,15));
                c49(c49sz,15)=c49(c49sz,15)-0.5;
                c50sz=size(c50,1);
                c50(c50sz-1:c50sz,15)=c50(c50sz-2,15)-(2*c50(c50sz-2,15));
                c50(c50sz,15)=c50(c50sz,15)-0.5;
                
                clustsizemat=zeros(50,196);
                clustsizemat(1,2)=c1sz;
                clustsizemat(2,2)=c2sz;
                clustsizemat(3,2)=c3sz;
                clustsizemat(4,2)=c4sz;
                clustsizemat(5,2)=c5sz;
                clustsizemat(6,2)=c6sz;
                clustsizemat(7,2)=c7sz;
                clustsizemat(8,2)=c8sz;
                clustsizemat(9,2)=c9sz;
                clustsizemat(10,2)=c10sz;
                clustsizemat(11,2)=c11sz;
                clustsizemat(12,2)=c12sz;
                clustsizemat(13,2)=c13sz;
                clustsizemat(14,2)=c14sz;
                clustsizemat(15,2)=c15sz;
                clustsizemat(16,2)=c16sz;
                clustsizemat(17,2)=c17sz;
                clustsizemat(18,2)=c18sz;
                clustsizemat(19,2)=c19sz;
                clustsizemat(20,2)=c20sz;
                clustsizemat(21,2)=c21sz;
                clustsizemat(22,2)=c22sz;
                clustsizemat(23,2)=c23sz;
                clustsizemat(24,2)=c24sz;
                clustsizemat(25,2)=c25sz;
                clustsizemat(26,2)=c26sz;
                clustsizemat(27,2)=c27sz;
                clustsizemat(28,2)=c28sz;
                clustsizemat(29,2)=c29sz;
                clustsizemat(30,2)=c30sz;
                clustsizemat(31,2)=c31sz;
                clustsizemat(32,2)=c32sz;
                clustsizemat(33,2)=c33sz;
                clustsizemat(34,2)=c34sz;
                clustsizemat(35,2)=c35sz;
                clustsizemat(36,2)=c36sz;
                clustsizemat(37,2)=c37sz;
                clustsizemat(38,2)=c38sz;
                clustsizemat(39,2)=c39sz;
                clustsizemat(40,2)=c40sz;
                clustsizemat(41,2)=c41sz;
                clustsizemat(42,2)=c42sz;
                clustsizemat(43,2)=c43sz;
                clustsizemat(44,2)=c44sz;
                clustsizemat(45,2)=c45sz;
                clustsizemat(46,2)=c46sz;
                clustsizemat(47,2)=c47sz;
                clustsizemat(48,2)=c48sz;
                clustsizemat(49,2)=c49sz;
                clustsizemat(50,2)=c50sz;
                
                
                for meanindx=2:size(subjasstblock,2)-1
                    c1(c1sz-1,meanindx)=(sum(c1(1:c1sz-2,meanindx)))/c1sz-2;
                    c1(c1sz,meanindx)=std(c1(1:c1sz-2,meanindx));
                    c2(c2sz-1,meanindx)=(sum(c2(1:c2sz-2,meanindx)))/c2sz-2;
                    c2(c2sz,meanindx)=std(c2(1:c2sz-2,meanindx));
                    c3(c3sz-1,meanindx)=(sum(c3(1:c3sz-2,meanindx)))/c3sz-2;
                    c3(c3sz,meanindx)=std(c3(1:c3sz-2,meanindx));
                    c4(c4sz-1,meanindx)=(sum(c4(1:c4sz-2,meanindx)))/c4sz-2;
                    c4(c4sz,meanindx)=std(c4(1:c4sz-2,meanindx));
                    c5(c5sz-1,meanindx)=(sum(c5(1:c5sz-2,meanindx)))/c5sz-2;
                    c5(c5sz,meanindx)=std(c5(1:c5sz-2,meanindx));
                    c6(c6sz-1,meanindx)=(sum(c6(1:c6sz-2,meanindx)))/c6sz-2;
                    c6(c6sz,meanindx)=std(c6(1:c6sz-2,meanindx));
                    c7(c7sz-1,meanindx)=(sum(c7(1:c7sz-2,meanindx)))/c7sz-2;
                    c7(c7sz,meanindx)=std(c7(1:c7sz-2,meanindx));
                    c8(c8sz-1,meanindx)=(sum(c8(1:c8sz-2,meanindx)))/c8sz-2;
                    c8(c8sz,meanindx)=std(c8(1:c8sz-2,meanindx));
                    c9(c9sz-1,meanindx)=(sum(c9(1:c9sz-2,meanindx)))/c9sz-2;
                    c9(c9sz,meanindx)=std(c9(1:c9sz-2,meanindx));
                    c10(c10sz-1,meanindx)=(sum(c10(1:c10sz-2,meanindx)))/c10sz-2;
                    c10(c10sz,meanindx)=std(c10(1:c10sz-2,meanindx));
                    c11(c11sz-1,meanindx)=(sum(c11(1:c11sz-2,meanindx)))/c11sz-2;
                    c11(c11sz,meanindx)=std(c11(1:c11sz-2,meanindx));
                    c12(c12sz-1,meanindx)=(sum(c12(1:c12sz-2,meanindx)))/c12sz-2;
                    c12(c12sz,meanindx)=std(c12(1:c12sz-2,meanindx));
                    c13(c13sz-1,meanindx)=(sum(c13(1:c13sz-2,meanindx)))/c13sz-2;
                    c13(c13sz,meanindx)=std(c13(1:c13sz-2,meanindx));
                    c14(c14sz-1,meanindx)=(sum(c14(1:c14sz-2,meanindx)))/c14sz-2;
                    c14(c14sz,meanindx)=std(c14(1:c14sz-2,meanindx));
                    c15(c15sz-1,meanindx)=(sum(c15(1:c15sz-2,meanindx)))/c15sz-2;
                    c15(c15sz,meanindx)=std(c15(1:c15sz-2,meanindx));
                    c16(c16sz-1,meanindx)=(sum(c16(1:c16sz-2,meanindx)))/c16sz-2;
                    c16(c16sz,meanindx)=std(c16(1:c16sz-2,meanindx));
                    c17(c17sz-1,meanindx)=(sum(c17(1:c17sz-2,meanindx)))/c17sz-2;
                    c17(c17sz,meanindx)=std(c17(1:c17sz-2,meanindx));
                    c18(c18sz-1,meanindx)=(sum(c18(1:c18sz-2,meanindx)))/c18sz-2;
                    c18(c18sz,meanindx)=std(c18(1:c18sz-2,meanindx));
                    c19(c19sz-1,meanindx)=(sum(c19(1:c19sz-2,meanindx)))/c19sz-2;
                    c19(c19sz,meanindx)=std(c19(1:c19sz-2,meanindx));
                    c20(c20sz-1,meanindx)=(sum(c20(1:c20sz-2,meanindx)))/c20sz-2;
                    c20(c20sz,meanindx)=std(c20(1:c20sz-2,meanindx));
                    c21(c21sz-1,meanindx)=(sum(c21(1:c21sz-2,meanindx)))/c21sz-2;
                    c21(c21sz,meanindx)=std(c21(1:c21sz-2,meanindx));
                    c22(c22sz-1,meanindx)=(sum(c22(1:c22sz-2,meanindx)))/c22sz-2;
                    c22(c22sz,meanindx)=std(c22(1:c22sz-2,meanindx));
                    c23(c23sz-1,meanindx)=(sum(c23(1:c23sz-2,meanindx)))/c23sz-2;
                    c23(c23sz,meanindx)=std(c23(1:c23sz-2,meanindx));
                    c24(c24sz-1,meanindx)=(sum(c24(1:c24sz-2,meanindx)))/c24sz-2;
                    c24(c24sz,meanindx)=std(c24(1:c24sz-2,meanindx));
                    c25(c25sz-1,meanindx)=(sum(c25(1:c25sz-2,meanindx)))/c25sz-2;
                    c25(c25sz,meanindx)=std(c25(1:c25sz-2,meanindx));
                    c26(c26sz-1,meanindx)=(sum(c26(1:c26sz-2,meanindx)))/c26sz-2;
                    c26(c26sz,meanindx)=std(c26(1:c26sz-2,meanindx));
                    c27(c27sz-1,meanindx)=(sum(c27(1:c27sz-2,meanindx)))/c27sz-2;
                    c27(c27sz,meanindx)=std(c27(1:c27sz-2,meanindx));
                    c28(c28sz-1,meanindx)=(sum(c28(1:c28sz-2,meanindx)))/c28sz-2;
                    c28(c28sz,meanindx)=std(c28(1:c28sz-2,meanindx));
                    c29(c29sz-1,meanindx)=(sum(c29(1:c29sz-2,meanindx)))/c29sz-2;
                    c29(c29sz,meanindx)=std(c29(1:c29sz-2,meanindx));
                    c30(c30sz-1,meanindx)=(sum(c30(1:c30sz-2,meanindx)))/c30sz-2;
                    c30(c30sz,meanindx)=std(c30(1:c30sz-2,meanindx));
                    c31(c31sz-1,meanindx)=(sum(c31(1:c31sz-2,meanindx)))/c31sz-2;
                    c31(c31sz,meanindx)=std(c31(1:c31sz-2,meanindx));
                    c32(c32sz-1,meanindx)=(sum(c32(1:c32sz-2,meanindx)))/c32sz-2;
                    c32(c32sz,meanindx)=std(c32(1:c32sz-2,meanindx));
                    c33(c33sz-1,meanindx)=(sum(c33(1:c33sz-2,meanindx)))/c33sz-2;
                    c33(c33sz,meanindx)=std(c33(1:c33sz-2,meanindx));
                    c34(c34sz-1,meanindx)=(sum(c34(1:c34sz-2,meanindx)))/c34sz-2;
                    c34(c34sz,meanindx)=std(c34(1:c34sz-2,meanindx));
                    c35(c35sz-1,meanindx)=(sum(c35(1:c35sz-2,meanindx)))/c35sz-2;
                    c35(c35sz,meanindx)=std(c35(1:c35sz-2,meanindx));
                    c36(c36sz-1,meanindx)=(sum(c36(1:c36sz-2,meanindx)))/c36sz-2;
                    c36(c36sz,meanindx)=std(c36(1:c36sz-2,meanindx));
                    c37(c37sz-1,meanindx)=(sum(c37(1:c37sz-2,meanindx)))/c37sz-2;
                    c37(c37sz,meanindx)=std(c37(1:c37sz-2,meanindx));
                    c38(c38sz-1,meanindx)=(sum(c38(1:c38sz-2,meanindx)))/c38sz-2;
                    c38(c38sz,meanindx)=std(c38(1:c38sz-2,meanindx));
                    c39(c39sz-1,meanindx)=(sum(c39(1:c39sz-2,meanindx)))/c39sz-2;
                    c39(c39sz,meanindx)=std(c39(1:c39sz-2,meanindx));
                    c40(c40sz-1,meanindx)=(sum(c40(1:c40sz-2,meanindx)))/c40sz-2;
                    c40(c40sz,meanindx)=std(c40(1:c40sz-2,meanindx));
                    c41(c41sz-1,meanindx)=(sum(c41(1:c41sz-2,meanindx)))/c41sz-2;
                    c41(c41sz,meanindx)=std(c41(1:c41sz-2,meanindx));
                    c42(c42sz-1,meanindx)=(sum(c42(1:c42sz-2,meanindx)))/c42sz-2;
                    c42(c42sz,meanindx)=std(c42(1:c42sz-2,meanindx));
                    c43(c43sz-1,meanindx)=(sum(c43(1:c43sz-2,meanindx)))/c43sz-2;
                    c43(c43sz,meanindx)=std(c43(1:c43sz-2,meanindx));
                    c44(c44sz-1,meanindx)=(sum(c44(1:c44sz-2,meanindx)))/c44sz-2;
                    c44(c44sz,meanindx)=std(c44(1:c44sz-2,meanindx));
                    c45(c45sz-1,meanindx)=(sum(c45(1:c45sz-2,meanindx)))/c45sz-2;
                    c45(c45sz,meanindx)=std(c45(1:c45sz-2,meanindx));
                    c46(c46sz-1,meanindx)=(sum(c46(1:c46sz-2,meanindx)))/c46sz-2;
                    c46(c46sz,meanindx)=std(c46(1:c46sz-2,meanindx));
                    c47(c47sz-1,meanindx)=(sum(c47(1:c47sz-2,meanindx)))/c47sz-2;
                    c47(c47sz,meanindx)=std(c47(1:c47sz-2,meanindx));
                    c48(c48sz-1,meanindx)=(sum(c48(1:c48sz-2,meanindx)))/c48sz-2;
                    c48(c48sz,meanindx)=std(c48(1:c48sz-2,meanindx));
                    c49(c49sz-1,meanindx)=(sum(c49(1:c49sz-2,meanindx)))/c49sz-2;
                    c49(c49sz,meanindx)=std(c49(1:c49sz-2,meanindx));
                    c50(c50sz-1,meanindx)=(sum(c50(1:c50sz-2,meanindx)))/c50sz-2;
                    c50(c50sz,meanindx)=std(c50(1:c50sz-2,meanindx));
                end
                
                colnum=1;
                gdnsmat1cnt=8; %Allows for the opening index column

                for frchckindx=1:13:size(fourierwavs,2)
                    tmpfourierwavs=fourierwavs(:,frchckindx:frchckindx+12);
                    tmpgoodnessmat=zeros(50,6);
                    goodnessmat=zeros(50,7);
                    goodnessmat(:,1)=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50];
                    
                    tmpc1=c1(1:end-2,2:14);
                    cmean=cmeanmat(1,colnum);
                    holdtmp=tmpc1;
                    for sstindx=1:size(tmpc1,2)
                        for sstrow=1:size(tmpc1,1)
                            tmpc1(sstrow,sstindx)=tmpc1(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc1=tmpc1.^2;
                    if size(tmpc1,1)==1
                        sst=sum(tmpc1(1,1:13));
                    else
                        sst=sum(tmpc1);
                    end
                    sst=sum(sst);
                    tmpc1=holdtmp;
                    dfe=(size(tmpc1,1)*size(tmpc1,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc1,2)
                        for sserow=1:size(tmpc1,1)
                            tmpc1(sserow,ssecol)=(tmpc1(sserow,ssecol)-tmpfourierwavs(1,ssecol));
                        end
                    end
                    tmpc1=tmpc1.^2;
                    sse=sum(tmpc1);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc1,1)*size(tmpc1,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(1,2)=rsquare;
                        goodnessmat(1,3)=adjrsquare;
                        goodnessmat(1,4)=sse;
                        goodnessmat(1,5)=rmse;
                        goodnessmat(1,6)=RIC;
                        goodnessmat(1,7)=sst;
                    else
                        tmpgoodnessmat(1,1)=rsquare;
                        tmpgoodnessmat(1,2)=adjrsquare;
                        tmpgoodnessmat(1,3)=sse;
                        tmpgoodnessmat(1,4)=rmse;
                        tmpgoodnessmat(1,5)=RIC;
                        tmpgoodnessmat(1,6)=sst;
                    end
                    
                    tmpc2=c2(1:end-2,2:14);
                    cmean=cmeanmat(2,colnum);
                    holdtmp=tmpc2;
                    for sstindx=1:size(tmpc2,2)
                        for sstrow=1:size(tmpc2,1)
                            tmpc2(sstrow,sstindx)=tmpc2(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc2=tmpc2.^2;
                    if size(tmpc2,1)==1
                        sst=sum(tmpc2(1,1:13));
                    else
                        sst=sum(tmpc2);
                    end
                    sst=sum(sst);
                    tmpc2=holdtmp;
                    dfe=(size(tmpc2,1)*size(tmpc2,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc2,2)
                        for sserow=1:size(tmpc2,1)
                            tmpc2(sserow,ssecol)=(tmpc2(sserow,ssecol)-tmpfourierwavs(2,ssecol));
                        end
                    end
                    tmpc2=tmpc2.^2;
                    sse=sum(tmpc2);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    adjrsquare=1-(1-rsquare)*((size(tmpc2,1)*size(tmpc2,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(2,2)=rsquare;
                        goodnessmat(2,3)=adjrsquare;
                        goodnessmat(2,4)=sse;
                        goodnessmat(2,5)=rmse;
                        goodnessmat(2,6)=RIC;
                        goodnessmat(2,7)=sst;
                    else
                        tmpgoodnessmat(2,1)=rsquare;
                        tmpgoodnessmat(2,2)=adjrsquare;
                        tmpgoodnessmat(2,3)=sse;
                        tmpgoodnessmat(2,4)=rmse;
                        tmpgoodnessmat(2,5)=RIC;
                        tmpgoodnessmat(2,6)=sst;
                    end
                    
                    tmpc3=c3(1:end-2,2:14);
                    cmean=cmeanmat(3,colnum);
                    holdtmp=tmpc3;
                    for sstindx=1:size(tmpc3,2)
                        for sstrow=1:size(tmpc3,1)
                            tmpc3(sstrow,sstindx)=tmpc3(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc3=tmpc3.^2;
                    if size(tmpc3,1)==1
                        sst=sum(tmpc3(1,1:13));
                    else
                        sst=sum(tmpc3);
                    end
                    sst=sum(sst);
                    tmpc3=holdtmp;
                    dfe=(size(tmpc3,1)*size(tmpc3,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc3,2)
                        for sserow=1:size(tmpc3,1)
                            tmpc3(sserow,ssecol)=(tmpc3(sserow,ssecol)-tmpfourierwavs(3,ssecol));
                        end
                    end
                    tmpc3=tmpc3.^2;
                    sse=sum(tmpc3);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc3,1)*size(tmpc3,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(3,2)=rsquare;
                        goodnessmat(3,3)=adjrsquare;
                        goodnessmat(3,4)=sse;
                        goodnessmat(3,5)=rmse;
                        goodnessmat(3,6)=RIC;
                        goodnessmat(3,7)=sst;
                    else
                        tmpgoodnessmat(3,1)=rsquare;
                        tmpgoodnessmat(3,2)=adjrsquare;
                        tmpgoodnessmat(3,3)=sse;
                        tmpgoodnessmat(3,4)=rmse;
                        tmpgoodnessmat(3,5)=RIC;
                        tmpgoodnessmat(3,6)=sst;
                    end
                    
                    tmpc4=c4(1:end-2,2:14);
                    cmean=cmeanmat(4,colnum);
                    holdtmp=tmpc4;
                    for sstindx=1:size(tmpc4,2)
                        for sstrow=1:size(tmpc4,1)
                            tmpc4(sstrow,sstindx)=tmpc4(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc4=tmpc4.^2;
                    if size(tmpc4,1)==1
                        sst=sum(tmpc4(1,1:13));
                    else
                        sst=sum(tmpc4);
                    end
                    sst=sum(sst);
                    tmpc4=holdtmp;
                    dfe=(size(tmpc4,1)*size(tmpc4,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc4,2)
                        for sserow=1:size(tmpc4,1)
                            tmpc4(sserow,ssecol)=(tmpc4(sserow,ssecol)-tmpfourierwavs(4,ssecol));
                        end
                    end
                    tmpc4=tmpc4.^2;
                    sse=sum(tmpc4);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc4,1)*size(tmpc4,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(4,2)=rsquare;
                        goodnessmat(4,3)=adjrsquare;
                        goodnessmat(4,4)=sse;
                        goodnessmat(4,5)=rmse;
                        goodnessmat(4,6)=RIC;
                        goodnessmat(4,7)=sst;
                    else
                        tmpgoodnessmat(4,1)=rsquare;
                        tmpgoodnessmat(4,2)=adjrsquare;
                        tmpgoodnessmat(4,3)=sse;
                        tmpgoodnessmat(4,4)=rmse;
                        tmpgoodnessmat(4,5)=RIC;
                        tmpgoodnessmat(4,6)=sst;
                    end
                    
                    tmpc5=c5(1:end-2,2:14);
                    cmean=cmeanmat(5,colnum);
                    holdtmp=tmpc5;
                    for sstindx=1:size(tmpc5,2)
                        for sstrow=1:size(tmpc5,1)
                            tmpc5(sstrow,sstindx)=tmpc5(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc5=tmpc5.^2;
                    if size(tmpc5,1)==1
                        sst=sum(tmpc5(1,1:13));
                    else
                        sst=sum(tmpc5);
                    end
                    sst=sum(sst);
                    tmpc5=holdtmp;
                    dfe=(size(tmpc5,1)*size(tmpc5,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc5,2)
                        for sserow=1:size(tmpc5,1)
                            tmpc5(sserow,ssecol)=(tmpc5(sserow,ssecol)-tmpfourierwavs(5,ssecol));
                        end
                    end
                    tmpc5=tmpc5.^2;
                    sse=sum(tmpc5);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc5,1)*size(tmpc5,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(5,2)=rsquare;
                        goodnessmat(5,3)=adjrsquare;
                        goodnessmat(5,4)=sse;
                        goodnessmat(5,5)=rmse;
                        goodnessmat(5,6)=RIC;
                        goodnessmat(5,7)=sst;
                    else
                        tmpgoodnessmat(5,1)=rsquare;
                        tmpgoodnessmat(5,2)=adjrsquare;
                        tmpgoodnessmat(5,3)=sse;
                        tmpgoodnessmat(5,4)=rmse;
                        tmpgoodnessmat(5,5)=RIC;
                        tmpgoodnessmat(5,6)=sst;
                    end
                    
                    tmpc6=c6(1:end-2,2:14);
                    cmean=cmeanmat(6,colnum);
                    holdtmp=tmpc6;
                    for sstindx=1:size(tmpc6,2)
                        for sstrow=1:size(tmpc6,1)
                            tmpc6(sstrow,sstindx)=tmpc6(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc6=tmpc6.^2;
                    if size(tmpc6,1)==1
                        sst=sum(tmpc6(1,1:13));
                    else
                        sst=sum(tmpc6);
                    end
                    sst=sum(sst);
                    tmpc6=holdtmp;
                    dfe=(size(tmpc6,1)*size(tmpc6,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc6,2)
                        for sserow=1:size(tmpc6,1)
                            tmpc6(sserow,ssecol)=(tmpc6(sserow,ssecol)-tmpfourierwavs(6,ssecol));
                        end
                    end
                    tmpc6=tmpc6.^2;
                    sse=sum(tmpc6);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc6,1)*size(tmpc6,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(6,2)=rsquare;
                        goodnessmat(6,3)=adjrsquare;
                        goodnessmat(6,4)=sse;
                        goodnessmat(6,5)=rmse;
                        goodnessmat(6,6)=RIC;
                        goodnessmat(6,7)=sst;
                    else
                        tmpgoodnessmat(6,1)=rsquare;
                        tmpgoodnessmat(6,2)=adjrsquare;
                        tmpgoodnessmat(6,3)=sse;
                        tmpgoodnessmat(6,4)=rmse;
                        tmpgoodnessmat(6,5)=RIC;
                        tmpgoodnessmat(6,6)=sst;
                    end
                    
                    tmpc7=c7(1:end-2,2:14);
                    cmean=cmeanmat(7,colnum);
                    holdtmp=tmpc7;
                    for sstindx=1:size(tmpc7,2)
                        for sstrow=1:size(tmpc7,1)
                            tmpc7(sstrow,sstindx)=tmpc7(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc7=tmpc7.^2;
                    if size(tmpc7,1)==1
                        sst=sum(tmpc7(1,1:13));
                    else
                        sst=sum(tmpc7);
                    end
                    sst=sum(sst);
                    tmpc7=holdtmp;
                    dfe=(size(tmpc7,1)*size(tmpc7,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc7,2)
                        for sserow=1:size(tmpc7,1)
                            tmpc7(sserow,ssecol)=(tmpc7(sserow,ssecol)-tmpfourierwavs(7,ssecol));
                        end
                    end
                    tmpc7=tmpc7.^2;
                    sse=sum(tmpc7);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc7,1)*size(tmpc7,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(7,2)=rsquare;
                        goodnessmat(7,3)=adjrsquare;
                        goodnessmat(7,4)=sse;
                        goodnessmat(7,5)=rmse;
                        goodnessmat(7,6)=RIC;
                        goodnessmat(7,7)=sst;
                    else
                        tmpgoodnessmat(7,1)=rsquare;
                        tmpgoodnessmat(7,2)=adjrsquare;
                        tmpgoodnessmat(7,3)=sse;
                        tmpgoodnessmat(7,4)=rmse;
                        tmpgoodnessmat(7,5)=RIC;
                        tmpgoodnessmat(7,6)=sst;
                    end
                    
                    tmpc8=c8(1:end-2,2:14);
                    cmean=cmeanmat(8,colnum);
                    holdtmp=tmpc8;
                    for sstindx=1:size(tmpc8,2)
                        for sstrow=1:size(tmpc8,1)
                            tmpc8(sstrow,sstindx)=tmpc8(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc8=tmpc8.^2;
                    if size(tmpc8,1)==1
                        sst=sum(tmpc8(1,1:13));
                    else
                        sst=sum(tmpc8);
                    end
                    sst=sum(sst);
                    tmpc8=holdtmp;
                    dfe=(size(tmpc8,1)*size(tmpc8,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc8,2)
                        for sserow=1:size(tmpc8,1)
                            tmpc8(sserow,ssecol)=(tmpc8(sserow,ssecol)-tmpfourierwavs(8,ssecol));
                        end
                    end
                    tmpc8=tmpc8.^2;
                    sse=sum(tmpc8);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc8,1)*size(tmpc8,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(8,2)=rsquare;
                        goodnessmat(8,3)=adjrsquare;
                        goodnessmat(8,4)=sse;
                        goodnessmat(8,5)=rmse;
                        goodnessmat(8,6)=RIC;
                        goodnessmat(8,7)=sst;
                    else
                        tmpgoodnessmat(8,1)=rsquare;
                        tmpgoodnessmat(8,2)=adjrsquare;
                        tmpgoodnessmat(8,3)=sse;
                        tmpgoodnessmat(8,4)=rmse;
                        tmpgoodnessmat(8,5)=RIC;
                        tmpgoodnessmat(8,6)=sst;
                    end
                    
                    tmpc9=c9(1:end-2,2:14);
                    cmean=cmeanmat(9,colnum);
                    holdtmp=tmpc9;
                    for sstindx=1:size(tmpc9,2)
                        for sstrow=1:size(tmpc9,1)
                            tmpc9(sstrow,sstindx)=tmpc9(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc9=tmpc9.^2;
                    if size(tmpc9,1)==1
                        sst=sum(tmpc9(1,1:13));
                    else
                        sst=sum(tmpc9);
                    end
                    sst=sum(sst);
                    tmpc9=holdtmp;
                    dfe=(size(tmpc9,1)*size(tmpc9,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc9,2)
                        for sserow=1:size(tmpc9,1)
                            tmpc9(sserow,ssecol)=(tmpc9(sserow,ssecol)-tmpfourierwavs(9,ssecol));
                        end
                    end
                    tmpc9=tmpc9.^2;
                    sse=sum(tmpc9);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc9,1)*size(tmpc9,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(9,2)=rsquare;
                        goodnessmat(9,3)=adjrsquare;
                        goodnessmat(9,4)=sse;
                        goodnessmat(9,5)=rmse;
                        goodnessmat(9,6)=RIC;
                        goodnessmat(9,7)=sst;
                    else
                        tmpgoodnessmat(9,1)=rsquare;
                        tmpgoodnessmat(9,2)=adjrsquare;
                        tmpgoodnessmat(9,3)=sse;
                        tmpgoodnessmat(9,4)=rmse;
                        tmpgoodnessmat(9,5)=RIC;
                        tmpgoodnessmat(9,6)=sst;
                    end
                    
                    tmpc10=c10(1:end-2,2:14);
                    cmean=cmeanmat(10,colnum);
                    holdtmp=tmpc10;
                    for sstindx=1:size(tmpc10,2)
                        for sstrow=1:size(tmpc10,1)
                            tmpc10(sstrow,sstindx)=tmpc10(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc10=tmpc10.^2;
                    if size(tmpc10,1)==1
                        sst=sum(tmpc10(1,1:13));
                    else
                        sst=sum(tmpc10);
                    end
                    sst=sum(sst);
                    tmpc10=holdtmp;
                    dfe=(size(tmpc10,1)*size(tmpc10,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc10,2)
                        for sserow=1:size(tmpc10,1)
                            tmpc10(sserow,ssecol)=(tmpc10(sserow,ssecol)-tmpfourierwavs(10,ssecol));
                        end
                    end
                    tmpc10=tmpc10.^2;
                    sse=sum(tmpc10);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc10,1)*size(tmpc10,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(10,2)=rsquare;
                        goodnessmat(10,3)=adjrsquare;
                        goodnessmat(10,4)=sse;
                        goodnessmat(10,5)=rmse;
                        goodnessmat(10,6)=RIC;
                        goodnessmat(10,7)=sst;
                    else
                        tmpgoodnessmat(10,1)=rsquare;
                        tmpgoodnessmat(10,2)=adjrsquare;
                        tmpgoodnessmat(10,3)=sse;
                        tmpgoodnessmat(10,4)=rmse;
                        tmpgoodnessmat(10,5)=RIC;
                        tmpgoodnessmat(10,6)=sst;
                    end
                    
                    tmpc11=c11(1:end-2,2:14);
                    cmean=cmeanmat(11,colnum);
                    holdtmp=tmpc11;
                    for sstindx=1:size(tmpc11,2)
                        for sstrow=1:size(tmpc11,1)
                            tmpc11(sstrow,sstindx)=tmpc11(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc11=tmpc11.^2;
                    if size(tmpc11,1)==1
                        sst=sum(tmpc11(1,1:13));
                    else
                        sst=sum(tmpc11);
                    end
                    sst=sum(sst);
                    tmpc11=holdtmp;
                    dfe=(size(tmpc11,1)*size(tmpc11,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc11,2)
                        for sserow=1:size(tmpc11,1)
                            tmpc11(sserow,ssecol)=(tmpc11(sserow,ssecol)-tmpfourierwavs(11,ssecol));
                        end
                    end
                    tmpc11=tmpc11.^2;
                    sse=sum(tmpc11);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc11,1)*size(tmpc11,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(11,2)=rsquare;
                        goodnessmat(11,3)=adjrsquare;
                        goodnessmat(11,4)=sse;
                        goodnessmat(11,5)=rmse;
                        goodnessmat(11,6)=RIC;
                        goodnessmat(11,7)=sst;
                    else
                        tmpgoodnessmat(11,1)=rsquare;
                        tmpgoodnessmat(11,2)=adjrsquare;
                        tmpgoodnessmat(11,3)=sse;
                        tmpgoodnessmat(11,4)=rmse;
                        tmpgoodnessmat(11,5)=RIC;
                        tmpgoodnessmat(11,6)=sst;
                    end
                    
                    tmpc12=c12(1:end-2,2:14);
                    cmean=cmeanmat(12,colnum);
                    holdtmp=tmpc12;
                    for sstindx=1:size(tmpc12,2)
                        for sstrow=1:size(tmpc12,1)
                            tmpc12(sstrow,sstindx)=tmpc12(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc12=tmpc12.^2;
                    if size(tmpc12,1)==1
                        sst=sum(tmpc12(1,1:13));
                    else
                        sst=sum(tmpc12);
                    end
                    sst=sum(sst);
                    tmpc12=holdtmp;
                    dfe=(size(tmpc12,1)*size(tmpc12,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc12,2)
                        for sserow=1:size(tmpc12,1)
                            tmpc12(sserow,ssecol)=(tmpc12(sserow,ssecol)-tmpfourierwavs(12,ssecol));
                        end
                    end
                    tmpc12=tmpc12.^2;
                    sse=sum(tmpc12);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc12,1)*size(tmpc12,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(12,2)=rsquare;
                        goodnessmat(12,3)=adjrsquare;
                        goodnessmat(12,4)=sse;
                        goodnessmat(12,5)=rmse;
                        goodnessmat(12,6)=RIC;
                        goodnessmat(12,7)=sst;
                    else
                        tmpgoodnessmat(12,1)=rsquare;
                        tmpgoodnessmat(12,2)=adjrsquare;
                        tmpgoodnessmat(12,3)=sse;
                        tmpgoodnessmat(12,4)=rmse;
                        tmpgoodnessmat(12,5)=RIC;
                        tmpgoodnessmat(12,6)=sst;
                    end
                    
                    tmpc13=c13(1:end-2,2:14);
                    cmean=cmeanmat(13,colnum);
                    holdtmp=tmpc13;
                    for sstindx=1:size(tmpc13,2)
                        for sstrow=1:size(tmpc13,1)
                            tmpc13(sstrow,sstindx)=tmpc13(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc13=tmpc13.^2;
                    if size(tmpc13,1)==1
                        sst=sum(tmpc13(1,1:13));
                    else
                        sst=sum(tmpc13);
                    end
                    sst=sum(sst);
                    tmpc13=holdtmp;
                    dfe=(size(tmpc13,1)*size(tmpc13,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc13,2)
                        for sserow=1:size(tmpc13,1)
                            tmpc13(sserow,ssecol)=(tmpc13(sserow,ssecol)-tmpfourierwavs(13,ssecol));
                        end
                    end
                    tmpc13=tmpc13.^2;
                    sse=sum(tmpc13);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc13,1)*size(tmpc13,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(13,2)=rsquare;
                        goodnessmat(13,3)=adjrsquare;
                        goodnessmat(13,4)=sse;
                        goodnessmat(13,5)=rmse;
                        goodnessmat(13,6)=RIC;
                        goodnessmat(13,7)=sst;
                    else
                        tmpgoodnessmat(13,1)=rsquare;
                        tmpgoodnessmat(13,2)=adjrsquare;
                        tmpgoodnessmat(13,3)=sse;
                        tmpgoodnessmat(13,4)=rmse;
                        tmpgoodnessmat(13,5)=RIC;
                        tmpgoodnessmat(13,6)=sst;
                    end
                    
                    tmpc14=c14(1:end-2,2:14);
                    cmean=cmeanmat(14,colnum);
                    holdtmp=tmpc14;
                    for sstindx=1:size(tmpc14,2)
                        for sstrow=1:size(tmpc14,1)
                            tmpc14(sstrow,sstindx)=tmpc14(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc14=tmpc14.^2;
                    if size(tmpc14,1)==1
                        sst=sum(tmpc14(1,1:13));
                    else
                        sst=sum(tmpc14);
                    end
                    sst=sum(sst);
                    tmpc14=holdtmp;
                    dfe=(size(tmpc14,1)*size(tmpc14,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc14,2)
                        for sserow=1:size(tmpc14,1)
                            tmpc14(sserow,ssecol)=(tmpc14(sserow,ssecol)-tmpfourierwavs(14,ssecol));
                        end
                    end
                    tmpc14=tmpc14.^2;
                    sse=sum(tmpc14);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc14,1)*size(tmpc14,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(14,2)=rsquare;
                        goodnessmat(14,3)=adjrsquare;
                        goodnessmat(14,4)=sse;
                        goodnessmat(14,5)=rmse;
                        goodnessmat(14,6)=RIC;
                        goodnessmat(14,7)=sst;
                    else
                        tmpgoodnessmat(14,1)=rsquare;
                        tmpgoodnessmat(14,2)=adjrsquare;
                        tmpgoodnessmat(14,3)=sse;
                        tmpgoodnessmat(14,4)=rmse;
                        tmpgoodnessmat(14,5)=RIC;
                        tmpgoodnessmat(14,6)=sst;
                    end
                    
                    tmpc15=c15(1:end-2,2:14);
                    cmean=cmeanmat(15,colnum);
                    holdtmp=tmpc15;
                    for sstindx=1:size(tmpc15,2)
                        for sstrow=1:size(tmpc15,1)
                            tmpc15(sstrow,sstindx)=tmpc15(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc15=tmpc15.^2;
                    if size(tmpc15,1)==1
                        sst=sum(tmpc15(1,1:13));
                    else
                        sst=sum(tmpc15);
                    end
                    sst=sum(sst);
                    tmpc15=holdtmp;
                    dfe=(size(tmpc15,1)*size(tmpc15,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc15,2)
                        for sserow=1:size(tmpc15,1)
                            tmpc15(sserow,ssecol)=(tmpc15(sserow,ssecol)-tmpfourierwavs(15,ssecol));
                        end
                    end
                    tmpc15=tmpc15.^2;
                    sse=sum(tmpc15);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc15,1)*size(tmpc15,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(15,2)=rsquare;
                        goodnessmat(15,3)=adjrsquare;
                        goodnessmat(15,4)=sse;
                        goodnessmat(15,5)=rmse;
                        goodnessmat(15,6)=RIC;
                        goodnessmat(15,7)=sst;
                    else
                        tmpgoodnessmat(15,1)=rsquare;
                        tmpgoodnessmat(15,2)=adjrsquare;
                        tmpgoodnessmat(15,3)=sse;
                        tmpgoodnessmat(15,4)=rmse;
                        tmpgoodnessmat(15,5)=RIC;
                        tmpgoodnessmat(15,6)=sst;
                    end
                    
                    tmpc16=c16(1:end-2,2:14);
                    cmean=cmeanmat(16,colnum);
                    holdtmp=tmpc16;
                    for sstindx=1:size(tmpc16,2)
                        for sstrow=1:size(tmpc16,1)
                            tmpc16(sstrow,sstindx)=tmpc16(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc16=tmpc16.^2;
                    if size(tmpc16,1)==1
                        sst=sum(tmpc16(1,1:13));
                    else
                        sst=sum(tmpc16);
                    end
                    sst=sum(sst);
                    tmpc16=holdtmp;
                    dfe=(size(tmpc16,1)*size(tmpc16,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc16,2)
                        for sserow=1:size(tmpc16,1)
                            tmpc16(sserow,ssecol)=(tmpc16(sserow,ssecol)-tmpfourierwavs(16,ssecol));
                        end
                    end
                    tmpc16=tmpc16.^2;
                    sse=sum(tmpc16);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc16,1)*size(tmpc16,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(16,2)=rsquare;
                        goodnessmat(16,3)=adjrsquare;
                        goodnessmat(16,4)=sse;
                        goodnessmat(16,5)=rmse;
                        goodnessmat(16,6)=RIC;
                        goodnessmat(16,7)=sst;
                    else
                        tmpgoodnessmat(16,1)=rsquare;
                        tmpgoodnessmat(16,2)=adjrsquare;
                        tmpgoodnessmat(16,3)=sse;
                        tmpgoodnessmat(16,4)=rmse;
                        tmpgoodnessmat(16,5)=RIC;
                        tmpgoodnessmat(16,6)=sst;
                    end
                    
                    tmpc17=c17(1:end-2,2:14);
                    cmean=cmeanmat(17,colnum);
                    holdtmp=tmpc17;
                    for sstindx=1:size(tmpc17,2)
                        for sstrow=1:size(tmpc17,1)
                            tmpc17(sstrow,sstindx)=tmpc17(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc17=tmpc17.^2;
                    if size(tmpc17,1)==1
                        sst=sum(tmpc17(1,1:13));
                    else
                        sst=sum(tmpc17);
                    end
                    sst=sum(sst);
                    tmpc17=holdtmp;
                    dfe=(size(tmpc17,1)*size(tmpc17,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc17,2)
                        for sserow=1:size(tmpc17,1)
                            tmpc17(sserow,ssecol)=(tmpc17(sserow,ssecol)-tmpfourierwavs(17,ssecol));
                        end
                    end
                    tmpc17=tmpc17.^2;
                    sse=sum(tmpc17);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc17,1)*size(tmpc17,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(17,2)=rsquare;
                        goodnessmat(17,3)=adjrsquare;
                        goodnessmat(17,4)=sse;
                        goodnessmat(17,5)=rmse;
                        goodnessmat(17,6)=RIC;
                        goodnessmat(17,7)=sst;
                    else
                        tmpgoodnessmat(17,1)=rsquare;
                        tmpgoodnessmat(17,2)=adjrsquare;
                        tmpgoodnessmat(17,3)=sse;
                        tmpgoodnessmat(17,4)=rmse;
                        tmpgoodnessmat(17,5)=RIC;
                        tmpgoodnessmat(17,6)=sst;
                    end
                    
                    tmpc18=c18(1:end-2,2:14);
                    cmean=cmeanmat(18,colnum);
                    holdtmp=tmpc18;
                    for sstindx=1:size(tmpc18,2)
                        for sstrow=1:size(tmpc18,1)
                            tmpc18(sstrow,sstindx)=tmpc18(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc18=tmpc18.^2;
                    if size(tmpc18,1)==1
                        sst=sum(tmpc18(1,1:13));
                    else
                        sst=sum(tmpc18);
                    end
                    sst=sum(sst);
                    tmpc18=holdtmp;
                    dfe=(size(tmpc18,1)*size(tmpc18,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc18,2)
                        for sserow=1:size(tmpc18,1)
                            tmpc18(sserow,ssecol)=(tmpc18(sserow,ssecol)-tmpfourierwavs(18,ssecol));
                        end
                    end
                    tmpc18=tmpc18.^2;
                    sse=sum(tmpc18);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc18,1)*size(tmpc18,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(18,2)=rsquare;
                        goodnessmat(18,3)=adjrsquare;
                        goodnessmat(18,4)=sse;
                        goodnessmat(18,5)=rmse;
                        goodnessmat(18,6)=RIC;
                        goodnessmat(18,7)=sst;
                    else
                        tmpgoodnessmat(18,1)=rsquare;
                        tmpgoodnessmat(18,2)=adjrsquare;
                        tmpgoodnessmat(18,3)=sse;
                        tmpgoodnessmat(18,4)=rmse;
                        tmpgoodnessmat(18,5)=RIC;
                        tmpgoodnessmat(18,6)=sst;
                    end
                    
                    tmpc19=c19(1:end-2,2:14);
                    cmean=cmeanmat(19,colnum);
                    holdtmp=tmpc19;
                    for sstindx=1:size(tmpc19,2)
                        for sstrow=1:size(tmpc19,1)
                            tmpc19(sstrow,sstindx)=tmpc19(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc19=tmpc19.^2;
                    if size(tmpc19,1)==1
                        sst=sum(tmpc19(1,1:13));
                    else
                        sst=sum(tmpc19);
                    end
                    sst=sum(sst);
                    tmpc19=holdtmp;
                    dfe=(size(tmpc19,1)*size(tmpc19,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc19,2)
                        for sserow=1:size(tmpc19,1)
                            tmpc19(sserow,ssecol)=(tmpc19(sserow,ssecol)-tmpfourierwavs(19,ssecol));
                        end
                    end
                    tmpc19=tmpc19.^2;
                    sse=sum(tmpc19);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc19,1)*size(tmpc19,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(19,2)=rsquare;
                        goodnessmat(19,3)=adjrsquare;
                        goodnessmat(19,4)=sse;
                        goodnessmat(19,5)=rmse;
                        goodnessmat(19,6)=RIC;
                        goodnessmat(19,7)=sst;
                    else
                        tmpgoodnessmat(19,1)=rsquare;
                        tmpgoodnessmat(19,2)=adjrsquare;
                        tmpgoodnessmat(19,3)=sse;
                        tmpgoodnessmat(19,4)=rmse;
                        tmpgoodnessmat(19,5)=RIC;
                        tmpgoodnessmat(19,6)=sst;
                    end
                    
                    tmpc20=c20(1:end-2,2:14);
                    cmean=cmeanmat(20,colnum);
                    holdtmp=tmpc20;
                    for sstindx=1:size(tmpc20,2)
                        for sstrow=1:size(tmpc20,1)
                            tmpc20(sstrow,sstindx)=tmpc20(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc20=tmpc20.^2;
                    if size(tmpc20,1)==1
                        sst=sum(tmpc20(1,1:13));
                    else
                        sst=sum(tmpc20);
                    end
                    sst=sum(sst);
                    tmpc20=holdtmp;
                    dfe=(size(tmpc20,1)*size(tmpc20,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc20,2)
                        for sserow=1:size(tmpc20,1)
                            tmpc20(sserow,ssecol)=(tmpc20(sserow,ssecol)-tmpfourierwavs(20,ssecol));
                        end
                    end
                    tmpc20=tmpc20.^2;
                    sse=sum(tmpc20);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc20,1)*size(tmpc20,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(20,2)=rsquare;
                        goodnessmat(20,3)=adjrsquare;
                        goodnessmat(20,4)=sse;
                        goodnessmat(20,5)=rmse;
                        goodnessmat(20,6)=RIC;
                        goodnessmat(20,7)=sst;
                    else
                        tmpgoodnessmat(20,1)=rsquare;
                        tmpgoodnessmat(20,2)=adjrsquare;
                        tmpgoodnessmat(20,3)=sse;
                        tmpgoodnessmat(20,4)=rmse;
                        tmpgoodnessmat(20,5)=RIC;
                        tmpgoodnessmat(20,6)=sst;
                    end
                    
                    tmpc21=c21(1:end-2,2:14);
                    cmean=cmeanmat(21,colnum);
                    holdtmp=tmpc21;
                    for sstindx=1:size(tmpc21,2)
                        for sstrow=1:size(tmpc21,1)
                            tmpc21(sstrow,sstindx)=tmpc21(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc21=tmpc21.^2;
                    if size(tmpc21,1)==1
                        sst=sum(tmpc21(1,1:13));
                    else
                        sst=sum(tmpc21);
                    end
                    sst=sum(sst);
                    tmpc21=holdtmp;
                    dfe=(size(tmpc21,1)*size(tmpc21,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc21,2)
                        for sserow=1:size(tmpc21,1)
                            tmpc21(sserow,ssecol)=(tmpc21(sserow,ssecol)-tmpfourierwavs(21,ssecol));
                        end
                    end
                    tmpc21=tmpc21.^2;
                    sse=sum(tmpc21);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc21,1)*size(tmpc21,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(21,2)=rsquare;
                        goodnessmat(21,3)=adjrsquare;
                        goodnessmat(21,4)=sse;
                        goodnessmat(21,5)=rmse;
                        goodnessmat(21,6)=RIC;
                        goodnessmat(21,7)=sst;
                    else
                        tmpgoodnessmat(21,1)=rsquare;
                        tmpgoodnessmat(21,2)=adjrsquare;
                        tmpgoodnessmat(21,3)=sse;
                        tmpgoodnessmat(21,4)=rmse;
                        tmpgoodnessmat(21,5)=RIC;
                        tmpgoodnessmat(21,6)=sst;
                    end
                    
                    tmpc22=c22(1:end-2,2:14);
                    cmean=cmeanmat(22,colnum);
                    holdtmp=tmpc22;
                    for sstindx=1:size(tmpc22,2)
                        for sstrow=1:size(tmpc22,1)
                            tmpc22(sstrow,sstindx)=tmpc22(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc22=tmpc22.^2;
                    if size(tmpc22,1)==1
                        sst=sum(tmpc22(1,1:13));
                    else
                        sst=sum(tmpc22);
                    end
                    sst=sum(sst);
                    tmpc22=holdtmp;
                    dfe=(size(tmpc22,1)*size(tmpc22,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc22,2)
                        for sserow=1:size(tmpc22,1)
                            tmpc22(sserow,ssecol)=(tmpc22(sserow,ssecol)-tmpfourierwavs(22,ssecol));
                        end
                    end
                    tmpc22=tmpc22.^2;
                    sse=sum(tmpc22);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc22,1)*size(tmpc22,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(22,2)=rsquare;
                        goodnessmat(22,3)=adjrsquare;
                        goodnessmat(22,4)=sse;
                        goodnessmat(22,5)=rmse;
                        goodnessmat(22,6)=RIC;
                        goodnessmat(22,7)=sst;
                    else
                        tmpgoodnessmat(22,1)=rsquare;
                        tmpgoodnessmat(22,2)=adjrsquare;
                        tmpgoodnessmat(22,3)=sse;
                        tmpgoodnessmat(22,4)=rmse;
                        tmpgoodnessmat(22,5)=RIC;
                        tmpgoodnessmat(22,6)=sst;
                    end
                    
                    tmpc23=c23(1:end-2,2:14);
                    cmean=cmeanmat(23,colnum);
                    holdtmp=tmpc23;
                    for sstindx=1:size(tmpc23,2)
                        for sstrow=1:size(tmpc23,1)
                            tmpc23(sstrow,sstindx)=tmpc23(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc23=tmpc23.^2;
                    if size(tmpc23,1)==1
                        sst=sum(tmpc23(1,1:13));
                    else
                        sst=sum(tmpc23);
                    end
                    sst=sum(sst);
                    tmpc23=holdtmp;
                    dfe=(size(tmpc23,1)*size(tmpc23,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc23,2)
                        for sserow=1:size(tmpc23,1)
                            tmpc23(sserow,ssecol)=(tmpc23(sserow,ssecol)-tmpfourierwavs(23,ssecol));
                        end
                    end
                    tmpc23=tmpc23.^2;
                    sse=sum(tmpc23);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc23,1)*size(tmpc23,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(23,2)=rsquare;
                        goodnessmat(23,3)=adjrsquare;
                        goodnessmat(23,4)=sse;
                        goodnessmat(23,5)=rmse;
                        goodnessmat(23,6)=RIC;
                        goodnessmat(23,7)=sst;
                    else
                        tmpgoodnessmat(23,1)=rsquare;
                        tmpgoodnessmat(23,2)=adjrsquare;
                        tmpgoodnessmat(23,3)=sse;
                        tmpgoodnessmat(23,4)=rmse;
                        tmpgoodnessmat(23,5)=RIC;
                        tmpgoodnessmat(23,6)=sst;
                    end
                    
                    tmpc24=c24(1:end-2,2:14);
                    cmean=cmeanmat(24,colnum);
                    holdtmp=tmpc24;
                    for sstindx=1:size(tmpc24,2)
                        for sstrow=1:size(tmpc24,1)
                            tmpc24(sstrow,sstindx)=tmpc24(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc24=tmpc24.^2;
                    if size(tmpc24,1)==1
                        sst=sum(tmpc24(1,1:13));
                    else
                        sst=sum(tmpc24);
                    end
                    sst=sum(sst);
                    tmpc24=holdtmp;
                    dfe=(size(tmpc24,1)*size(tmpc24,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc24,2)
                        for sserow=1:size(tmpc24,1)
                            tmpc24(sserow,ssecol)=(tmpc24(sserow,ssecol)-tmpfourierwavs(24,ssecol));
                        end
                    end
                    tmpc24=tmpc24.^2;
                    sse=sum(tmpc24);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc24,1)*size(tmpc24,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(24,2)=rsquare;
                        goodnessmat(24,3)=adjrsquare;
                        goodnessmat(24,4)=sse;
                        goodnessmat(24,5)=rmse;
                        goodnessmat(24,6)=RIC;
                        goodnessmat(24,7)=sst;
                    else
                        tmpgoodnessmat(24,1)=rsquare;
                        tmpgoodnessmat(24,2)=adjrsquare;
                        tmpgoodnessmat(24,3)=sse;
                        tmpgoodnessmat(24,4)=rmse;
                        tmpgoodnessmat(24,5)=RIC;
                        tmpgoodnessmat(24,6)=sst;
                    end
                    
                    tmpc25=c25(1:end-2,2:14);
                    cmean=cmeanmat(25,colnum);
                    holdtmp=tmpc25;
                    for sstindx=1:size(tmpc25,2)
                        for sstrow=1:size(tmpc25,1)
                            tmpc25(sstrow,sstindx)=tmpc25(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc25=tmpc25.^2;
                    if size(tmpc25,1)==1
                        sst=sum(tmpc25(1,1:13));
                    else
                        sst=sum(tmpc25);
                    end
                    sst=sum(sst);
                    tmpc25=holdtmp;
                    dfe=(size(tmpc25,1)*size(tmpc25,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc25,2)
                        for sserow=1:size(tmpc25,1)
                            tmpc25(sserow,ssecol)=(tmpc25(sserow,ssecol)-tmpfourierwavs(25,ssecol));
                        end
                    end
                    tmpc25=tmpc25.^2;
                    sse=sum(tmpc25);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc25,1)*size(tmpc25,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(25,2)=rsquare;
                        goodnessmat(25,3)=adjrsquare;
                        goodnessmat(25,4)=sse;
                        goodnessmat(25,5)=rmse;
                        goodnessmat(25,6)=RIC;
                        goodnessmat(25,7)=sst;
                    else
                        tmpgoodnessmat(25,1)=rsquare;
                        tmpgoodnessmat(25,2)=adjrsquare;
                        tmpgoodnessmat(25,3)=sse;
                        tmpgoodnessmat(25,4)=rmse;
                        tmpgoodnessmat(25,5)=RIC;
                        tmpgoodnessmat(25,6)=sst;
                    end
                    
                    tmpc26=c26(1:end-2,2:14);
                    cmean=cmeanmat(26,colnum);
                    holdtmp=tmpc26;
                    for sstindx=1:size(tmpc26,2)
                        for sstrow=1:size(tmpc26,1)
                            tmpc26(sstrow,sstindx)=tmpc26(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc26=tmpc26.^2;
                    if size(tmpc26,1)==1
                        sst=sum(tmpc26(1,1:13));
                    else
                        sst=sum(tmpc26);
                    end
                    sst=sum(sst);
                    tmpc26=holdtmp;
                    dfe=(size(tmpc26,1)*size(tmpc26,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc26,2)
                        for sserow=1:size(tmpc26,1)
                            tmpc26(sserow,ssecol)=(tmpc26(sserow,ssecol)-tmpfourierwavs(26,ssecol));
                        end
                    end
                    tmpc26=tmpc26.^2;
                    sse=sum(tmpc26);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc26,1)*size(tmpc26,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(26,2)=rsquare;
                        goodnessmat(26,3)=adjrsquare;
                        goodnessmat(26,4)=sse;
                        goodnessmat(26,5)=rmse;
                        goodnessmat(26,6)=RIC;
                        goodnessmat(26,7)=sst;
                    else
                        tmpgoodnessmat(26,1)=rsquare;
                        tmpgoodnessmat(26,2)=adjrsquare;
                        tmpgoodnessmat(26,3)=sse;
                        tmpgoodnessmat(26,4)=rmse;
                        tmpgoodnessmat(26,5)=RIC;
                        tmpgoodnessmat(26,6)=sst;
                    end
                    
                    tmpc27=c27(1:end-2,2:14);
                    cmean=cmeanmat(27,colnum);
                    holdtmp=tmpc27;
                    for sstindx=1:size(tmpc27,2)
                        for sstrow=1:size(tmpc27,1)
                            tmpc27(sstrow,sstindx)=tmpc27(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc27=tmpc27.^2;
                    if size(tmpc27,1)==1
                        sst=sum(tmpc27(1,1:13));
                    else
                        sst=sum(tmpc27);
                    end
                    sst=sum(sst);
                    tmpc27=holdtmp;
                    dfe=(size(tmpc27,1)*size(tmpc27,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc27,2)
                        for sserow=1:size(tmpc27,1)
                            tmpc27(sserow,ssecol)=(tmpc27(sserow,ssecol)-tmpfourierwavs(27,ssecol));
                        end
                    end
                    tmpc27=tmpc27.^2;
                    sse=sum(tmpc27);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc27,1)*size(tmpc27,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(27,2)=rsquare;
                        goodnessmat(27,3)=adjrsquare;
                        goodnessmat(27,4)=sse;
                        goodnessmat(27,5)=rmse;
                        goodnessmat(27,6)=RIC;
                        goodnessmat(27,7)=sst;
                    else
                        tmpgoodnessmat(27,1)=rsquare;
                        tmpgoodnessmat(27,2)=adjrsquare;
                        tmpgoodnessmat(27,3)=sse;
                        tmpgoodnessmat(27,4)=rmse;
                        tmpgoodnessmat(27,5)=RIC;
                        tmpgoodnessmat(27,6)=sst;
                    end
                    
                    tmpc28=c28(1:end-2,2:14);
                    cmean=cmeanmat(28,colnum);
                    holdtmp=tmpc28;
                    for sstindx=1:size(tmpc28,2)
                        for sstrow=1:size(tmpc28,1)
                            tmpc28(sstrow,sstindx)=tmpc28(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc28=tmpc28.^2;
                    if size(tmpc28,1)==1
                        sst=sum(tmpc28(1,1:13));
                    else
                        sst=sum(tmpc28);
                    end
                    sst=sum(sst);
                    tmpc28=holdtmp;
                    dfe=(size(tmpc28,1)*size(tmpc28,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc28,2)
                        for sserow=1:size(tmpc28,1)
                            tmpc28(sserow,ssecol)=(tmpc28(sserow,ssecol)-tmpfourierwavs(28,ssecol));
                        end
                    end
                    tmpc28=tmpc28.^2;
                    sse=sum(tmpc28);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc28,1)*size(tmpc28,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(28,2)=rsquare;
                        goodnessmat(28,3)=adjrsquare;
                        goodnessmat(28,4)=sse;
                        goodnessmat(28,5)=rmse;
                        goodnessmat(28,6)=RIC;
                        goodnessmat(28,7)=sst;
                    else
                        tmpgoodnessmat(28,1)=rsquare;
                        tmpgoodnessmat(28,2)=adjrsquare;
                        tmpgoodnessmat(28,3)=sse;
                        tmpgoodnessmat(28,4)=rmse;
                        tmpgoodnessmat(28,5)=RIC;
                        tmpgoodnessmat(28,6)=sst;
                    end
                    
                    tmpc29=c29(1:end-2,2:14);
                    cmean=cmeanmat(29,colnum);
                    holdtmp=tmpc29;
                    for sstindx=1:size(tmpc29,2)
                        for sstrow=1:size(tmpc29,1)
                            tmpc29(sstrow,sstindx)=tmpc29(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc29=tmpc29.^2;
                    if size(tmpc29,1)==1
                        sst=sum(tmpc29(1,1:13));
                    else
                        sst=sum(tmpc29);
                    end
                    sst=sum(sst);
                    tmpc29=holdtmp;
                    dfe=(size(tmpc29,1)*size(tmpc29,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc29,2)
                        for sserow=1:size(tmpc29,1)
                            tmpc29(sserow,ssecol)=(tmpc29(sserow,ssecol)-tmpfourierwavs(29,ssecol));
                        end
                    end
                    tmpc29=tmpc29.^2;
                    sse=sum(tmpc29);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc29,1)*size(tmpc29,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(29,2)=rsquare;
                        goodnessmat(29,3)=adjrsquare;
                        goodnessmat(29,4)=sse;
                        goodnessmat(29,5)=rmse;
                        goodnessmat(29,6)=RIC;
                        goodnessmat(29,7)=sst;
                    else
                        tmpgoodnessmat(29,1)=rsquare;
                        tmpgoodnessmat(29,2)=adjrsquare;
                        tmpgoodnessmat(29,3)=sse;
                        tmpgoodnessmat(29,4)=rmse;
                        tmpgoodnessmat(29,5)=RIC;
                        tmpgoodnessmat(29,6)=sst;
                    end
                    
                    tmpc30=c30(1:end-2,2:14);
                    cmean=cmeanmat(30,colnum);
                    holdtmp=tmpc30;
                    for sstindx=1:size(tmpc30,2)
                        for sstrow=1:size(tmpc30,1)
                            tmpc30(sstrow,sstindx)=tmpc30(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc30=tmpc30.^2;
                    if size(tmpc30,1)==1
                        sst=sum(tmpc30(1,1:13));
                    else
                        sst=sum(tmpc30);
                    end
                    sst=sum(sst);
                    tmpc30=holdtmp;
                    dfe=(size(tmpc30,1)*size(tmpc30,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc30,2)
                        for sserow=1:size(tmpc30,1)
                            tmpc30(sserow,ssecol)=(tmpc30(sserow,ssecol)-tmpfourierwavs(30,ssecol));
                        end
                    end
                    tmpc30=tmpc30.^2;
                    sse=sum(tmpc30);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc30,1)*size(tmpc30,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(30,2)=rsquare;
                        goodnessmat(30,3)=adjrsquare;
                        goodnessmat(30,4)=sse;
                        goodnessmat(30,5)=rmse;
                        goodnessmat(30,6)=RIC;
                        goodnessmat(30,7)=sst;
                    else
                        tmpgoodnessmat(30,1)=rsquare;
                        tmpgoodnessmat(30,2)=adjrsquare;
                        tmpgoodnessmat(30,3)=sse;
                        tmpgoodnessmat(30,4)=rmse;
                        tmpgoodnessmat(30,5)=RIC;
                        tmpgoodnessmat(30,6)=sst;
                    end
                    
                    tmpc31=c31(1:end-2,2:14);
                    cmean=cmeanmat(31,colnum);
                    holdtmp=tmpc31;
                    for sstindx=1:size(tmpc31,2)
                        for sstrow=1:size(tmpc31,1)
                            tmpc31(sstrow,sstindx)=tmpc31(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc31=tmpc31.^2;
                    if size(tmpc31,1)==1
                        sst=sum(tmpc31(1,1:13));
                    else
                        sst=sum(tmpc31);
                    end
                    sst=sum(sst);
                    tmpc31=holdtmp;
                    dfe=(size(tmpc31,1)*size(tmpc31,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc31,2)
                        for sserow=1:size(tmpc31,1)
                            tmpc31(sserow,ssecol)=(tmpc31(sserow,ssecol)-tmpfourierwavs(31,ssecol));
                        end
                    end
                    tmpc31=tmpc31.^2;
                    sse=sum(tmpc31);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc31,1)*size(tmpc31,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(31,2)=rsquare;
                        goodnessmat(31,3)=adjrsquare;
                        goodnessmat(31,4)=sse;
                        goodnessmat(31,5)=rmse;
                        goodnessmat(31,6)=RIC;
                        goodnessmat(31,7)=sst;
                    else
                        tmpgoodnessmat(31,1)=rsquare;
                        tmpgoodnessmat(31,2)=adjrsquare;
                        tmpgoodnessmat(31,3)=sse;
                        tmpgoodnessmat(31,4)=rmse;
                        tmpgoodnessmat(31,5)=RIC;
                        tmpgoodnessmat(31,6)=sst;
                    end
                    
                    tmpc32=c32(1:end-2,2:14);
                    cmean=cmeanmat(32,colnum);
                    holdtmp=tmpc32;
                    for sstindx=1:size(tmpc32,2)
                        for sstrow=1:size(tmpc32,1)
                            tmpc32(sstrow,sstindx)=tmpc32(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc32=tmpc32.^2;
                    if size(tmpc32,1)==1
                        sst=sum(tmpc32(1,1:13));
                    else
                        sst=sum(tmpc32);
                    end
                    sst=sum(sst);
                    tmpc32=holdtmp;
                    dfe=(size(tmpc32,1)*size(tmpc32,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc32,2)
                        for sserow=1:size(tmpc32,1)
                            tmpc32(sserow,ssecol)=(tmpc32(sserow,ssecol)-tmpfourierwavs(32,ssecol));
                        end
                    end
                    tmpc32=tmpc32.^2;
                    sse=sum(tmpc32);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc32,1)*size(tmpc32,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(32,2)=rsquare;
                        goodnessmat(32,3)=adjrsquare;
                        goodnessmat(32,4)=sse;
                        goodnessmat(32,5)=rmse;
                        goodnessmat(32,6)=RIC;
                        goodnessmat(32,7)=sst;
                    else
                        tmpgoodnessmat(32,1)=rsquare;
                        tmpgoodnessmat(32,2)=adjrsquare;
                        tmpgoodnessmat(32,3)=sse;
                        tmpgoodnessmat(32,4)=rmse;
                        tmpgoodnessmat(32,5)=RIC;
                        tmpgoodnessmat(32,6)=sst;
                    end
                    
                    tmpc33=c33(1:end-2,2:14);
                    cmean=cmeanmat(33,colnum);
                    holdtmp=tmpc33;
                    for sstindx=1:size(tmpc33,2)
                        for sstrow=1:size(tmpc33,1)
                            tmpc33(sstrow,sstindx)=tmpc33(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc33=tmpc33.^2;
                    if size(tmpc33,1)==1
                        sst=sum(tmpc33(1,1:13));
                    else
                        sst=sum(tmpc33);
                    end
                    sst=sum(sst);
                    tmpc33=holdtmp;
                    dfe=(size(tmpc33,1)*size(tmpc33,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc33,2)
                        for sserow=1:size(tmpc33,1)
                            tmpc33(sserow,ssecol)=(tmpc33(sserow,ssecol)-tmpfourierwavs(33,ssecol));
                        end
                    end
                    tmpc33=tmpc33.^2;
                    sse=sum(tmpc33);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc33,1)*size(tmpc33,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(33,2)=rsquare;
                        goodnessmat(33,3)=adjrsquare;
                        goodnessmat(33,4)=sse;
                        goodnessmat(33,5)=rmse;
                        goodnessmat(33,6)=RIC;
                        goodnessmat(33,7)=sst;
                    else
                        tmpgoodnessmat(33,1)=rsquare;
                        tmpgoodnessmat(33,2)=adjrsquare;
                        tmpgoodnessmat(33,3)=sse;
                        tmpgoodnessmat(33,4)=rmse;
                        tmpgoodnessmat(33,5)=RIC;
                        tmpgoodnessmat(33,6)=sst;
                    end
                    
                    tmpc34=c34(1:end-2,2:14);
                    cmean=cmeanmat(34,colnum);
                    holdtmp=tmpc34;
                    for sstindx=1:size(tmpc34,2)
                        for sstrow=1:size(tmpc34,1)
                            tmpc34(sstrow,sstindx)=tmpc34(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc34=tmpc34.^2;
                    if size(tmpc34,1)==1
                        sst=sum(tmpc34(1,1:13));
                    else
                        sst=sum(tmpc34);
                    end
                    sst=sum(sst);
                    tmpc34=holdtmp;
                    dfe=(size(tmpc34,1)*size(tmpc34,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc34,2)
                        for sserow=1:size(tmpc34,1)
                            tmpc34(sserow,ssecol)=(tmpc34(sserow,ssecol)-tmpfourierwavs(34,ssecol));
                        end
                    end
                    tmpc34=tmpc34.^2;
                    sse=sum(tmpc34);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc34,1)*size(tmpc34,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(34,2)=rsquare;
                        goodnessmat(34,3)=adjrsquare;
                        goodnessmat(34,4)=sse;
                        goodnessmat(34,5)=rmse;
                        goodnessmat(34,6)=RIC;
                        goodnessmat(34,7)=sst;
                    else
                        tmpgoodnessmat(34,1)=rsquare;
                        tmpgoodnessmat(34,2)=adjrsquare;
                        tmpgoodnessmat(34,3)=sse;
                        tmpgoodnessmat(34,4)=rmse;
                        tmpgoodnessmat(34,5)=RIC;
                        tmpgoodnessmat(34,6)=sst;
                    end
                    
                    tmpc35=c35(1:end-2,2:14);
                    cmean=cmeanmat(35,colnum);
                    holdtmp=tmpc35;
                    for sstindx=1:size(tmpc35,2)
                        for sstrow=1:size(tmpc35,1)
                            tmpc35(sstrow,sstindx)=tmpc35(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc35=tmpc35.^2;
                    if size(tmpc35,1)==1
                        sst=sum(tmpc35(1,1:13));
                    else
                        sst=sum(tmpc35);
                    end
                    sst=sum(sst);
                    tmpc35=holdtmp;
                    dfe=(size(tmpc35,1)*size(tmpc35,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc35,2)
                        for sserow=1:size(tmpc35,1)
                            tmpc35(sserow,ssecol)=(tmpc35(sserow,ssecol)-tmpfourierwavs(35,ssecol));
                        end
                    end
                    tmpc35=tmpc35.^2;
                    sse=sum(tmpc35);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc35,1)*size(tmpc35,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(35,2)=rsquare;
                        goodnessmat(35,3)=adjrsquare;
                        goodnessmat(35,4)=sse;
                        goodnessmat(35,5)=rmse;
                        goodnessmat(35,6)=RIC;
                        goodnessmat(35,7)=sst;
                    else
                        tmpgoodnessmat(35,1)=rsquare;
                        tmpgoodnessmat(35,2)=adjrsquare;
                        tmpgoodnessmat(35,3)=sse;
                        tmpgoodnessmat(35,4)=rmse;
                        tmpgoodnessmat(35,5)=RIC;
                        tmpgoodnessmat(35,6)=sst;
                    end
                    
                    tmpc36=c36(1:end-2,2:14);
                    cmean=cmeanmat(36,colnum);
                    holdtmp=tmpc36;
                    for sstindx=1:size(tmpc36,2)
                        for sstrow=1:size(tmpc36,1)
                            tmpc36(sstrow,sstindx)=tmpc36(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc36=tmpc36.^2;
                    if size(tmpc36,1)==1
                        sst=sum(tmpc36(1,1:13));
                    else
                        sst=sum(tmpc36);
                    end
                    sst=sum(sst);
                    tmpc36=holdtmp;
                    dfe=(size(tmpc36,1)*size(tmpc36,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc36,2)
                        for sserow=1:size(tmpc36,1)
                            tmpc36(sserow,ssecol)=(tmpc36(sserow,ssecol)-tmpfourierwavs(36,ssecol));
                        end
                    end
                    tmpc36=tmpc36.^2;
                    sse=sum(tmpc36);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc36,1)*size(tmpc36,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(36,2)=rsquare;
                        goodnessmat(36,3)=adjrsquare;
                        goodnessmat(36,4)=sse;
                        goodnessmat(36,5)=rmse;
                        goodnessmat(36,6)=RIC;
                        goodnessmat(36,7)=sst;
                    else
                        tmpgoodnessmat(36,1)=rsquare;
                        tmpgoodnessmat(36,2)=adjrsquare;
                        tmpgoodnessmat(36,3)=sse;
                        tmpgoodnessmat(36,4)=rmse;
                        tmpgoodnessmat(36,5)=RIC;
                        tmpgoodnessmat(36,6)=sst;
                    end
                    
                    tmpc37=c37(1:end-2,2:14);
                    cmean=cmeanmat(37,colnum);
                    holdtmp=tmpc37;
                    for sstindx=1:size(tmpc37,2)
                        for sstrow=1:size(tmpc37,1)
                            tmpc37(sstrow,sstindx)=tmpc37(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc37=tmpc37.^2;
                    if size(tmpc37,1)==1
                        sst=sum(tmpc37(1,1:13));
                    else
                        sst=sum(tmpc37);
                    end
                    sst=sum(sst);
                    tmpc37=holdtmp;
                    dfe=(size(tmpc37,1)*size(tmpc37,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc37,2)
                        for sserow=1:size(tmpc37,1)
                            tmpc37(sserow,ssecol)=(tmpc37(sserow,ssecol)-tmpfourierwavs(37,ssecol));
                        end
                    end
                    tmpc37=tmpc37.^2;
                    sse=sum(tmpc37);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc37,1)*size(tmpc37,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(37,2)=rsquare;
                        goodnessmat(37,3)=adjrsquare;
                        goodnessmat(37,4)=sse;
                        goodnessmat(37,5)=rmse;
                        goodnessmat(37,6)=RIC;
                        goodnessmat(37,7)=sst;
                    else
                        tmpgoodnessmat(37,1)=rsquare;
                        tmpgoodnessmat(37,2)=adjrsquare;
                        tmpgoodnessmat(37,3)=sse;
                        tmpgoodnessmat(37,4)=rmse;
                        tmpgoodnessmat(37,5)=RIC;
                        tmpgoodnessmat(37,6)=sst;
                    end
                    
                    tmpc38=c38(1:end-2,2:14);
                    cmean=cmeanmat(38,colnum);
                    holdtmp=tmpc38;
                    for sstindx=1:size(tmpc38,2)
                        for sstrow=1:size(tmpc38,1)
                            tmpc38(sstrow,sstindx)=tmpc38(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc38=tmpc38.^2;
                    if size(tmpc38,1)==1
                        sst=sum(tmpc38(1,1:13));
                    else
                        sst=sum(tmpc38);
                    end
                    sst=sum(sst);
                    tmpc38=holdtmp;
                    dfe=(size(tmpc38,1)*size(tmpc38,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc38,2)
                        for sserow=1:size(tmpc38,1)
                            tmpc38(sserow,ssecol)=(tmpc38(sserow,ssecol)-tmpfourierwavs(38,ssecol));
                        end
                    end
                    tmpc38=tmpc38.^2;
                    sse=sum(tmpc38);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc38,1)*size(tmpc38,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(38,2)=rsquare;
                        goodnessmat(38,3)=adjrsquare;
                        goodnessmat(38,4)=sse;
                        goodnessmat(38,5)=rmse;
                        goodnessmat(38,6)=RIC;
                        goodnessmat(38,7)=sst;
                    else
                        tmpgoodnessmat(38,1)=rsquare;
                        tmpgoodnessmat(38,2)=adjrsquare;
                        tmpgoodnessmat(38,3)=sse;
                        tmpgoodnessmat(38,4)=rmse;
                        tmpgoodnessmat(38,5)=RIC;
                        tmpgoodnessmat(38,6)=sst;
                    end
                    
                    tmpc39=c39(1:end-2,2:14);
                    cmean=cmeanmat(39,colnum);
                    holdtmp=tmpc39;
                    for sstindx=1:size(tmpc39,2)
                        for sstrow=1:size(tmpc39,1)
                            tmpc39(sstrow,sstindx)=tmpc39(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc39=tmpc39.^2;
                    if size(tmpc39,1)==1
                        sst=sum(tmpc39(1,1:13));
                    else
                        sst=sum(tmpc39);
                    end
                    sst=sum(sst);
                    tmpc39=holdtmp;
                    dfe=(size(tmpc39,1)*size(tmpc39,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc39,2)
                        for sserow=1:size(tmpc39,1)
                            tmpc39(sserow,ssecol)=(tmpc39(sserow,ssecol)-tmpfourierwavs(39,ssecol));
                        end
                    end
                    tmpc39=tmpc39.^2;
                    sse=sum(tmpc39);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc39,1)*size(tmpc39,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(39,2)=rsquare;
                        goodnessmat(39,3)=adjrsquare;
                        goodnessmat(39,4)=sse;
                        goodnessmat(39,5)=rmse;
                        goodnessmat(39,6)=RIC;
                        goodnessmat(39,7)=sst;
                    else
                        tmpgoodnessmat(39,1)=rsquare;
                        tmpgoodnessmat(39,2)=adjrsquare;
                        tmpgoodnessmat(39,3)=sse;
                        tmpgoodnessmat(39,4)=rmse;
                        tmpgoodnessmat(39,5)=RIC;
                        tmpgoodnessmat(39,6)=sst;
                    end
                    
                    tmpc40=c40(1:end-2,2:14);
                    cmean=cmeanmat(40,colnum);
                    holdtmp=tmpc40;
                    for sstindx=1:size(tmpc40,2)
                        for sstrow=1:size(tmpc40,1)
                            tmpc40(sstrow,sstindx)=tmpc40(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc40=tmpc40.^2;
                    if size(tmpc40,1)==1
                        sst=sum(tmpc40(1,1:13));
                    else
                        sst=sum(tmpc40);
                    end
                    sst=sum(sst);
                    tmpc40=holdtmp;
                    dfe=(size(tmpc40,1)*size(tmpc40,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc40,2)
                        for sserow=1:size(tmpc40,1)
                            tmpc40(sserow,ssecol)=(tmpc40(sserow,ssecol)-tmpfourierwavs(40,ssecol));
                        end
                    end
                    tmpc40=tmpc40.^2;
                    sse=sum(tmpc40);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc40,1)*size(tmpc40,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(40,2)=rsquare;
                        goodnessmat(40,3)=adjrsquare;
                        goodnessmat(40,4)=sse;
                        goodnessmat(40,5)=rmse;
                        goodnessmat(40,6)=RIC;
                        goodnessmat(40,7)=sst;
                    else
                        tmpgoodnessmat(40,1)=rsquare;
                        tmpgoodnessmat(40,2)=adjrsquare;
                        tmpgoodnessmat(40,3)=sse;
                        tmpgoodnessmat(40,4)=rmse;
                        tmpgoodnessmat(40,5)=RIC;
                        tmpgoodnessmat(40,6)=sst;
                    end
                    
                    tmpc41=c41(1:end-2,2:14);
                    cmean=cmeanmat(41,colnum);
                    holdtmp=tmpc41;
                    for sstindx=1:size(tmpc41,2)
                        for sstrow=1:size(tmpc41,1)
                            tmpc41(sstrow,sstindx)=tmpc41(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc41=tmpc41.^2;
                    if size(tmpc41,1)==1
                        sst=sum(tmpc41(1,1:13));
                    else
                        sst=sum(tmpc41);
                    end
                    sst=sum(sst);
                    tmpc41=holdtmp;
                    dfe=(size(tmpc41,1)*size(tmpc41,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc41,2)
                        for sserow=1:size(tmpc41,1)
                            tmpc41(sserow,ssecol)=(tmpc41(sserow,ssecol)-tmpfourierwavs(41,ssecol));
                        end
                    end
                    tmpc41=tmpc41.^2;
                    sse=sum(tmpc41);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc41,1)*size(tmpc41,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(41,2)=rsquare;
                        goodnessmat(41,3)=adjrsquare;
                        goodnessmat(41,4)=sse;
                        goodnessmat(41,5)=rmse;
                        goodnessmat(41,6)=RIC;
                        goodnessmat(41,7)=sst;
                    else
                        tmpgoodnessmat(41,1)=rsquare;
                        tmpgoodnessmat(41,2)=adjrsquare;
                        tmpgoodnessmat(41,3)=sse;
                        tmpgoodnessmat(41,4)=rmse;
                        tmpgoodnessmat(41,5)=RIC;
                        tmpgoodnessmat(41,6)=sst;
                    end
                    
                    tmpc42=c42(1:end-2,2:14);
                    cmean=cmeanmat(42,colnum);
                    holdtmp=tmpc42;
                    for sstindx=1:size(tmpc42,2)
                        for sstrow=1:size(tmpc42,1)
                            tmpc42(sstrow,sstindx)=tmpc42(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc42=tmpc42.^2;
                    if size(tmpc42,1)==1
                        sst=sum(tmpc42(1,1:13));
                    else
                        sst=sum(tmpc42);
                    end
                    sst=sum(sst);
                    tmpc42=holdtmp;
                    dfe=(size(tmpc42,1)*size(tmpc42,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc42,2)
                        for sserow=1:size(tmpc42,1)
                            tmpc42(sserow,ssecol)=(tmpc42(sserow,ssecol)-tmpfourierwavs(42,ssecol));
                        end
                    end
                    tmpc42=tmpc42.^2;
                    sse=sum(tmpc42);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc42,1)*size(tmpc42,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(42,2)=rsquare;
                        goodnessmat(42,3)=adjrsquare;
                        goodnessmat(42,4)=sse;
                        goodnessmat(42,5)=rmse;
                        goodnessmat(42,6)=RIC;
                        goodnessmat(42,7)=sst;
                    else
                        tmpgoodnessmat(42,1)=rsquare;
                        tmpgoodnessmat(42,2)=adjrsquare;
                        tmpgoodnessmat(42,3)=sse;
                        tmpgoodnessmat(42,4)=rmse;
                        tmpgoodnessmat(42,5)=RIC;
                        tmpgoodnessmat(42,6)=sst;
                    end
                    
                    tmpc43=c43(1:end-2,2:14);
                    cmean=cmeanmat(43,colnum);
                    holdtmp=tmpc43;
                    for sstindx=1:size(tmpc43,2)
                        for sstrow=1:size(tmpc43,1)
                            tmpc43(sstrow,sstindx)=tmpc43(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc43=tmpc43.^2;
                    if size(tmpc43,1)==1
                        sst=sum(tmpc43(1,1:13));
                    else
                        sst=sum(tmpc43);
                    end
                    sst=sum(sst);
                    tmpc43=holdtmp;
                    dfe=(size(tmpc43,1)*size(tmpc43,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc43,2)
                        for sserow=1:size(tmpc43,1)
                            tmpc43(sserow,ssecol)=(tmpc43(sserow,ssecol)-tmpfourierwavs(43,ssecol));
                        end
                    end
                    tmpc43=tmpc43.^2;
                    sse=sum(tmpc43);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc43,1)*size(tmpc43,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(43,2)=rsquare;
                        goodnessmat(43,3)=adjrsquare;
                        goodnessmat(43,4)=sse;
                        goodnessmat(43,5)=rmse;
                        goodnessmat(43,6)=RIC;
                        goodnessmat(43,7)=sst;
                    else
                        tmpgoodnessmat(43,1)=rsquare;
                        tmpgoodnessmat(43,2)=adjrsquare;
                        tmpgoodnessmat(43,3)=sse;
                        tmpgoodnessmat(43,4)=rmse;
                        tmpgoodnessmat(43,5)=RIC;
                        tmpgoodnessmat(43,6)=sst;
                    end
                    
                    tmpc44=c44(1:end-2,2:14);
                    cmean=cmeanmat(44,colnum);
                    holdtmp=tmpc44;
                    for sstindx=1:size(tmpc44,2)
                        for sstrow=1:size(tmpc44,1)
                            tmpc44(sstrow,sstindx)=tmpc44(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc44=tmpc44.^2;
                    if size(tmpc44,1)==1
                        sst=sum(tmpc44(1,1:13));
                    else
                        sst=sum(tmpc44);
                    end
                    sst=sum(sst);
                    tmpc44=holdtmp;
                    dfe=(size(tmpc44,1)*size(tmpc44,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc44,2)
                        for sserow=1:size(tmpc44,1)
                            tmpc44(sserow,ssecol)=(tmpc44(sserow,ssecol)-tmpfourierwavs(44,ssecol));
                        end
                    end
                    tmpc44=tmpc44.^2;
                    sse=sum(tmpc44);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc44,1)*size(tmpc44,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(44,2)=rsquare;
                        goodnessmat(44,3)=adjrsquare;
                        goodnessmat(44,4)=sse;
                        goodnessmat(44,5)=rmse;
                        goodnessmat(44,6)=RIC;
                        goodnessmat(44,7)=sst;
                    else
                        tmpgoodnessmat(44,1)=rsquare;
                        tmpgoodnessmat(44,2)=adjrsquare;
                        tmpgoodnessmat(44,3)=sse;
                        tmpgoodnessmat(44,4)=rmse;
                        tmpgoodnessmat(44,5)=RIC;
                        tmpgoodnessmat(44,6)=sst;
                    end
                    
                    tmpc45=c45(1:end-2,2:14);
                    cmean=cmeanmat(45,colnum);
                    holdtmp=tmpc45;
                    for sstindx=1:size(tmpc45,2)
                        for sstrow=1:size(tmpc45,1)
                            tmpc45(sstrow,sstindx)=tmpc45(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc45=tmpc45.^2;
                    if size(tmpc45,1)==1
                        sst=sum(tmpc45(1,1:13));
                    else
                        sst=sum(tmpc45);
                    end
                    sst=sum(sst);
                    tmpc45=holdtmp;
                    dfe=(size(tmpc45,1)*size(tmpc45,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc45,2)
                        for sserow=1:size(tmpc45,1)
                            tmpc45(sserow,ssecol)=(tmpc45(sserow,ssecol)-tmpfourierwavs(45,ssecol));
                        end
                    end
                    tmpc45=tmpc45.^2;
                    sse=sum(tmpc45);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc45,1)*size(tmpc45,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(45,2)=rsquare;
                        goodnessmat(45,3)=adjrsquare;
                        goodnessmat(45,4)=sse;
                        goodnessmat(45,5)=rmse;
                        goodnessmat(45,6)=RIC;
                        goodnessmat(45,7)=sst;
                    else
                        tmpgoodnessmat(45,1)=rsquare;
                        tmpgoodnessmat(45,2)=adjrsquare;
                        tmpgoodnessmat(45,3)=sse;
                        tmpgoodnessmat(45,4)=rmse;
                        tmpgoodnessmat(45,5)=RIC;
                        tmpgoodnessmat(45,6)=sst;
                    end
                    
                    tmpc46=c46(1:end-2,2:14);
                    cmean=cmeanmat(46,colnum);
                    holdtmp=tmpc46;
                    for sstindx=1:size(tmpc46,2)
                        for sstrow=1:size(tmpc46,1)
                            tmpc46(sstrow,sstindx)=tmpc46(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc46=tmpc46.^2;
                    if size(tmpc46,1)==1
                        sst=sum(tmpc46(1,1:13));
                    else
                        sst=sum(tmpc46);
                    end
                    sst=sum(sst);
                    tmpc46=holdtmp;
                    dfe=(size(tmpc46,1)*size(tmpc46,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc46,2)
                        for sserow=1:size(tmpc46,1)
                            tmpc46(sserow,ssecol)=(tmpc46(sserow,ssecol)-tmpfourierwavs(46,ssecol));
                        end
                    end
                    tmpc46=tmpc46.^2;
                    sse=sum(tmpc46);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc46,1)*size(tmpc46,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(46,2)=rsquare;
                        goodnessmat(46,3)=adjrsquare;
                        goodnessmat(46,4)=sse;
                        goodnessmat(46,5)=rmse;
                        goodnessmat(46,6)=RIC;
                        goodnessmat(46,7)=sst;
                    else
                        tmpgoodnessmat(46,1)=rsquare;
                        tmpgoodnessmat(46,2)=adjrsquare;
                        tmpgoodnessmat(46,3)=sse;
                        tmpgoodnessmat(46,4)=rmse;
                        tmpgoodnessmat(46,5)=RIC;
                        tmpgoodnessmat(46,6)=sst;
                    end
                    
                    tmpc47=c47(1:end-2,2:14);
                    cmean=cmeanmat(47,colnum);
                    holdtmp=tmpc47;
                    for sstindx=1:size(tmpc47,2)
                        for sstrow=1:size(tmpc47,1)
                            tmpc47(sstrow,sstindx)=tmpc47(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc47=tmpc47.^2;
                    if size(tmpc47,1)==1
                        sst=sum(tmpc47(1,1:13));
                    else
                        sst=sum(tmpc47);
                    end
                    sst=sum(sst);
                    tmpc47=holdtmp;
                    dfe=(size(tmpc47,1)*size(tmpc47,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc47,2)
                        for sserow=1:size(tmpc47,1)
                            tmpc47(sserow,ssecol)=(tmpc47(sserow,ssecol)-tmpfourierwavs(47,ssecol));
                        end
                    end
                    tmpc47=tmpc47.^2;
                    sse=sum(tmpc47);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc47,1)*size(tmpc47,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(47,2)=rsquare;
                        goodnessmat(47,3)=adjrsquare;
                        goodnessmat(47,4)=sse;
                        goodnessmat(47,5)=rmse;
                        goodnessmat(47,6)=RIC;
                        goodnessmat(47,7)=sst;
                    else
                        tmpgoodnessmat(47,1)=rsquare;
                        tmpgoodnessmat(47,2)=adjrsquare;
                        tmpgoodnessmat(47,3)=sse;
                        tmpgoodnessmat(47,4)=rmse;
                        tmpgoodnessmat(47,5)=RIC;
                        tmpgoodnessmat(47,6)=sst;
                    end
                    
                    tmpc48=c48(1:end-2,2:14);
                    cmean=cmeanmat(48,colnum);
                    holdtmp=tmpc48;
                    for sstindx=1:size(tmpc48,2)
                        for sstrow=1:size(tmpc48,1)
                            tmpc48(sstrow,sstindx)=tmpc48(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc48=tmpc48.^2;
                    if size(tmpc48,1)==1
                        sst=sum(tmpc48(1,1:13));
                    else
                        sst=sum(tmpc48);
                    end
                    sst=sum(sst);
                    tmpc48=holdtmp;
                    dfe=(size(tmpc48,1)*size(tmpc48,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc48,2)
                        for sserow=1:size(tmpc48,1)
                            tmpc48(sserow,ssecol)=(tmpc48(sserow,ssecol)-tmpfourierwavs(48,ssecol));
                        end
                    end
                    tmpc48=tmpc48.^2;
                    sse=sum(tmpc48);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc48,1)*size(tmpc48,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(48,2)=rsquare;
                        goodnessmat(48,3)=adjrsquare;
                        goodnessmat(48,4)=sse;
                        goodnessmat(48,5)=rmse;
                        goodnessmat(48,6)=RIC;
                        goodnessmat(48,7)=sst;
                    else
                        tmpgoodnessmat(48,1)=rsquare;
                        tmpgoodnessmat(48,2)=adjrsquare;
                        tmpgoodnessmat(48,3)=sse;
                        tmpgoodnessmat(48,4)=rmse;
                        tmpgoodnessmat(48,5)=RIC;
                        tmpgoodnessmat(48,6)=sst;
                    end
                    
                    tmpc49=c49(1:end-2,2:14);
                    cmean=cmeanmat(49,colnum);
                    holdtmp=tmpc49;
                    for sstindx=1:size(tmpc49,2)
                        for sstrow=1:size(tmpc49,1)
                            tmpc49(sstrow,sstindx)=tmpc49(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc49=tmpc49.^2;
                    if size(tmpc49,1)==1
                        sst=sum(tmpc49(1,1:13));
                    else
                        sst=sum(tmpc49);
                    end
                    sst=sum(sst);
                    tmpc49=holdtmp;
                    dfe=(size(tmpc49,1)*size(tmpc49,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc49,2)
                        for sserow=1:size(tmpc49,1)
                            tmpc49(sserow,ssecol)=(tmpc49(sserow,ssecol)-tmpfourierwavs(49,ssecol));
                        end
                    end
                    tmpc49=tmpc49.^2;
                    sse=sum(tmpc49);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    
                    adjrsquare=1-(1-rsquare)*((size(tmpc49,1)*size(tmpc49,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(49,2)=rsquare;
                        goodnessmat(49,3)=adjrsquare;
                        goodnessmat(49,4)=sse;
                        goodnessmat(49,5)=rmse;
                        goodnessmat(49,6)=RIC;
                        goodnessmat(49,7)=sst;
                    else
                        tmpgoodnessmat(49,1)=rsquare;
                        tmpgoodnessmat(49,2)=adjrsquare;
                        tmpgoodnessmat(49,3)=sse;
                        tmpgoodnessmat(49,4)=rmse;
                        tmpgoodnessmat(49,5)=RIC;
                        tmpgoodnessmat(49,6)=sst;
                    end
                    
                    tmpc50=c50(1:end-2,2:14);
                    cmean=cmeanmat(50,colnum);
                    holdtmp=tmpc50;
                    for sstindx=1:size(tmpc50,2)
                        for sstrow=1:size(tmpc50,1)
                            tmpc50(sstrow,sstindx)=tmpc50(sstrow,sstindx)-cmean;
                        end
                    end
                    tmpc50=tmpc50.^2;
                    if size(tmpc50,1)==1
                        sst=sum(tmpc50(1,1:13));
                    else
                        sst=sum(tmpc50);
                    end
                    sst=sum(sst);
                    tmpc50=holdtmp;
                    dfe=(size(tmpc50,1)*size(tmpc50,2))-12;
                    RIC=0;
                    for ssecol=1:size(tmpc50,2)
                        for sserow=1:size(tmpc50,1)
                            tmpc50(sserow,ssecol)=(tmpc50(sserow,ssecol)-tmpfourierwavs(50,ssecol));
                        end
                    end
                    tmpc50=tmpc50.^2;
                    sse=sum(tmpc50);
                    sse=sum(sse)/dfe;
                    rsquare=1-(sse/sst);
                    adjrsquare=1-(1-rsquare)*((size(tmpc50,1)*size(tmpc50,2))-1)/dfe;
                    rmse=sqrt(sse);
                    
                    if frchckindx==1
                        goodnessmat(50,2)=rsquare;
                        goodnessmat(50,3)=adjrsquare;
                        goodnessmat(50,4)=sse;
                        goodnessmat(50,5)=rmse;
                        goodnessmat(50,6)=RIC;
                        goodnessmat(50,7)=sst;
                    else
                        tmpgoodnessmat(50,1)=rsquare;
                        tmpgoodnessmat(50,2)=adjrsquare;
                        tmpgoodnessmat(50,3)=sse;
                        tmpgoodnessmat(50,4)=rmse;
                        tmpgoodnessmat(50,5)=RIC;
                        tmpgoodnessmat(50,6)=sst;
                    end

                    
                    if frchckindx==1
                        goodnessmat1(:,1:7)=goodnessmat;
                    else
                        goodnessmat1(:,gdnsmat1cnt:gdnsmat1cnt+5)=tmpgoodnessmat;
                        gdnsmat1cnt=gdnsmat1cnt+6;
                    end
                end
            end
            
            goodnessmat1=goodnessmat1(:,:);
            gdnsmatfincnt=361;
                
            if intindx1==11
                goodnessmatfin(:,1:361)=goodnessmat1;
            else
                goodnessmatfin(:,gdnsmatfincnt:gdnsmatfincnt+360)=goodnessmat1;
                gdnsmatfincnt=gdnsmatfincnt+361;
            end
        end
    end
    subtractmat=zeros(50,16245);
    
    for subtractindx=1:361:16245
        subtractmat(:,subtractindx:subtractindx+360)=submatGOF;
    end
    goodnessmatfin=abs(goodnessmatfin-subtractmat);
    %TEST
    for clnindx=1:361:16245
        goodnessmatfin(:,clnindx)=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50];
    end
    %ENDTEST
    %mex_WriteMatrix('p1.csvx.gofcheck.csv', goodnessmatfin,'%2.2f',',','w+');
    switch switchindex
        case 1
            dlmwrite('p51.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 2
            dlmwrite('p52.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 3
            dlmwrite('p53.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 4
            dlmwrite('p54.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 5
            dlmwrite('p55.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 6
            dlmwrite('p56.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 7
            dlmwrite('p57.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 8
            dlmwrite('p58.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 9
            dlmwrite('p59.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 10
            dlmwrite('p60.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 11
            dlmwrite('p61.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 12
            dlmwrite('p62.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 13
            dlmwrite('p63.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 14
            dlmwrite('p64.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 15
            dlmwrite('p65.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 16
            dlmwrite('p66.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 17
            dlmwrite('p67.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 18
            dlmwrite('p68.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 19
            dlmwrite('p69.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 20
            dlmwrite('p70.csvx.gofcheck.csv', goodnessmatfin);
           goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 21
            dlmwrite('p71.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 22
            dlmwrite('p72.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 23
            dlmwrite('p73.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 24
            dlmwrite('p74.csvx.gofcheck.csv', goodnessmatfin);
            goodnessmatfin(:,:)=0;
            subtractmat(:,:)=0;
        case 25
            dlmwrite('p75.csvx.gofcheck.csv', goodnessmatfin);
 
    end
end
