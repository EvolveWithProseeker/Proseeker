%Function:          Extract physical motifs.
%inputs:            csv consisting of only numeric values and BLANKS from
%                   p*.csvx as output by profilecreator.m
%outputs:           Physical profile of user selected sites to serve as
%                   search motifs.
%other              Numbering is in serial but seperately tracked for
%                   structural and functional profiles. i.e. You can have a
%                   p1 struct and a p1 func but not 2 x p1 func.
disp('I number everything according to file sequence numbers IN SERIAL. I have two files where I keep the last value used so please do not delete them or I will get very cross!')
disp('I also do this in loops with minimal user input so I am going to be asking you the same questions repeated for each profile in their natural order.')
sophies=input('Is this going to be a structural profile (enter 1) or a functional profile (enter 2)     ');
files = dir( fullfile('*bp.csvx') );   
files = {files.name}'; 
files = sort_nat(files);

if sophies==2
    for sdindx=1:size(files,1)
    aaug=importdata(files{sdindx,1});
    disp(files{sdindx,1})
    fseqnum=importdata('seqnumfunc.dat');
    bres1=input('Binding residue 1:     ');
    bres2=input('Binding residue 2:     ');
    textfilename1 = ['p' num2str(fseqnum) '.bres1.x.csv'];
    textfilename2 = ['p' num2str(fseqnum) '.bres2.x.csv'];
    a1=aaug(:,bres1-6:bres1+6);
    b1=aaug(:,bres2-6:bres2+6);
    csvwrite(textfilename1,a1);
    csvwrite(textfilename2,b1);
    fseqnum=fseqnum+1;
    delete seqnumfunc.dat;
    dlmwrite('seqnumfunc.dat',fseqnum);
    end
elseif sophies==1
    disp('Files from which structural profiles are taken should NOT be circularized. This prediction process is different from the functional prediction process.')
    disp('Please DO NOT take multiple profiles of the same site labelling them as different secondary structures. This confuses the system and renders the scoring ineffective')
    disp('Profiles should encompass only a SINGLE secondary structure. Please adhere to this or the wave functions generated will be worse than useless.')
    disp('In addition you need to be absolutely sure of the secondary structure of any area to be profiled within the cellular environment. Simulation data is NOT good enough as it needs a crystal structure as ambiguity results in inaccurate profiles.')
    for sdindx=1:size(files,1)
    aaug=importdata(files{sdindx,1});
    disp(files{sdindx,1})
    fseqnum=importdata('seqnumstruct.dat');
    disp('I now need to know the starting residue of the secondary structure that you want a profile of and the residue on which the structure ends.')
    stbresstart=input('Starting residue:     ');
    stbresend=input('Ending residue:     ');
    disp('I now need to know what kind of secondary structure this profile refers to:')
    disp('1 = Alpha helix')
    disp('2 = Beta sheet')
    disp('3 = Bend')
    disp('4 = 310 helix')
    disp('5 = Pi helix')
    disp('6 = Hydrogen bonded turn')
    disp('7 = Other (meaning not coil but some unlisted secondary structure')
    struct=input('What kind of structure does this refer to (only ONE selection may be made)?     ');
    stwin=aaug(:,stbresstart:stbresend);
    denom=zeros(1,size(stwin,2));
    for denomindx=1:size(denom,2)
        denom(1,demonindx)=struct;
    end
    stwin=[denom stwin];
    textfilename1 = ['p' num2str(fseqnum) '.bres1.st.csv'];
    csvwrite(textfilename1,stwin);
    fseqnum=fseqnum+1;
    dlmwrite('seqnumstruct.dat',fseqnum);
    end
end
