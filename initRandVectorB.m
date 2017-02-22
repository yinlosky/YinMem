function initRandVectorB(NumOfNodes,NumOfProcessors, NumOfMachines)
%%
%%
%%
myDB;
machines = getMachines(NumOfMachines);
fname = ('benchmark/version5_stat.txt');
fstat = fopen(fname,'a+');

%B=DB(['B' num2str(NumOfNodes)]);
%delete(B);
%B=DB(['B' num2str(NumOfNodes)]);

%disp(['Initializing the random vector b in table B' num2str(NumOfNodes)]);
%	this = tic;eval(pRUN('initB',NumOfProcessors,machines));that = toc(this);
%	disp(['InitB takes: ' num2str(that)]);
%	fwrite(fstat,['InitB takes: ' num2str(that) sprintf('\n')]);
    
    %% initialize the lz_q1 for the for loop
    eval(pRUN('parallel_lz_norm_B_p1',NumOfProcessors,machines));
    parallel_lz_norm_B_p2(NumOfProcessors);
    eval(pRUN('parallel_scalarmult_B',NumOfProcessors,machines));
