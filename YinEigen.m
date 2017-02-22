function YinEigen(NumOfMachines,  NumOfProcessors, NumOfNodes, initMat, EdgesPerVertex, max_iteration, eig_k, KeepB, Run_schedule, StoreTFS, StoreLHD, TFS)
%%
%% Usage: This function is used for calculating the eigen values and eigen vectors for input symmetric matrix 
%% with the size of NumOfNodes using Alluxio in-memory file system to store the input matrix 
%%     Input: 1.NumOfMachines is the number of machines to be running in the cluster (16)
%%            2.NumOfProcessors is the number of processes to be running in the cluster (32)
%%            3.NumOfNodes is the matrix dimension (2^18)
%%	      4.initMat is set to initialize the matrix or not typically (0)
%%	      5.EdgesPerVertex is the number of edges per each vertex (105 for 2^18)
%%            4.max_iteration is the number of iterations for the Lanczos-SO algorithm (20)
%%            5.eig_k is the number of eigenvalues to be calculated (10)
%%            6.KeepB is to keep the same random B vector or not (1)
%%            7.Run_schedule is to scan the input table and divide the task evenly among processes (0) 
%%            8.StoreTFS is whether to store the input matrix into TFS (1)
%%            9.StoreLHD is whether to store the input matrix into local hard disk
%%            10.TFS is to run the for loop within TFS or LHD (1 for TFS, otherwise LHD)
%% Note 1: the main process can read the variables in m files.
%% Note 2: the parallel version should not delete the temporary table, it will mess up other processes' opertaions. So I move the delete temporary table in the main process.
%% Note 3: The inputmatrix will be set as automatically as 'M{NumOfNodes}' say M4096
% Note 4: The random vector B will be set as automatically 'B{NumOfNodes}' say B4096
%% Note 5: The first lz_q1 will be named as {NumOfNodes}lz_q1

%% Author: Yin Huang
%% Date: Mar, 15, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% Version by Mar 15: assumes that input matrix has been initialized. Vector B has not been initialized
%%

%%% Connect to the DB to access the global variables among multiple
%%% processes in the whole cluster
myDB;


%%%% Create a folder benchmark to store the debugging information
if ~exist('benchmark','dir')
        mkdir('benchmark');
end

fname = ('benchmark/stat.txt');
fstat = fopen(fname,'a+');

disp(['Start time: ' sprintf('\n')]);
StartTime = datestr(now);
fwrite(fstat,['***********************************************' sprintf('\n') 'Begin time: ' StartTime sprintf('\n*******************************************')]);
diary (['YinEigen: ' num2str(NumOfNodes) '_Machines' num2str(NumOfMachines) '_Proc' num2str(NumOfProcessors) '_Iter' num2str(max_iteration) '_logs.txt']);
fwrite(fstat,['**Commands: YinEigen( ' num2str(NumOfMachines) ',' num2str(NumOfProcessors) ',' num2str(NumOfNodes) ',' num2str(initMat) ',' num2str(EdgesPerVertex) ',' num2str(max_iteration) ',' num2str(eig_k) ',' num2str(KeepB) ',' num2str(Run_schedule) ',' num2str(StoreTFS) ',' num2str(StoreLHD)  ',' num2str(TFS) ')' ]);
disp(['TFS is ' num2str(TFS)]);
lz_allTime = tic;

machines=getMachines(NumOfMachines);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL variables need be accessed by all processors. I store them in the table every processor will read from the table.
NumOfMachines; % num of machines for computation
NumOfNodes; %nodes in the graph
max_iteration; % iteration times 

machines_t = DB('NumOfMachines');
nodes_t = DB('NumOfNodes');
proc_t=DB('NumOfProcessors');

cur_it = DB('cur_it');
alpha_t = DB('alpha'); %% store the alpha array in accumulo table 'alpha'
beta_t = DB('beta'); %% store the beta array in accumulo table 'beta'
parallel_sax_alpha_output = DB('alpha_sax_temp'); % delete temp tables in main note2
parallel_sax_beta_output = DB('beta_sax_temp');   % delete temp tables in main note2
norm_v_temp = DB(['lz_norm_v' num2str(NumOfNodes) '_temp']);
norm_b_temp = DB(['lz_norm_B' num2str(NumOfNodes) '_temp']);
so_rpath = DB('so_rpath');  %% selective orthogonalize intermidate output table
cur_loop_j = DB('cur_loop_j'); %% so inside loop identifier j every process need to know this value to computeR
rtv_temp = DB('rtv_temp'); %% so inside we need calculate the dotproduct of rtv, this table is used to save the temp result
so_rrtv = DB('so_rrtv'); %% so to store the vector 'rrtv' which is used to update lz_vpath, lz_vpath = lz_vpath - so_rrtv;
temp_lz_vpath = DB([num2str(NumOfNodes) 'lz_vpath']);
temp_mv_temp=DB('mv_temp');
temp_dot_temp=DB('dot_temp');


delete(alpha_t);
delete(beta_t);
delete(parallel_sax_alpha_output);
delete(parallel_sax_beta_output);
delete(norm_v_temp);
delete(norm_b_temp);
delete(so_rpath);
delete(cur_loop_j);
delete(rtv_temp);
delete(so_rrtv);
delete(temp_lz_vpath);
delete(temp_mv_temp);
delete(temp_dot_temp);

alpha_t = DB('alpha');
beta_t = DB('beta');
parallel_sax_alpha_output = DB('alpha_sax_temp');
parallel_sax_beta_output = DB('beta_sax_temp');
norm_v_temp = DB(['lz_norm_v' num2str(NumOfNodes) '_temp']);
norm_b_temp = DB(['lz_norm_B' num2str(NumOfNodes) '_temp']);
so_rpath = DB('so_rpath');
cur_loop_j = DB('cur_loop_j');
rtv_temp = DB('rtv_temp');
so_rrtv = DB('so_rrtv');
temp_dot_temp=DB('dot_temp');
temp_lz_vpath = DB([num2str(NumOfNodes) 'lz_vpath']);

m_assoc = Assoc('1,','1,',sprintf('%d,',NumOfMachines));
put(machines_t,m_assoc);
n_assoc = Assoc('1,','1,',sprintf('%d,',NumOfNodes));
put(nodes_t,n_assoc);
p_assoc = Assoc('1,','1,',sprintf('%d,',NumOfProcessors));
put(proc_t,p_assoc);


% 'scalar_b' is the norm of the random vector B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local variables to construct the Tridigonal matrix%%%%%%%%
alpha = zeros(1,max_iteration);
bet = zeros(1,max_iteration);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hard coded variables
v_prefix = [num2str(NumOfNodes) 'lz_q'];   %% v_prefix is lz_q to retrieve the tables named from lz_q{1:row}
q_path = cell(max_iteration+1,1);
scalar_b_path = 'scalar_b';
B_path = ['B' num2str(NumOfNodes)];

%%% initialize q_path array with the name lz_q{i}%%%%%%%%%%
for i = 1:max_iteration+1
	q_path{i} = [v_prefix num2str(i)];
end
for i = 2:max_iteration+1
	tempary = DB(q_path{i});
	delete(tempary);
end
for i = 1:max_iteration+1
        tempary =DB(q_path{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% assume InputMatrix has already been initilized and stored in 'InputMatrix' in my case yes it is the test case.%%%%%

disp([sprintf('\tRunning YinEigen with the following configuration:\n')]);
disp([num2str(NumOfMachines) 'machines:' machines sprintf('\n')]);
disp([num2str(NumOfProcessors) sprintf(' processors\t')]);
disp([num2str(NumOfNodes) sprintf(' nodes\t')])
disp([num2str(max_iteration) sprintf(' max iterations\t')]);
disp([num2str(eig_k) sprintf(' top eigen values')]);
disp(['KeepB yes?no: ' num2str(KeepB)]);

disp(['Run_schedule yes?no: ' num2str(Run_schedule)]);

disp(['StoreTFS yes?no: ' num2str(StoreTFS)]);
disp(['StoreLHD yes?no: ' num2str(StoreLHD)]);
disp(['TFS yes?no: ' num2str(TFS)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing the input matrix%


if (initMat == 1)
initM_edges = DB('edges');
put(initM_edges,Assoc('1,','1,',sprintf('%d,',EdgesPerVertex)));
disp(['Now initializing the input matrix in ' 'M' num2str(NumOfNodes)]);
tic;
eval(pRUN('SaveGraphData',NumOfProcessors/2,machines));
total_time = toc;
disp(['Total time to initialize M' num2str(NumOfNodes) ' is ' num2str(total_time)]);
fwrite(fstat,['Total time to initialize M' num2str(NumOfNodes) ' is ' num2str(total_time) ]);
end


%%% initialize random vector b stored in table 'B{NumOfNodes}'
disp(['Initializing the random vector b in table B' num2str(NumOfNodes)]);
% initB.m :
%	1.create the random vector B in table B{NumOfNodes}
%   2.calculate the norm of vector B in scalar_b table
%	3.save the normalized vector B in {NumOfNodes}lz_q1 

%% 
%%  Initialize vector B and store it into the RAM of each process
%%  Random vector B will be stored at B{NumOfNodes}
%%
if (KeepB ~= 1)
	this = tic;
	eval(pRUN('initB',NumOfProcessors,machines));
    that = toc(this);
	disp(['InitB takes: ' num2str(that)]);
	fwrite(fstat,['InitB takes: ' num2str(that) sprintf('\n')]);
    
    %% initialize the lz_q1 for the for loop
     eval(pRUN('parallel_lz_norm_B_p1',NumOfProcessors,machines));
    parallel_lz_norm_B_p2(NumOfProcessors);
    eval(pRUN('parallel_scalarmult_B',NumOfProcessors,machines));
    
    %this = tic;
    %eval(pRUN('storeBToTFS',NumOfProcessors,machines));
    %that = toc(this);
    %disp(['Store B to TFS takes: ' num2str(that)]);
    %fwrite(fstat,['StoreB takes: ' num2str(that) sprintf('\n')]);
end

%% TEST of initB passed.

%% Normalize B and save the scalar*B into Lz_q1 only need to do once.


 
%%%% Below is to schedule the tasks evenly among all processors %%%%%%%
%%%% The cut of input matrix will be stored in the table Cut{NumOfNodes}
%%%% For parallel part, each processor will read through cut('previous cut',:)+1 until cut('current cut',:)

scheduler = DB(['Cut' num2str(NumOfNodes)]);
totalentries = DB(['Entries' num2str(NumOfNodes)]);
if (Run_schedule == 1)
    delete(scheduler);
	delete(totalentries);
	scheduler = DB(['Cut' num2str(NumOfNodes)]);
	totalentries = DB(['Entries' num2str(NumOfNodes)]);
    this = tic;
    eval(pRUN('myscheduler',NumOfProcessors,machines));
    that =toc(this);
    disp(['Scan table running time: ' num2str(that) 's' sprintf('\n')]);
    fwrite(fstat,['Scan table running time: ' num2str(that) 's' sprintf('\n')]);

		this = tic;
		myscheduler_p2(NumOfMachines,NumOfProcessors);
		that = toc(this);
		disp(['Scheduler 2 running time: ' num2str(that) 's']);
		fwrite(fstat, ['Scheduler 2 running time: ' num2str(that) 's' sprintf('\n')]);
end


%%%%% Now we store the input matrices to corresponding work node %%%%
if(StoreTFS == 1)  
    this = tic;
    eval(pRUN('storeDataToTFS',NumOfProcessors,machines));
    that=toc(this);
    disp(['Store to TFS running time: ' num2str(that) 's']);
    fwrite(fstat, ['Store to TFS running time: ' num2str(that) 's' sprintf('\n')]);
end

if(StoreLHD == 1)
	this = tic;
	eval(pRUN('LHD_store',NumOfProcessors,machines));
	that = toc(this);
    disp(['Store to LHD running time: ' num2str(that) 's']);
    fwrite(fstat, ['Store to LHD running time: ' num2str(that) 's' sprintf('\n')]);
end


%%Now start the for loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = ('benchmark/stat.txt');
fstat = fopen(fname,'a+');

for it = 1:max_iteration
	thistic=tic;
	disp('**************myEigen iterations***********************');
	
	 %%%%%%%%%%%%%%%%%%%%%%  matrix * vector begin **********************
       %% update the cur_it first
        it_assoc = Assoc('1,','1,',sprintf('%d,',it));
        put(cur_it,it_assoc); %% globalize the current iteration so all processors will be able to read the right lz_q{it}
	 disp('Now saving vector to the global alluxio file');
	disp('First check if the file already exists from previous experiment so we delete it');
	system(['alluxio fs rmr /mytest/' num2str(it) 'v*']);
        this = tic;
        saveVectorToGTFS();
        saveT = toc(this);
        disp([ 'Saving vector to global file costs ' num2str(saveT) 's']);

        disp(['Now each machine makes its own copy of vector']);
        this = tic;
        eval(pRUN('saveVectorToTFS', NumOfMachines,machines));
        savelocal = toc(this);
        disp(['Machine copy vector costs ' num2str(savelocal) 's']);

        disp(['computing v=Aq ' num2str(it) ' ...']);
	temp = DB('mv_temp'); delete(temp);temp = DB('mv_temp');  %% remove the temp table from previous operation for paralell_mv_p1.m
  if(TFS == 1)
        disp('Running TFS version of matrix multipilcation');
        this = tic;
        %% Version 1: when vector not stored in Alluxio
        %% eval(pRUN('Alluxio_Row_mv',NumOfProcessors,machines));
        %% Version 2: when vector saved in Alluxio
        eval(pRUN('Alluxio_Row_mv_version2',NumOfProcessors,machines));
        that = toc(this);
        fstat = fopen(fname,'a+');
        disp(['Iteration ' num2str(it) ' Alluxio_Row_mv_version2 takes: '  num2str(that)]);
        fwrite(fstat,['Iteration ' num2str(it) ' Alluxio_Row_mv_version2 takes: '  num2str(that) sprintf('\n')]);
  else
        disp(['Running the local disk version of matrix*vector']);
         this = tic;
        eval(pRUN('LHD_Row_mv',NumOfProcessors,machines));
        that = toc(this);
        disp(['Iteration ' num2str(it) ' LHD_Row_mv takes: '  num2str(that)]);
        fwrite(fstat,['Iteration ' num2str(it) ' LHD_Row_mv takes: '  num2str(that) sprintf('\n')]);

        end


        disp(['Result of v = Aq ' num2str(it) ' is saved in table ' num2str(NumOfNodes) 'lz_vpath']);

    
    %%%%%%%%%%%%%  matrix * vector done! ***************************************   

	%% alpha(it) = dotproduct(aqnpath,q_path{it}, NumOfNodes, NumOfMachines);
	%% dotproduct should save the output in table 'dot_output' ('1,','1,'), also result could be read from dot_product;
	%% alpha(it) will read from the dot_output table.

    % num2str(NumOfNodes)
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vi * v begin **********************************************
    %%  Begin Version 2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	disp('Now Prepare for vi*v storing the elemetns into the tfs');
	disp('Removing existing file in tfs if exists');
	%% lz_vpath is located as /mytest/1048576lz_vpath_global
	system(['alluxio fs rmr /mytest/' num2str(NumOfNodes) 'lz_vpath*']);
	
	str = ['Now saving lz_vpath to global file ... ' sprintf('\n')];
	disp(str); fwrite(fstat, str);
	this = tic;
	saveVpathToGTFS();  
	saveVtime = toc(this);
	
	str = ['Saving lz_vpath to global file takes: ' num2str(saveVtime) 's' sprintf('\n')];
	disp(str); fwrite(fstat, str);
	
	str = ['Saving lz_vpath to local machines ...' sprintf('\n')];
	disp(str); fwrite(fstat, str);
	this = tic;
	eval(pRUN('saveVpathToTFS', NumOfMachines,machines));
	saveLocalV = toc(this);
	str = ['Saving lz_vpath to local machines takes ' num2str(saveLocalV) 's' sprintf('\n')];
	disp(str); fwrite(fstat, str);

	disp(['Computing dotproduct of vi * v ... and saving the result in alpha(' num2str(it) ')']);

        parallel_dotproduct_p1_temp = DB('dot_temp');
        this = tic;
        eval(pRUN('parallel_dotproduct_p1_version2',NumOfProcessors,machines));
        alpha(it) = parallel_dotproduct_p2();
        that = toc(this);
        disp(['Iteration ' num2str(it) ' Calculating alpha takes: '  num2str(that)]);
        fwrite(fstat,['Iteration ' num2str(it) ' Calculating alpha takes: '  num2str(that) sprintf('\n')]);
        delete(parallel_dotproduct_p1_temp);
        disp('Saving alpha to alpha_t');
        alpha_temp_Assoc = Assoc(sprintf('%d,',it),'1,',sprintf('%.15f,',alpha(it)));
        put(alpha_t, alpha_temp_Assoc);
        disp(['Result of alpha[' num2str(it) '] =' num2str(alpha(it)) ' is saved.']);

		
%%	save vi vector to Global file       
%%	saveViToGTFS()
%%	saveViToTFS()

%%      save lz_vpath to Global file. first need remove the previous one 
%% 	remove_previous_lz_vpath
%%	saveVToGTFS()
%% 	saveVToTFS()
%%	

%%      parallel_dotproduct_p1_version2()
%%	parallel_dotproduct_p2_version2() %% should be the same as version 1. 
	
	%%%%%%%%%
	%% Below is 1st version of vi * v
	%{
	disp(['Computing dotproduct of vi * v ... and saving the result in alpha(' num2str(it) ')']);
	
	parallel_dotproduct_p1_temp = DB('dot_temp');
	this = tic;
	eval(pRUN('parallel_dotproduct_p1',NumOfProcessors,machines));
	alpha(it) = parallel_dotproduct_p2();
	that = toc(this);
	disp(['Iteration ' num2str(it) ' Calculating alpha takes: '  num2str(that)]); 
	fwrite(fstat,['Iteration ' num2str(it) ' Calculating alpha takes: '  num2str(that) sprintf('\n')]);
	delete(parallel_dotproduct_p1_temp);
	disp('Saving alpha to alpha_t');
	alpha_temp_Assoc = Assoc(sprintf('%d,',it),'1,',sprintf('%.15f,',alpha(it)));
	put(alpha_t, alpha_temp_Assoc);
	disp(['Result of alpha[' num2str(it) '] =' num2str(alpha(it)) ' is saved.']);
	
       %}
	%%%%%%%%%%%%%%%%%%% vi * v done! *****************************************************

    %%%%%%%%%%%%%%%%%%% Calculating v = v - beta{i-1}*v{i-1} - alpha{i}*v{i} **********************
	
	 this = tic;
        eval(pRUN('onetime_saxv',NumOfProcessors,machines));
        that = toc(this);
        disp(['Iteration ' num2str(it) ' onetime_saxv: '  num2str(that)]);
        fwrite(fstat,['Iteration ' num2str(it) ' onetime_saxv: '  num2str(that) sprintf('\n')]);

	disp(['v is saved in ' num2str(NumOfNodes) 'lz_vpath table']);
     

	%%%%%%%%%%%%%%%%%%% Calculating v = v - beta{i-1}*v{i-1} - alpha{i}*v{i}  Done!**********************

	%*************  Calculating beta{i} = ||v|| *************************************************************

    disp(['Computing beta[' num2str(it) ']...']);
	parallel_lz_norm_v_tempt = DB(['lz_norm_v' num2str(NumOfNodes) '_temp']);
	this = tic;
	eval(pRUN('parallel_lz_norm_v_p1',NumOfProcessors,machines));
	parallel_lz_norm_v_p2; %% scalar_v is written to beta_i in the table beta_t('i,','1,')
	that = toc(this);
	disp(['Iteration ' num2str(it) ' beta takes: '  num2str(that)]);
	fwrite(fstat,['Iteration ' num2str(it) ' beta takes: '  num2str(that) sprintf('\n')]);	
	bet(it) = scalar_v;
	delete(parallel_lz_norm_v_tempt);
	disp(['beta[' num2str(it) '] = ' num2str(bet(it))]);
	
	%******** Calculating beta{i} Done locally %%%%%%%%%%%%%%%%%%%


	disp(['Constructing the Tridigonal matrix...']);
	num_ortho = 0;
	tempTmatrix = constructT(it, alpha, bet); 
	[Q,D] = eig(tempTmatrix);
        D = diag(D);
	%% Do selective_orthogonalize locally%%%%
	v_path = ([num2str(NumOfNodes) 'lz_vpath']);
	disp(['NumOfMachines in SO: ' num2str(NumOfMachines) 'Starting so, iterations # is ' num2str(it) ' beta_it value is: ' num2str(bet(it))]);
	this = tic;
	num_ortho = parallel_selective_orthogonalize(it, bet(it), v_path, Q,D, NumOfNodes, NumOfMachines,NumOfProcessors);
	that = toc(this);
	disp(['Iteration ' num2str(it) ' SO takes: '  num2str(that)]);
	fwrite(fstat,['Iteration ' num2str(it) ' SO takes: '  num2str(that) sprintf('\n')]);
	disp(['Number of orthongalization: ' num2str(num_ortho)]);

   	

	if(num_ortho > 0)
		disp(['Recomputing beta[' num2str(it) ']']);
		 parallel_lz_norm_v_tempt = DB(['lz_norm_v' num2str(NumOfNodes) '_temp']);
		eval(pRUN('parallel_lz_norm_v_p1',NumOfProcessors,machines));
		parallel_lz_norm_v_p2; %% scalar_v is written to beta_i in the table beta_t('i,','1,')
		bet(it) = scalar_v;
		disp(['beta[' num2str(it) ']=' num2str(bet(it))]);
		delete(parallel_lz_norm_v_tempt);
	end

	if(num_ortho > it - 1)
	disp('The new vector converged. Finishing ...');
        compute_eigval(it, alpha, bet, eig_k);
	save_tridiagonal_matrix(alpha, bet, it);
	break
	end 
	if(bet(it) == 0.0)
	disp(['beta[' num2str(it) ']=0. finishing']);
	disp('Saving the tridiagonal matrix');
	compute_eigval(it, alpha, bet, eig_k);
	save_tridiagonal_matrix(alpha, bet, it);
	break
	end
	
	disp(['Computing q' num2str(it+1) '...']);

	%%%%%%%%%%%%%%  Update {NumOfNodes}lz_vpath %%%%%%%%%%%%%%%%%%%%%%
	this = tic;
	eval(pRUN('parallel_update_q',NumOfProcessors,machines));
	that = toc(this);
	disp(['Iteration ' num2str(it) ' Update Q takes: '  num2str(that)]);
	fwrite(fstat,['Iteration ' num2str(it) ' Update Q takes: '  num2str(that) sprintf('\n')]);
	disp(['q_{' num2str(it+1) '} is calcualted']);

oneIterationTime=toc(thistic);
    disp(['Iteration: ' num2str(it) ': ' num2str(oneIterationTime) 's']);
        fwrite(fstat,['Iteration: ' num2str(it) ': ' num2str(oneIterationTime) 's' sprintf('\n')]);

	compute_eigval(it, alpha, bet, eig_k);
	disp('Saving the tridiagonal matrix');
	save_tridiagonal_matrix(alpha, bet, it);
	
    %oneIterationTime=toc(thistic);
   % disp(['Iteration: ' num2str(it) ': ' num2str(oneIterationTime) 's']);
%	fwrite(fstat,['Iteration: ' num2str(it) ': ' num2str(oneIterationTime) 's' sprintf('\n')]);
end  %% end for loop
%}	
	if ( TFS ~= 1) %% LHD mode we remove the local disk files to save space.
    eval(pRUN('deletefiles',NumOfMachines,machines));
    end
	disp('!!!!!!Reached the max iterations. Finishing...');
	
	disp('Summarizing alpha[] and bet[]...');
	disp(sprintf('\n\talpha\tbeta'));
	for n = 1:max_iteration
	disp([num2str(n) sprintf('\t') num2str(alpha(n)) sprintf('\t\t') num2str(bet(n))]);
	end
	
	alltime = toc(lz_allTime);
	disp(['Total running time is: ' num2str(alltime)]);
	  endtime = datestr(now);

	disp(['Ending time: ' endtime  sprintf('\n')]);
	
	fwrite(fstat,['Ending time: ' endtime  sprintf('\n')]);
	disp(['Begin time: ' StartTime]);
	diary off;
	fclose(fstat);
end %end function

