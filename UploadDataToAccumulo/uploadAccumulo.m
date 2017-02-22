function uploadAccumulo()
%% Usage:
%%	This function is used to upload data from *.edge into Accumulo table
%%	The reason is because parallel uploading data into accumulo stops in middle 
%%	This function will read one file at one time to see if it proceeds to upload data
%%
%% input: 
%%	1. the location of where .edge files located 
%%   	2. 
%%%
%% Author: Yin Huang
%% Date: Apr 23, 2016

%%connect to DB
myDB;

out = DB('TEST');
%% 

contents = dir('*.edge'); % only look for .edge filesnum
for i = 1:numel(contents)
  filename = contents(i).name;
  fileID = fopen(filename,'r');
  formatSpec = '%d %d';
  sizeA = [2 Inf];
  A = fscanf(fileID,formatSpec,sizeA);
  A = A';
  %disp(['First vector is:' sprintf('\n')])   
  %A(:,1);
  %disp(['Second vector is:' sprintf('\n')])	
  %A(:,2);
  Arow = sprintf('%d,',A(:,1));
  Acol = sprintf('%d,',A(:,2));
  AssocArr = Assoc(Arow, Acol, '1,');
  clear A;
  disp('Now uploading');
  put(out, AssocArr);
  disp('One file has been uploaded');
end
