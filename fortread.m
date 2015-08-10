function array = fortread(fid, M,N)
%
% Read fortran binary files usage is
% array = fortread(fid, M,N)
% where 
% M is the row size of the array
% N is the column size of the array
% fid is the fileid that is supposed to have been opened as binary.
% fid = fopen(filename,'r','b');
%
 [Bogus,count] = fread(fid,1,'int32');
 array = fread(fid,[M,N],'double','s');
 [Bogus,count] = fread(fid,1,'int32');
