%%*********************************************************
%% smat: compute the matrix smat(x).
%%
%%   M = smat(blk,x,isspM); 
%%
%% SDPNAL: 
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh 
%%**********************************************************

function M = smat(blk,xvec,isspM)
   if (nargin < 3); isspM = zeros(size(blk,1),1); end 
   if ~iscell(xvec) 
      if strcmp(blk{1},'s')
         M = mexsmat(blk,xvec,isspM);      
      else
         M = xvec; 
      end   
   else 
      M = cell(size(blk,1),1);     
      if (length(isspM)==1)
         isspM = isspM*ones(size(blk,1),1); 
      end
      for p=1:size(blk,1)
         pblk = blk(p,:);
         if strcmp(pblk{1},'s')
            M{p} = mexsmat(pblk,xvec{p},isspM(p));
         else
            M{p} = xvec{p}; 
         end   
      end
   end
end
%%*********************************************************

