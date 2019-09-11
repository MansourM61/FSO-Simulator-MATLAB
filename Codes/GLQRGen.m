function [ zerosave,w,reccalls ] = GLQRGen( nmax,tol,numpoint,reccalls,maxrec)

%this function tested to work for N < 300, larger values suffer from
%underflow issues in the weights and exponential.

%this function generates a gauss Laguerre quadrature rule for a given value n, it
%does so by generating via reccursion a table of the Laguerre functions 
% (i.e. a Laguerre polynomial times the root of a weight function) in
%the positive domain and then looking for zeros via a change of sign and 
%interpolating (quadratically) between the two values (or finding an exact zero)
%Using more points will give better accuracy, a warning is output if zeros
%are missed due to lack of points
%The function will then focus in on the zeros until a tolerance is reached,
%this tolerance is based on absolute uncertainty in each root

%note such rules are used to approximate integrals of the form
% int_{0}^{inf} dx f(x) exp(-x) to  sum_{i=1}^N w_i*f(x_i)

%recalls follows how many times the function has called itself

%Author David Holdaway: year 2012 
%% check inputs
if nargin==1
    tol(1) = 1.0e-12; %default
end
if nargin <4
    reccalls = 0;
end
if nargin <5
    maxrec = 20; %maximum number of recalls to achieve convergence
end

if tol(1) < 5*eps && reccalls ==0
    warning('this position tolerance is of order eps, function may fail'); %#ok<*WNTAG>
end

if nargin < 3
    numpoint = 200*(nmax+5); %usually an adequate amount 
end

if nmax ==0
   zerosave=nan; w=nan; return
end
if nmax ==1
   zerosave=0; w=1; return
end
if length(numpoint) == 1 %initial course zero findings 
x = linspace(0,2*sqrt(nmax+1),numpoint); %needs lots of points for resolution
x = x.^2; %want increasing spacing towards inf, reduced towards zero
else  %previous zeros fed back in
    x = numpoint;    
end
%tic
%% generate Laguerre function table
ltable = zeros(3,length(x)); %need to store three to generate
ltable(1,:) = exp(-x/2);
underflowchk = ltable(1,:)==0;
if sum(underflowchk)~=0
    warning('some x values cut off due to underflow, may still be fine or may cut off zeros')
    %this is a problem for large N
    x = x(ltable(1,:)~=0); %cut x
    ltable = zeros(3,length(x));
    ltable(1,:) = exp(-x/2);
end
ltable(2,:) = (1-x).*ltable(1,:);

for k = 2:nmax
    %L_{n} = (1/n) [ (2n-1 - x) L_{n-1} -  (n-1) L_{n-2} ]

    ltable(3,:) = ((2*k-1 - x).*ltable(2,:) -(k-1).*ltable(1,:))/k;
   ltable([1,2,3],:) =  ltable([2,3,1],:); %permute to use again
   
    
end

lfin = ltable(2,:);

% if reccalls == 0
%    figure
%     plot(sqrt(x),lfin )
%   hold on
% end

clear ltable
%% find zeros only in +ve

zerosave = zeros(1,nmax); 
zcount = 0; zerosave2 = zerosave;
rng = 2:length(x); 

for lp = rng
    if sign(lfin(lp)) == 0 %found "exact" zero!
        zcount = zcount + 1;
       zerosave(zcount) = x(lp);
       zerosave2(zcount) = lp; 
      %skip loop an extra one
    elseif sign(lfin(lp-1)) == 0
        %skip this as it will screw up next logic
    elseif sign(lfin(lp-1)) ~= sign(lfin(lp))
        %change of sign implies zero
    zcount = zcount + 1;
    %linearly interpolate between two x points to approximate zero
    xpsave = x(lp)-x(lp-1); %spacing
    root = xpsave*abs(lfin(lp-1))/abs(lfin(lp)-lfin(lp-1));

    zerosave(zcount) = x(lp-1) + root;
    zerosave2(zcount) = lp; 
    end
    
end
clear hgrad hfin 



if zcount < nmax
    
   warning('missed zeros, doubling number of points')
    if reccalls == 0
    x = linspace(0,nmax^2*20+4*(1+nmax),numpoint*2); 
    x = sqrt(x);
    [ zerosave,w,reccalls ] = GLQRGen( nmax,tol,x,reccalls,maxrec );
    return
    end %call again with double number of points
     warning('failed: not enough zeros found use more points or nmax too large')
    
    w=nan; zerosave=zcount ;
    return %fail
end
if zcount > nmax
    warning('failed to distinguish zeros, tolerance likely set low, removing extra zeros')
    %zcount 
tmp  = diff(sqrt(zerosave)) > 1e-5; %find repeated zeros.. low x has too much resolution for dp
if sum(tmp) ==nmax;
    
    zerosave = zerosave(tmp); %remove repeated
    zerosave2 = zerosave2(tmp); %remove repeated
else
    zerosave=zcount; w = nan; return %fail
    
end
end

        xpsave = zerosave2;
        for k=1:length(zerosave2)
           xpsave(k) = x(zerosave2(k)) - x(zerosave2(k)-1);            
        end
         converged = xpsave <= tol;
         notconv = logical(1-converged);
         convcond = sum(notconv) ==0;
        %spacings vary

clear x %no longer required

%toc
%% used to refine results by feeding back zero estimate

if convcond==false %feed back x to function, create refined x range
if reccalls < maxrec

%     finex = kron(zerosave,ones(1,9));
%     finex = sort([0,finex + kron(xpsave.*ones(1,length(zerosave)),...
%        [-1,-1/4,-1/16,-1/64,0,1/64,1/16,1/4,1])]);
%    
%     finex = finex( diff(finex)>eps); %if these points are not distinct, remove them
    
    tmp1 = kron(zerosave(notconv),ones(1,9)) + ...
            kron(xpsave(notconv).*ones(1,sum(notconv)),...
            [-1,-1/4,-1/16,-1/32,0,1/32,1/16,1/4,1]);
    tmp2 =  kron(zerosave(converged),ones(1,3)) + ...
            kron(xpsave(converged).*ones(1,sum(converged)),...
            [-1,0,1]); %converged so don't both refining much   
    finex = sort([0,tmp1,tmp2]);
    clear tmp1 tmp2

    reccalls = reccalls+1;

   [ zerosave,w,reccalls ] = GLQRGen( nmax,tol,finex,reccalls,maxrec ); 
   return
else

    warning('failed to converge to desired tolerence, after max reccursions');
    disp(max(xpsave));

end
end
%% find weights if zeros to tolerance
% find weights, regenerate legtable using zero positions only and no
% exp(-x/2) factor

%w_i = x_i/(n+1)^2(L_{n+1}(x_i))^2

ltable = zeros(nmax+2,length(zerosave));
ltable(1,:) = ones(size(zerosave));
ltable(2,:) = (1-zerosave).*ltable(1,:);

for k = 2:nmax+1 %up to n+1th, note again indexes start from zero
        %L_{n} = (1/n) [ (2n-1 - x) L_{n-1} -  (n-1) L_{n-2} ]
    ltable(k+1,:) = ((2*k-1 - zerosave).*ltable(k,:) -(k-1).*ltable(k-1,:))/k;
    
end

w = zerosave;

for k = 1:length(zerosave)

 w(k) = zerosave(k)./ltable(nmax+2,k).^2;
    
end

w = w/(nmax+1)^2;

% format longe
% for k =1:20
% diffsave(k)=(sum(x.^k.*w)-factorial(k))/factorial(k);
% end