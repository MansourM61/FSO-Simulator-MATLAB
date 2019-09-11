function [ x,w,reccalls ] = GHQRGen( nmax,tol,numpoint,reccalls,maxrec)

%this function generates a gaussian quadrature rule for a given value n, it
%does so by generating via reccursion a table of the Hermite functions in
%the positive domain and then looking for zeros via a change of sign and
%interpolating (quadratically) between the two values (or finding an exact zero)
%Using more points will give better accuracy, a warning is output if zeros
%are missed due to lack of points
%The function will then focus in on the zeros until a tolerance is reached,
%this tolerance is two part (though only 1 is needed)
% 1) absolute uncertainty in each position
% 2) measured by how well the zeros satisfy
%       sum_{k \ne m) 1/(x_k - x_m) = x_m
% to just use condition 2 use something very high for tol(1) e.g. [1,1e-14]
% as an input

%note such rules are used to approximate integrals of the form
% int_{-inf}^{inf} dx f(x) exp(-x^2) to  sum_{i-1}^N w_i*f(x_i)
%which are exact if f(x) is a polynomial of degree < 2N-1

%recalls follows how many times the function has called itself

%Author David Holdaway: year 2012
%% check inputs
if nargin==1
    tol(1) = 1.0e-12; %default
    tol(2) = nmax*1.0e-12; %default
end
if nargin <4
    reccalls = 0;
end
if nargin <5
    maxrec = 5; %maximum number of recalls to achieve convergence
end

if tol(1) < 5*eps && reccalls ==0
    warning('this position tolerance is of order eps, function may fail') %#ok<WNTAG>
end
if length(tol)==1
    tol = [tol,1];%no condition on tol 2
end
if tol(2) < nmax*eps && reccalls ==0
    warning('second convergence condition unlikely to be achieved') %#ok<WNTAG>
end
if nargin < 3
    numpoint = 50*(nmax+5); %usually an adequate amount
end

if nmax ==0
    x=nan; w=nan; return
end
if nmax ==1
    x=0; w=1; return
end
if length(numpoint) == 1 %initial course zero findings
    x = linspace(0,sqrt(nmax)*1.5+5,numpoint); %needs lots of points for resolution
else  %previous zeros fed back in
    x = numpoint;
end
%tic
%% generate Hermite function table
htable = zeros(3,length(x)); %need to store three to generate
htable(1,:) = exp(-x.^2/2)/pi^(1/4);
underflowchk = htable(1,:)==0;
if sum(underflowchk)~=0
    warning('some x values cut off due to underflow of psi_0(x), may still be fine or may cut off zeros') %#ok<WNTAG>
    %this is a problem for large N
    x = x(htable(1,:)~=0); %cut x
    htable = zeros(3,length(x));
    htable(1,:) = exp(-x.^2/2)/pi^(1/4);
end
htable(2,:) = x.*htable(1,:)*sqrt(2);

for k = 2:nmax
    %phi_{k}(x) = sqrt(2/k) x phi_{k-1}(x) - sqrt((k-1)/k) x phi_{k-2}(x)
    
    htable(3,:) = x.*sqrt(2/k).*htable(2,:) ...
        -sqrt((k-1)/k).*htable(1,:);
    htable([1,2,3],:) =  htable([2,3,1],:); %permute to use again
    
    
end

hfin = htable(2,:);
hgrad = htable(1,:).*sqrt(2*nmax) - x.*htable(2,:);
%  figure
%   plot(x,hfin )
%   hold on
%  plot(x,hgrad,'r')
%psi'_n = sqrt(2*n) psi_{n-1} - x psi_{n}
clear htable
%% find zeros only in +ve

zerosave = zeros(ceil(nmax/2),1);
zcount = 0; zerosave2 = zerosave;
rng = 2:length(x); %idvalchk = diff(hfin)~=0;
%when errors of order eps some of the fine x grids can give the same value

if hfin(1) == 0 %for odd
    zerosave(1)=0;
    zcount =1;
end

for lp = rng
    if sign(hfin(lp)) == 0 %found "exact" zero!
        zcount = zcount + 1;
        zerosave(zcount) = x(lp);
        zerosave2(zcount) = lp;
        %skip loop an extra one
    elseif sign(hfin(lp-1)) == 0
        %skip this as it will screw up next logic
    elseif sign(hfin(lp-1)) ~= sign(hfin(lp))
        %change of sign implies zero
        zcount = zcount + 1;
        %interpolate between two x points to approximate zero via cubic spline
        % f=a(x-x_1)^3+b(x-x_1)^2+c(x-x_1)+d  x_1 being left point x(lp-1)
        xpsave = x(lp)-x(lp-1); %spacing
        d = hfin(lp-1);
        c = hgrad(lp-1);
        b = 3*( hfin(lp)-d)/xpsave^2 - (hgrad(lp)+2*c)/xpsave;
        a = (hgrad(lp)-c)/(3*xpsave.^2)-2*b/(3*xpsave);
        term1 = 2*b^3 - 9*a*b*c +27*a^2*d;
        term2 = sqrt((2*b^3-9*a*b*c+27*a^2*d)^2-4*(b^2-3*a*c)^3);
        %     root= [-1/(3*a),(1+1i*sqrt(3))/(6*a),(1-1i*sqrt(3))/(6*a)]*((term1+term2)/2)^(1/3);
        %     root= root + [-1/(3*a),(1-1i*sqrt(3))/(6*a),(1+1i*sqrt(3))/(6*a)]*((term1-term2)/2)^(1/3);
        %     root= root -b/(3*a);
        %
        %     gdroot = (isreal(root) & root>0 & root < xpsave)
        %     root = root(gdroot);
        %middle root is always the one required
        root= (1+1i*sqrt(3))/(6*a)*((term1+term2)/2)^(1/3)...
            + (1-1i*sqrt(3))/(6*a)*((term1-term2)/2)^(1/3) - b/(3*a);
        
        if (isreal(root)~=1 || root<0 || root > xpsave)
            %this has failed just linear interp
            root = xpsave*abs(hfin(lp-1))/abs(hfin(lp)-hfin(lp-1));
        end
        
        zerosave(zcount) = x(lp-1) + root;
        zerosave2(zcount) = lp;
    end
    
end
clear hgrad hfin

if zcount < ceil(nmax/2)
    warning('not enough zeros found use more points or nmax too large') %#ok<WNTAG>
    disp(zcount);
    w=nan;
    return
end
if zcount > ceil(nmax/2)
    warning('failed to distinguish zeros for some reason') %#ok<WNTAG>
    disp(zcount);
    w=nan;
    return
end

if reccalls == 0;
    xpsave = abs(x(2) - x(1)); %all spacings the same
    convcond1 = xpsave<tol(1);
else
    xpsave = zerosave2';
    if zerosave(1) ~= 0 %even
        xpsave(1) = x(zerosave2(1)) - x(zerosave2(1)-1);
    end
    for k=2:length(zerosave2)
        xpsave(k) = x(zerosave2(k)) - x(zerosave2(k)-1);
    end
    convcond1 = max(xpsave)<tol(1);
    %spacings vary
end
x = zerosave';
if x(1) ==0 %odd
    x = [-x(length(x):-1:2),x];
else %even
    x = [-x(length(x):-1:1),x];
end
%toc
%% used to refine results
%test using the identity sum_{k \ne m) 1/(x_k - x_m) = x_m and
if length(tol) == 2
    Ftest=x;
    for k = 1:nmax
        Ftest(k) = sum(1./(x(k)-x([1:k-1,k+1:nmax])))-x(k);
        %since x can only be accurate to eps, errors will propogate here
        %could instead look at convergence as uncertainty in x
    end
    convcond2 = max(Ftest) < tol(2);
else
    convcond2 =true;
end


if convcond2==false || convcond1==false %feed back x to function, create refined x range
    if reccalls < maxrec
        if zerosave(1) == 0
            finex = [0,kron(zerosave(2:length(zerosave))',ones(1,11))];
            if length(xpsave)==1
                finex = sort(finex + [0,kron(xpsave.*ones(1,length(zerosave)-1),...
                    [-1,-1/4,-1/16,-1/64,-1/256,0,1/256,1/64,1/16,1/4,1])]);
            else
                finex = sort(finex + [0,kron(xpsave(2:length(xpsave)).*ones(1,length(zerosave)-1),...
                    [-1,-1/4,-1/16,-1/64,-1/256,0,1/256,1/64,1/16,1/4,1])]);
            end
        else
            finex = kron(zerosave',ones(1,11));
            finex = sort([0,finex + kron(xpsave.*ones(1,length(zerosave)),...
                [-1,-1/4,-1/16,-1/64,-1/256,0,1/256,1/64,1/16,1/4,1])]);
        end
        
        reccalls = reccalls+1;
        
        [ x,w,reccalls ] = GHQRGen( nmax,tol,finex,reccalls,maxrec );
        return
    else
        
        warning('failed to converge to desired tolerence, after max reccursions') %#ok<WNTAG>
        bestachieved1 = max(xpsave);
        disp(bestachieved1);
        if length(tol)==2
            bestachieved2 = max(Ftest);
            disp(bestachieved2);
        end
    end
end
%% find weights if zeros to tolerance
% find weights, regenerate hermtable using zero positions only
%  and hermite polys NOT functions but with same norm factor, i.e.
%  psi_k(x)/exp(-x^2/2) and only at zeros of n + 1 th function

%w_i = 1/n(H_{n-1}(x_i))^2

htable = zeros(nmax,length(x));
htable(1,:) = ones(size(x))/pi^(1/4);
htable(2,:) = x.*htable(1,:)*sqrt(2);

for k = 2:nmax-1 %up to n-1th
    %phi_{k}(x) = sqrt(2/k) x phi_{k-1}(x) - sqrt((k-1)/k) x phi_{k-2}(x)
    
    htable(k+1,:) = x.*sqrt(2/k).*htable(k,:) ...
        -sqrt((k-1)/k).*htable(k-1,:);
    
end

w = x; %weights equal 2^{n-1} n! sqrt(pi) /n^2 H_{n-1}(x_i)^2
% with normed polys  = 1 / n Hnorm_{n-1}(x_i)^2
for k = 1:length(x)
    
    w(k) = 1./htable(nmax,k).^2;
    
end

w = w/(nmax);
