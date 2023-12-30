function symbaviol(h,tmax)
global GNEWT CASE MSCRIPT YEAR
% Good parameters: h = 0.03, t = 1e3 (close encounter before t~300), M=3?
% Bulirsch--Stoer section did a part on centering of mass-- check
CASE = 1; %Mercurius is always WHD
% GNEWT = 39.4845;
GNEWT = 0.00029591220823221284; %gaussian
% YEAR = 365.25;
YEAR = 1;
MSCRIPT = 2; %WH choice is 1, RT is otherwise for WHJ

str = sprintf("data/fcons.txt");
fcons = fopen(str,'w');
t = 0;
[m,x,v] = insolardunc(); 
[Psum0] = calcm(v,m);
[Q,junk] = convertcart(m,x,v);
h0 = h;
dE = 0;

[rlev,hlev,steps,hsub] = initv(h);
levv = zeros(size(rlev,2),1);
loop = 0;
[levc,vfunc] = calclevv(Q,v,m,rlev,h);
% [levc] = calclevsimple(Q,v,m,rlev,hlev(1));
% levc
[E0,L0] = consqv(m,Q,v);   
nstep = 0;
fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
    t,E0,norm(Q(:,2)-Q(:,3)),norm(Q(:,2)-Q(:,4)),norm(Q(:,2)-Q(:,5)),norm(Q(:,3)-Q(:,4)),norm(Q(:,3)-Q(:,5)),norm(Q(:,4)-Q(:,5)),levc(2,3),levc(2,4),levc(2,5),levc(3,4),levc(3,5),levc(4,5),0);
% fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
%     t,E0,norm(Q(:,3)-Q(:,4)),0,0,0,0,0,0,0,0,0,0);
t1 = cputime;
inc=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(t < tmax)
    t/tmax
    nstep = nstep+1;
    %before step 9415 redoing
%     if(nstep == 9415+0)
%         break
%     end
%     [Q,v] = mapv(Q,v,h,m);
%     [Q,v,levco] = SyMBAand(Q,v,m,levco,rlev,hlev,steps);
%     [Q,v,levc] = SyMBAsmooth(Q,v,m,rlev,hlev,hsub,levmax);
%     [Q,v,levc] = SyMBAsmoothfull(Q,v,m,rlev,hlev,hsub,rhill0);
%     [Q,v,levc,f1,vfunc,finc] = SyMBAdiscrete(Q,v,m,rlev,hlev,hsub,levc);
% [Q,v,levc,loop] = globalstepirr(Q,v,m,hsub,levc,1,rlev,hlev);
[Q,v,levc,loop] = globalstep(Q,v,m,hsub,levc,1,rlev,hlev);
%     [Q,v,levc,loop] = globalstepp(Q,v,m,hsub,levc,1,rlev,hlev);
%     [Q,v,levc,levv,levu] = adstep(Q,v,m,hsub,levc,rlev,hlev,levv);
%     levc
%     blergh
    t = t+h;
%     t = t + hlev(levu);
    [E,L] = consqv(m,Q,v);
    dE = (E - E0)/E0;
% fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
%     t,E,norm(Q(:,3)-Q(:,4)),0,0,0,0,0,0,0,0,0,0);
fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
    t,E,norm(Q(:,2)-Q(:,3)),norm(Q(:,2)-Q(:,4)),norm(Q(:,2)-Q(:,5)),norm(Q(:,3)-Q(:,4)),norm(Q(:,3)-Q(:,5)),norm(Q(:,4)-Q(:,5)),levc(2,3),levc(2,4),levc(2,5),levc(3,4),levc(3,5),levc(4,5),loop);

% display('step')
% break
end
t2 = cputime;
dt = t2-t1;
dE = (E - E0)/E0;
dL = (L - L0)./L0;
dL(3)
dE
dt
nstep
fclose('all');

    

end









% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Start SyMBA algorithm



function [Q,v] = SyMBA(Q,v,m,levcn,rlev,hlev,steps)
n = size(m,2);
rlevs = size(rlev,2);
p0 = calcm(v,m);
[Q] = mapSun(Q,v,hlev(1)/2,m);
for i=1:rlevs
    for j=2:n 
        for k=j+1:n
            if(levcn == i) %do steps at the level required, start with three-body problem
                stepsn = steps(levcn);
                hlevn = hlev(levcn);
                for l=1:stepsn
                    [v] = interpair(Q,v,hlevn/2,m,j,k); %combine interaction Kepler steps
                    [Q,v] = Keppair(Q,v,hlevn,m,j,k);
                    [v] = interpair(Q,v,hlevn/2,m,j,k);
                end
            end
        end
    end
end
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = adjustSun(Q,v,m,p0);

end





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Start SyMBA algorithm



function [Q,v] = SyMBAsimple(Q,v,m,hcurrent)
global GNEWT

n = size(m,2);
p0 = calcm(v,m);
[Q] = mapSun(Q,v,hcurrent/2,m);
for j=2:n 
    for k=j+1:n
        [v] = interpair(Q,v,hcurrent/2,m,j,k);
    end
end


for j=2:n
    gm = GNEWT*m(1);
    a = Q(:,j);
    b = v(:,j);
    [a,b] = kepler_stepxv(gm,a,b,hcurrent);
    v(:,j) = b(:);
    Q(:,j) = a(:);
end

for j=2:n 
    for k=j+1:n
        [v] = interpair(Q,v,hcurrent/2,m,j,k);
    end
end
        
[Q] = mapSun(Q,v,hcurrent/2,m);
[v] = adjustSun(Q,v,m,p0);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



function [Q,v,levc] = SyMBAand(Q,v,m,levc,rlev,hlev,steps)
% This will only work for three bodies so far

[Q1,v1] = SyMBA(Q,v,m,levc,rlev,hlev,steps);
[levcp] = calclevorig(Q1,rlev);

if(levcp == levc) %step is the same after the fact
    levc = levcp;
    Q = Q1;
    v = v1;
    return
elseif(levcp > levc ) %step is now bigger; keep result
    levc = levcp;
    Q = Q1;
    v = v1;
    return
elseif(levcp < levc) %step is smaller; and means always smaller
    [Q1,v1] = SyMBA(Q,v,m,levcp,rlev,hlev,steps);
    levc = levcp;
    Q = Q1;
    v = v1;
    return

end




end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v,levc] = SyMBAsmooth(Q,v,m,rlev,hlev,hsub,levmax)
n = size(m,2);
rlevs = size(rlev,2);
p0 = calcm(v,m);
[levc] = calclev(Q,rlev);
% levmax = levc+1; %for safety, go two levels further

[Q] = mapSun(Q,v,hlev(1)/2,m);
for j=2:n 
    for k=j+1:n
        [v] = interpairsmooth(Q,v,hlev(1)/2,m,rlev,1,j,k);
        [Q,v] = driftop(Q,v,m,2,rlev,hlev,hsub,levmax,j,k);
        [v] = interpairsmooth(Q,v,hlev(1)/2,m,rlev,1,j,k);
    end
end
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = adjustSun(Q,v,m,p0); 







end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v,levc] = SyMBAsmoothfull(Q,v,m,rlev,hlev,hsub,rhill0)
n = size(m,2);
rlevs = size(rlev,2);
p0 = calcm(v,m);
[levc,junk] = calclevhillp(Q,v,m,rlev,rhill0);
% [levc,junk] = calclevhill(Q,v,m,rlev);
levmax = max(max(levc));
% levmax
levmax = 11;
kepdone = zeros(levmax,n);

[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = interpairsmoothfull(Q,v,hlev(1)/2,m,rlev,1);
[Q,v,kepdone] = driftopfull(Q,v,m,2,rlev,hlev,hsub,levmax,levc,kepdone);
[v] = interpairsmoothfull(Q,v,hlev(1)/2,m,rlev,1);
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = adjustSun(Q,v,m,p0); 







end





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v,levc,f1,vfunc,finc] = SyMBAdiscrete(Q,v,m,rlev,hlev,hsub,levc)
f1 = 0;
finc=0;
n = size(m,2);
[Q1,v1] = SyMBAdiscretemap(Q,v,m,hlev,hsub,levc);
% [levcp,rhill] = calclevhillp(Q1,v1,m,rlev,rhill0);
[levcp,vfunc] = calclevv(Q1,v1,m,rlev,hlev(1));
diff = levcp - levc;

levcpp = levc;
flg = 0;
for i=2:n
    for j=i+1:n
        if(diff(i,j) < 0) %otherwise, same as original levels
            levcpp(i,j) = levcp(i,j);
            levcpp(j,i) = levcp(j,i);
            flg = 1;
        end
    end
end
% flg = 0;
if(flg == 1) %if any steps should be reduced, then redo
    [Q1,v1] = SyMBAdiscretemap(Q,v,m,hlev,hsub,levcpp);
%     [levcp1,rhill] = calclevhillp(Q1,v1,m,rlev,rhill0);
    [levcp1,vfunc] = calclevv(Q1,v1,m,rlev,hlev(1));
    for i=2:n
        for j=i+1:n
            if(levcp1(i,j) ~= levcp(i,j))
                display('inconsistent error')
                finc = 1;
%                 levcp1
%                 levcp
%                 levc
%                 blergh
            end
        end
    end
    levc = levcp; %also increasing steps
    Q = Q1;
    v = v1;
    f1 = 1;
else
    Q = Q1;
    v = v1;
    levc = levcp;
end






end

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [Q,v,levc1,levh1] = globalstep_h(Q,v,m,hsub,levc,collect_levels,rlev,hlev,frac,dE,levh,hsub,levv)
% loop = 0;
% levc0 = levc;
% %sum of steps will still be h0
% [Q1,v1,levels] = SyMBAdiscretemap_m_h(Q,v,m,hsub,levc,1,rlev,hlev,levh); %levels are max values
% [repeatl,check_againl,repeath,check_againh] = repcheck_h(levels,levc,levh);
% loop = 0;
% repeatl0 = repeatl;
% repeath0 = repeath;
% if(repeatl == 1 || repeath == 1)  %%make sure one doesn't cause other to change...
%     while(repeatl == 1 || repeath == 1)
%         loop = loop+1;
%         if(repeatl == 1)
%             levc = max(levels,[],3); %matrix values
%         end
%         if(repeath == 1)
%             levh = min(find(levels > 0),[],"all"); %just one value, minimum which is not 0
%         end
%         [Q1,v1,levelstmp] = SyMBAdiscretemap_m_h(Q,v,m,hsub,levc,1,rlev,hlev,levh);
%         if(check_againl == 1 || check_againh == 1)
%             levels = levelstmp;
%         end
%         [repeatl,check_againl,repeath,check_againh] = repcheck(levels,levc,levh); %hopefully just once
%     end
%     if(loop > 1)
%         display('loop')
%         loop
% %         blergh
%     end
%     if(repeatl0 == 1)
%         levc1 = findlastlev(levelstmp);
%     else
%         levc1
% 
% 
%     end
%     if(repeath0 == 1)
%         levh1 = levhn;
%         [levv] = updatelevv1(levv,levhn,hsub); %timestep may have gone lower
%     end
%     levh1 = levhn;
%     Q = Q1;
%     v = v1;
% else
%     [levv,levh1] = updatelevv2(levv,levels,hsub);
%     [levc1] = findlastlev(levels);
%     Q = Q1;
%     v = v1;
% end
% 
% 
% 
% 
% end
% 
% 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v,levc1,loop] = globalstep(Q,v,m,hsub,levc,collect_levels,rlev,hlev)
loop = 0;
levc0 = levc;
%sum of steps will still be h0
[Q1,v1,levels] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev); %levels are max values
[repeat,check_again] = repcheck(levels,levc);
loop = 0;
if(repeat == 1)
    while(repeat == 1)
        loop = loop+1;
%         levels
%         levc
        levc = max(levels,[],3);
%         levc
        [Q1,v1,levelstmp] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev);
        if(check_again == 1)
            levels = levelstmp;
        end
        [repeat,check_again] = repcheck(levels,levc); %hopefully just once
    end
    if(loop > 1)
        display('loop;LAKSDFLA;KJDF;LAKJSD;LAKSJDF')
        loop
%         blergh
    end
    levc1 = findlastlev(levelstmp);
    Q = Q1;
    v = v1;
%     loop
else
    [levc1] = findlastlev(levels);
    Q = Q1;
    v = v1;
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v,levc1,loop] = globalstepirr(Q,v,m,hsub,levc,collect_levels,rlev,hlev)
loop = 0;
levc0 = levc;
%sum of steps will still be h0
[Q1,v1,levels] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev); %levels are max values
[repeat,check_again] = repcheck(levels,levc);
loop = 0;
% if(repeat == 1)
%     while(repeat == 1)
%         loop = loop+1;
% %         levels
% %         levc
%         levc = max(levels,[],3);
% %         levc
%         [Q1,v1,levelstmp] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev);
%         if(check_again == 1)
%             levels = levelstmp;
%         end
%         [repeat,check_again] = repcheck(levels,levc); %hopefully just once
%     end
%     if(loop > 1)
%         display('loop;LAKSDFLA;KJDF;LAKJSD;LAKSJDF')
%         loop
% %         blergh
%     end
%     levc1 = findlastlev(levelstmp);
%     Q = Q1;
%     v = v1;
% %     loop
% else
    [levc1] = findlastlev(levels);
    Q = Q1;
    v = v1;
% end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v,levc1,loop] = globalstepp(Q,v,m,hsub,levc,collect_levels,rlev,hlev)
loop = 0;
levc0 = levc;
%sum of steps will still be h0
[Q1,v1,levels] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev); %levels are max values
[repeat,check_again] = repcheck(levels,levc);
loop = 0;
if(repeat == 1)
    display('here')
    loop = loop+1;
    levc
    levc = max(levels,[],3);
    levc
    [Q1,v1,levelstmp] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,1,rlev,hlev);
    levc1 = findlastlev(levelstmp);
    levc1
    Q = Q1;
    v = v1;
else
    [levc1] = findlastlev(levels);
    Q = Q1;
    v = v1;
end




end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v,levc,levv,levu] = adstep(Q,v,m,hsub,levc,rlev,hlev,levv)
%sum of steps will still be h0

levu = levc;
[Q1,v1] = SyMBAsimple(Q,v,m,hlev(levc));
[levcp] = calclevsimple(Q1,v1,m,rlev,hlev(1));
rd = 0;
if(levcp > levc)
    rd = 1;
    [Q1,v1] = SyMBAsimple(Q,v,m,hlev(levcp));
    levu = levcp;
    levv(levcp) = levv(levcp) + 1; %add step after attempt
    if(levv(levcp) == hsub)
        levv(levcp) = 0;
    end
    Q = Q1;
    v = v1;
    levc = levcp;
else %case when lower or equal level
    levv(levc) = levv(levc) + 1;
    if(levv(levc) == hsub)
        levv(levc) = 0;
    end
    if(levcp < levc) 
        levt = levc;
        while(levt > levcp)
            if(levv(levt) == 0)
                levt = levt - 1;
            else
                levt = levt;
                break;
            end
        end
        levcp = levt; %go as low as possible in level
    end %endif
    Q = Q1;
    v = v1;
    levc = levcp; %may be levt or levcp from first lines (same as initial)
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [levc] = findlastlev(levels)
n = size(levels,1);
sub = size(levels,3);
levc = zeros(n,n);
for i=1:n
    for j=i+1:n
        for k=1:sub
            if(levels(i,j,k)~= 0 )
                levc(i,j) = levels(i,j,k);
                levc(j,i) = levc(i,j);
            end
        end
    end
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [repeat,check_again] = updatelevv1(levelm,levc,levh)
n = size(levc,2);
repeat = 0;
check_again = 0;
for i=2:n
    for j=i+1:n
        if(max(levelm(i,j,:)) > levc(i,j))
            repeat = 1;
            if(max(levelm(i,j,:)) > levc(i,j) + 1)
                check_again = 1;
            end
        end
    end
end

%does global timestep need to get adapted?
levn = max(1,min(levelm,[],"all"));
if(levn > levh)
    repeat = 1;
    if(levn > levh + 1)
        check_again = 1;
    end
end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [repeatl,check_againl,repeath,check_againh] = repcheck_h(levelm,levc,levh)
n = size(levc,2);
repeatl = 0;
check_againl = 0;
repeath = 0;
check_againh = 0;
for i=2:n
    for j=i+1:n
        if(max(levelm(i,j,:)) > levc(i,j) )
            repeatl = 1;
            if(max(levelm(i,j,:)) > levc(i,j) + 1 )
                check_againl = 1;
            end
        end
    end
end

levelmp = min(find(levelm > 0),[],"all"); 
if(levelmp > levh)
    repeath = 1;
    if(levelmp > levh + 1)
        check_againh = 1;
    end
end

if(repeath == 1 && repeatl == 1)
    display('both repeat...')
    blergh
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [repeat,check_again] = repcheck(levelm,levc)
n = size(levc,2);
repeat = 0;
check_again = 0;
for i=2:n
    for j=i+1:n
        if(max(levelm(i,j,:)) > levc(i,j))
            repeat = 1;
            if(max(levelm(i,j,:)) > levc(i,j) + 1)
                check_again = 1;
            end
        end
    end
end




end






% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v] = SyMBAdiscretemap(Q,v,m,hlev,hsub,levc)

n = size(m,2);
p0 = calcm(v,m);
levmax = max(max(levc));
kepdone = zeros(levmax,n);
[E0,L0] = consqv(m,Q,v); 

[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = interfull(Q,v,hlev(1)/2,m,1,levc);
if(levmax > 1)
    [Q,v,kepdone] = driftopdiscrete(Q,v,m,2,hlev,hsub,levmax,levc,kepdone); %if levmax==1, this will return nothing
end
[Q,v,kepdone] = Kepfull(Q,v,hlev(1),m,1,levc,kepdone,levmax); %if levmax=1, do Kepler solver here

[v] = interfull(Q,v,hlev(1)/2,m,1,levc);
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = adjustSun(Q,v,m,p0); 


end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v,levsub] = SyMBAdiscretemap_m_h(Q,v,m,hsub,levc,collect_levels,rlev,hlev,levh)

n = size(m,2);
p0 = calcm(v,m);
levmax = max(max(levc));
kepdone = zeros(levmax,n);
levsub = zeros(n,n,1);
[E0,L0] = consqv(m,Q,v); 


[levsub] = calclevvsub(Q,v,m,rlev,hlev(1),1,levc,levsub); %beginning, end of step, and in between step check
[Q] = mapSun(Q,v,hlev(levh)/2,m);
[v] = interfull(Q,v,hlev(levh)/2,m,1,levc);
if(levmax > 1)
    [Q,v,kepdone,levsub] = driftopdiscrete_m(Q,v,m,2,hlev,hsub,levmax,levc,kepdone,levsub,rlev); %if levmax==1, this will return nothing
end
[Q,v,kepdone] = Kepfull(Q,v,hlev(levh),m,1,levc,kepdone,levmax); %if levmax=1, do Kepler solver here
[v] = interfull(Q,v,hlev(levh)/2,m,1,levc);
[levsub] = calclevvsub(Q,v,m,rlev,hlev(1),1,levc,levsub);
[Q] = mapSun(Q,v,hlev(levh)/2,m);
[v] = adjustSun(Q,v,m,p0); 


end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Q,v,levsub] = SyMBAdiscretemap_m(Q,v,m,hsub,levc,collect_levels,rlev,hlev)

n = size(m,2);
p0 = calcm(v,m);
levmax = max(max(levc));
kepdone = zeros(levmax,n);
levsub = zeros(n,n,1);
[E0,L0] = consqv(m,Q,v); 


[levsub] = calclevvsub(Q,v,m,rlev,hlev(1),1,levc,levsub); %beginning, end of step, and in between step check
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = interfull(Q,v,hlev(1)/2,m,1,levc);
if(levmax > 1)
    [Q,v,kepdone,levsub] = driftopdiscrete_m(Q,v,m,2,hlev,hsub,levmax,levc,kepdone,levsub,rlev); %if levmax==1, this will return nothing
end
[Q,v,kepdone] = Kepfull(Q,v,hlev(1),m,1,levc,kepdone,levmax); %if levmax=1, do Kepler solver here
[v] = interfull(Q,v,hlev(1)/2,m,1,levc);
[levsub] = calclevvsub(Q,v,m,rlev,hlev(1),1,levc,levsub);
[Q] = mapSun(Q,v,hlev(1)/2,m);
[v] = adjustSun(Q,v,m,p0); 


end





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levc] = calclev(Q,rlev)

si = size(rlev,2);
n = size(Q,2);
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
for j=2:n
    for k=j+1:n
        r = sqrt(sum((Q(:,j) - Q(:,k)).^2));
        if(r > rlev(1)) %if greater than first curoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            return
        end
        for i=1:si
            if((r <= rlev(i)) && (r > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end
    end
end



end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levc,rhill] = calclevhill(Q,v,m,rlev)
global GNEWT
si = size(rlev,2);
n = size(Q,2);
rhill = zeros(n,n);
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
% rhmax = 0.0929;
for j=2:n
    flgj = 0;
    muj = GNEWT*(m(1) + m(j));
    aj = 1/(2/norm(Q(:,j)) - sum((v(:,j)- v(:,1)).^2)/muj);
    if(aj < 0)
%         display('ejection')
        flgj = 1;
    end
    for k=j+1:n
        flgk = 0;
        muk = GNEWT*(m(1) + m(k));
        ak = 1/(2/norm(Q(:,k)) - sum((v(:,k)- v(:,1)).^2)/muk);
        if(ak < 0)
%             display('ejection')
            flgk = 1;
        end
        if( m(j) > 0 && m(k) > 0)
            rmut = ((m(j) + m(k))/(3*m(1)))^(1/3)*1/2*(aj+ak);
        elseif( m(j) > 0)
            rmut = (m(j)/(3*m(1)))^(1/3)*aj;
        elseif(m(k) > 0)
            rmut = (m(k)/(3*m(1)))^(1/3)*ak;
        else
            display('two zero mass colliding, bad')
            blergh
        end
        r = sqrt(sum((Q(:,j) - Q(:,k)).^2));
%         if(rmut > rhmax || rmut < 0)
%             rmut = rhmax;
%         end
        rhill(j,k) = rmut;
        rhill(k,j) = rmut;
%         rmut

        if(r/rmut > rlev(1)) %if greater than first curoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            continue
        end
        for i=1:si
            if((r/rmut <= rlev(i)) && (r/rmut > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end
    end
end
% rhill



end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levc,vfuncarr] = calclevv(Q,v,m,rlev,hglob)
global GNEWT
% rhill = 0.355*(50)^(1/3);
rhill = 0.4120*(50)^(1/3); %Saturn
% rlevel is really a timestep level
si = size(rlev,2);
n = size(Q,2);
vfunc = zeros(n,n);
vfuncarr = 0;
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
for j=2:n
    for k=j+1:n
        Qjk = norm(Q(:,j)-Q(:,k));
        sig = dot(v(:,j),v(:,k));
        Vjk2 = sum((v(:,j)-v(:,k)).^2);
        t1 = Qjk/sqrt(Vjk2);
        t2 = sqrt((Qjk*Qjk*Qjk)/(GNEWT*(m(j)+m(k))));
        if(sig < 0)
            Vfunc = min(t1,t2);
        else
            Vfunc = t2;
        end
        Vfunc = t2;
        ratio = Vfunc/hglob;
        vfuncarr(j,k) = Vfunc/hglob;
        vfuncarr(k,j) = Vfunc/hglob;
        ratio = Qjk/rhill;

        if(ratio > rlev(1)) %if greater than first cutoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            continue
        end
        for i=1:si
            if((ratio <= rlev(i)) && (ratio > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end
%         levc(j,k) = 1;
%         levc(k,j) = 1;
    end
end



end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levreturn] = calclevsimple(Q,v,m,rlev,hglob)
global GNEWT
rhill = 0.355;
rhill = 0.4120*(50)^(1/3); %Saturn
% rlevel is really a timestep level
si = size(rlev,2);
n = size(Q,2);
vfunc = zeros(n,n);
vfuncarr = 0;
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
for j=2:n
    for k=j+1:n
        Qjk = norm(Q(:,j)-Q(:,k));
        sig = dot(v(:,j),v(:,k));
        Vjk2 = sum((v(:,j)-v(:,k)).^2);
        t1 = Qjk/sqrt(Vjk2);
        t2 = sqrt((Qjk*Qjk*Qjk)/(GNEWT*(m(j)+m(k))));
        if(sig < 0)
            Vfunc = min(t1,t2);
        else
            Vfunc = t2;
        end
        Vfunc = t2;
        ratio = Vfunc/hglob;
        vfuncarr(j,k) = Vfunc/hglob;
        vfuncarr(k,j) = Vfunc/hglob;
        ratio = Qjk/rhill;

        if(ratio > rlev(1)) %if greater than first cutoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            continue
        end
        for i=1:si
            if((ratio <= rlev(i)) && (ratio > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end
%         levc(j,k) = 1;
%         levc(k,j) = 1;
    end
end

levreturn = max(levc,[],'all');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [levsub] = calclevvsub(Q,v,m,rlev,hglob,lev,levc,levsub)
global GNEWT
% rlevel is really a timestep level
si = size(rlev,2);
n = size(Q,2);
vfunc = zeros(n,n);
% rhill = 0.355*(50)^(1/3);
rhill = 0.4120*(50)^(1/3); %Saturn
%need to initialize levsub outside function to level 1

lasta = size(levsub,3); %adds a substep time, first important one is element 2
for j=2:n
    for k=j+1:n
            if(levc(j,k) == lev)
                Qjk = norm(Q(:,j)-Q(:,k));
                sig = dot(v(:,j),v(:,k));
                Vjk2 = sum((v(:,j)-v(:,k)).^2);
                t1 = Qjk/sqrt(Vjk2);
                t2 = sqrt((Qjk*Qjk*Qjk)/(GNEWT*(m(j)+m(k))));
                if(sig < 0)
                    Vfunc = min(t1,t2);
                else
                    Vfunc = t2;
                end
                Vfunc = t2;
                ratio = Vfunc/hglob;
                ratio = Qjk/rhill;
        
                if(ratio > rlev(1)) %if greater than first cutoff, then done
                    levsub(j,k,lasta+1) = 1;
                    levsub(k,j,lasta+1) = 1;
                    continue
                end
                for i=1:si
                    if((ratio <= rlev(i)) && (ratio > rlev(i+1)) )
                        levsub(j,k,lasta+1) = i;
                        levsub(k,j,lasta+1) = i;
                        break
                    end
                end
            end
    end %end particle loop
end



end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levc,vfuncarr] = checklevv(Q,v,m,rlev,hglob)
global GNEWT
% rlevel is really a timestep level
si = size(rlev,2);
n = size(Q,2);
vfunc = zeros(n,n);
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
for j=2:n
    for k=j+1:n
        Qjk = norm(Q(:,j)-Q(:,k));
        sig = dot(v(:,j),v(:,k));
        Vjk2 = sum((v(:,j)-v(:,k)).^2);
        t1 = Qjk/sqrt(Vjk2);
        t2 = sqrt((Qjk*Qjk*Qjk)/(GNEWT*(m(j)+m(k))));
        if(sig < 0)
            Vfunc = min(t1,t2);
        else
            Vfunc = t2;
        end
        Vfunc = t2;
        ratio = Vfunc/hglob;
        vfuncarr(j,k) = Vfunc/hglob;
        vfuncarr(k,j) = Vfunc/hglob;

        if(ratio > rlev(1)) %if greater than first cutoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            continue
        end
        for i=1:si
            if((ratio <= rlev(i)) && (ratio > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end

    end
end



end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





function [levc,rhill] = calclevhillp(Q,v,m,rlev,rhill0)
global GNEWT
si = size(rlev,2);
n = size(Q,2);
rhill = zeros(n,n);
levc = zeros(n,n); %0 column, row, unused, diagonal is 0
rhmax = 5.0;
for j=2:n
    flgj = 0;
    muj = GNEWT*(m(1) + m(j));
    aj = 1/(2/norm(Q(:,j)) - sum((v(:,j)- v(:,1)).^2)/muj);
    if(aj < 0)
%         display('ejection')
        flgj = 1;
    end
    for k=j+1:n
        flgk = 0;
        muk = GNEWT*(m(1) + m(k));
        ak = 1/(2/norm(Q(:,k)) - sum((v(:,k)- v(:,1)).^2)/muk);
        if(ak < 0)
%             display('ejection')
            flgk = 1;
        end
        if( m(j) > 0 && m(k) > 0)
            rmut = ((m(j) + m(k))/(3*m(1)))^(1/3)*1/2*(aj+ak);
        elseif( m(j) > 0)
            rmut = (m(j)/(3*m(1)))^(1/3)*aj;
        elseif(m(k) > 0)
            rmut = (m(k)/(3*m(1)))^(1/3)*ak;
        else
            display('two zero mass colliding, bad')
            blergh
        end
        r = sqrt(sum((Q(:,j) - Q(:,k)).^2));
        if(rmut > rhmax || rmut < 0)
            rmut = rhmax;
        end
        rhill(j,k) = rhill0(j,k);
        rhill(k,j) = rhill0(j,k);

        if(r/rmut > rlev(1)) %if greater than first curoff, then done
            levc(j,k) = 1;
            levc(k,j) = 1;
            continue
        end
        for i=1:si
            if((r/rmut <= rlev(i)) && (r/rmut > rlev(i+1)) )
                levc(j,k) = i;
                levc(k,j) = i;
                break
            end
        end
    end
end



end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [levc] = calclevorig(Q,rlev)

si = size(rlev,2);
n = size(Q,2);
for j=2:n
    for k=j+1:n
        r = sqrt(sum((Q(:,j) - Q(:,k)).^2));
        if(r > rlev(1)) %if greater than first curoff, then done
            levc = 1;
            return
        end
        for i=1:si
            if((r <= rlev(i)) && (r > rlev(i+1)) )
                levc = i;
                break
            end
        end
    end
end



end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [rlev,hlev,steps,hsub] = init(rad,h)
levtot = 30;
rmax = 3*rad; %Outside this value, maximum time step
rmin = 0;
delr = rmax-rmin;
rsub = 2;
hsub = 3;

rlev(1) = rmax;
hlev(1) = h;
steps(1) = 1;
for i=2:levtot
    rlev(i) = delr/(rsub^(i-1));
    hlev(i) = h/(hsub^(i-1)); 
    steps(i) = steps(i-1)*hsub;
end





end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [rlev,hlev,steps,hsub] = inithill(h)
levtot = 30;
rmax = 3; %Units of hill radius, outside this value is maximum time step
rmin = 0;
delr = rmax-rmin;
rsub = 2;
hsub = 3;

rlev(1) = rmax;
hlev(1) = h;
steps(1) = 1;
for i=2:levtot
    rlev(i) = delr/(rsub^(i-1));
    hlev(i) = h/(hsub^(i-1)); 
    steps(i) = steps(i-1)*hsub;
end





end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [rlev,hlev,steps,hsub] = initv(h)
levtot = 30;
% rmax = 30; %Units of TIME STEP FUNCTION/h_{global}, free fall and close encounter times; if small, need reduce steps, to keep ratio constant...
% rmax = 5;
rmax = 3;
rmin = 0;
delr = rmax-rmin;
rsub = 2; 
hsub = 4;

rlev(1) = rmax;
hlev(1) = h;
steps(1) = 1;
for i=2:levtot
    rlev(i) = delr/(rsub^(i-1)); %ratio levels of funcv/h
    hlev(i) = h/(hsub^(i-1));  %keeps track of timestep at each level
    steps(i) = steps(i-1)*hsub;
end





end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [hp] = initializeplot(q,p)
    lim = 6;
    figure()
    set(gca,'FontSize',25)
    hp=plot(q,p,'o');
    hold on
    xmin = -lim;
    xmax = lim;
    ymin = -lim;
    ymax = lim;
    n = 100;
    axis([xmin xmax ymin ymax])
    axis square
    grid off
    set(gca,'FontSize',25);
end











function runtime_output(h,q,p)
set(h,'XData',q,'YData',p);
set(h,'MarkerSize',6,'MarkerFaceColor','k');
drawnow

end






















function [x,v] = in2bod()
e = 0.5;
mu = 1;
M = 1;
x(1) = (1 + e);
x(2) = 0;
v(1) = 0;
v(2) = sqrt(abs((1-e)/(1+e)));

end


function [m,x0,v0,pair] = in2bodkep(nbin)
m0 = 1/2;
a = 1; 
e = 0.5;
m = [m0 m0];
n = 2;
mu = (m(1)*m(2))/(m(1)+m(2));
mt = m(1)+m(2);
xp = a*(1+e);
vp = sqrt(mt/a*(1-e)/(1+e));
x0 = [-xp/2 0 0; xp/2 0 0];
x0 = x0';
v0 = [0 -vp/2 0; 0 vp/2 0];
v0 = v0';

a = 1e-2;
e = .9;
[x0,v0,m] = binadd(x0,v0,m,nbin,a,e);
n = size(m,2);
for i=1:n
    for j=1:n
%         Kepler solver group
        pair(i,j) = 0;
        pair(j,i) = 0;
    end
end

end






function [m,xout,vout] = insimple()
global YEAR GNEWT
%Below data from Hairer pg. 13-14.  Convert time to years.
%Distance is in AU, mass in solar masses. 

m = [1.00000597682 .000954786104043 .000285583733151 ...
    .0000437273164546 .0000517759138449];

n = size(m,2);
xout = [0 0 0; ...
    -3.5023653 -3.8169847 -1.5507963; ...
    9.0755314 -3.0458353 -1.6483708; ...
    8.3101420 -16.2901086 -7.2521278; ...
    11.4707666 -25.7294829 -10.8169456];
xout = xout';
vout = [0 0 0; ...
    .00565429 -0.00412490 -0.00190589; 
    .00168318 .00483525 .00192462; ...
    .00354178 .00137102 .00055029; ...
    .00288930 .00114527 .00039677];
vout = vout*YEAR;
vout = vout';
% Only go up to Saturn
n = 3;
xout = xout(:,1:n);
vout = vout(:,1:n);
m = m(1:n);





vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vout(:,j);
    xcm(:) = xcm(:)+m(j)*xout(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary; barycentric velocities
for j=1:n
    vout(:,j) = vout(:,j)-vcm(:); 
    xout(:,j) = xout(:,j)-xcm(:);
end


% E0 = -0.004294287207322 in solar mass, yr, AU

end





    




function [m,xout,vout] = insolardunc()
global YEAR GNEWT
%Below data from Hairer pg. 13-14.  Convert time to years.
%Distance is in AU, mass in solar masses. 
fac = 50;
m = fac*[1.00000597682 .000954786104043 .000285583733151 ...
    .0000437273164546 .0000517759138449];
m(1) = 1.00000597682;
m = m(1:5);

n = size(m,2);

xout = [0 0 0; ...
    -3.5023653 -3.8169847 -1.5507963; ...
    9.0755314 -3.0458353 -1.6483708; ...
    8.3101420 -16.2901086 -7.2521278; ...
    11.4707666 -25.7294829 -10.8169456];
xout = xout';
vout = [0 0 0; ...
    .00565429 -0.00412490 -0.00190589; 
    .00168318 .00483525 .00192462; ...
    .00354178 .00137102 .00055029; ...
    .00288930 .00114527 .00039677];
vout = vout*YEAR;
vout = vout';
xout = xout(:,1:n);
vout = vout(:,1:n);


vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vout(:,j);
    xcm(:) = xcm(:)+m(j)*xout(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary
for j=1:n
    vout(:,j) = vout(:,j)-vcm(:); 
    xout(:,j) = xout(:,j)-xcm(:);
end



% E0 = -0.004294287207322 in solar mass, yr, AU

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [m,xout,vout,pair] = insolarpluto()
global YEAR GNEWT
%Below data from Hairer pg. 13-14.  Convert time to years.
%Distance is in AU, mass in solar masses.  
fac = 1;
m = fac*[1.00000597682 .000954786104043 .000285583733151 ...
    .0000437273164546 .0000517759138449 6.58086572e-9];
m(1) = 1.00000;
n = size(m,2);
for i=1:n
    for j=1:n
%         Kepler solver group
        if(i == 1 || j == 1) 
            pair(i,j) = 0;
            pair(j,i) = 0;
        else
            pair(i,j) = 1;
            pair(j,i) = 1;
        end
    end
end
xout = [-2.079997415328555E-04  7.127853194812450E-03 -1.352450694676177E-05; ...
    -3.502576700516146E+00 -4.111754741095586E+00  9.546978009906396E-02;...
    9.075323061767737E+00 -3.443060862268533E+00 -3.008002403885198E-01;...
    8.309900066449559E+00 -1.782348877489204E+01 -1.738826162402036E-01;...
     1.147049510166812E+01 -2.790203169301273E+01  3.102324955757055E-01;...
     -1.553841709421204E+01 -2.440295115792555E+01  7.105854443660053E+00];
xout = xout';
vout = [-6.227982601533108E-06  2.641634501527718E-06  1.564697381040213E-07; ...
    5.647185656190083E-03 -4.540768041260330E-03 -1.077099720398784E-04; ... 
    1.677252499111402E-03  5.205044577942047E-03 -1.577215030049337E-04; ...
    3.535508197097127E-03  1.479452678720917E-03 -4.019422185567764E-05; ...
    2.882592399188369E-03  1.211095412047072E-03 -9.118527716949448E-05; ...
     2.754640676017983E-03 -2.105690992946069E-03 -5.607958889969929E-04];
 vout = vout';
vout = vout*YEAR;

xout = xout(:,1:n);
vout = vout(:,1:n);


vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vout(:,j);
    xcm(:) = xcm(:)+m(j)*xout(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary
for j=1:n
    vout(:,j) = vout(:,j)-vcm(:); 
    xout(:,j) = xout(:,j)-xcm(:);
end

% Add CoM velocity
vcm = 100*vcm;
vcm = 0*vcm;

for j=1:n
    vout(:,j) = vout(:,j) + vcm(:); 
end


end




function [m,xout,vout,pair] = inthreebod()
global GNEWT
%  au, yr, M_{\odot} units
n = 3;
for i=1:n
    for j=1:n
%         Kepler solver group
        if(i == 1 || j == 1) 
            pair(i,j) = 0;
            pair(j,i) = 0;
        else
            pair(i,j) = 1;
            pair(j,i) = 1;
        end
    end
end

m0 = 1e-4;
m = [1 m0 m0];

a1 = 5.2;
e1 = 0.9;
% e1 = 0;
mu1 = (m(1)*m(2))/(m(1)+m(2));
mt1 = m(1)+m(2);
x1 = a1*(1+e1);
v1 = sqrt(GNEWT*mt1/a1*(1-e1)/(1+e1));
fa1 = m(2)/mt1;
fb1 = m(1)/mt1;

a2 = 10.0;
e2 = 0;
e2 = 0.9;
mu2 = (m(1)*m(3))/(m(1)+m(3));
mt2 = m(1)+m(3);
x2 = a2*(1+e2);
v2 = sqrt(GNEWT*mt2/a2*(1-e2)/(1+e2));
fa2 = m(3)/mt2;
fb2 = m(1)/mt2;

xout = [0 0 0; -fb1*x1 0 0; -fb2*x2 0 0];
xout = xout';
vout = [0 0 0; 0 fb1*v1 0; 0 fb2*v2 0];
vout = vout';



vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vout(:,j);
    xcm(:) = xcm(:)+m(j)*xout(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary
for j=1:n
    vout(:,j) = vout(:,j)-vcm(:); 
    xout(:,j) = xout(:,j)-xcm(:);
end

% Add CoM velocity
vcm = 0*vcm;

for j=1:n
    vout(:,j) = vout(:,j) + vcm(:); 
end



end




function [Q,P] = convertcart(m,q,v)
global CASE
n = size(m,2);
M = sum(m);
Q = zeros(3,n);
P = zeros(3,n);
for i=1:n
    p(:,i) = m(i)*v(:,i);
end

if(CASE == 1 || CASE == 2 || CASE == 3)
    if(CASE == 1 || CASE == 2)
        a = m/M;
    else
       a(1) = 1;
       a(2:n) = 0;
    end
    for i=1:n
        Q(:,1) = Q(:,1) + a(i)*q(:,i); 
        P(:,1) = P(:,1) + p(:,i);
    end
    for i=2:n
       P(:,i) = p(:,i) - a(i)*P(:,1); 
       Q(:,i) = q(:,i) - q(:,1); 
    end
elseif(CASE == 4) %WHJ case
    g = q(:,1);
    G = p(:,1);
    M = m(1);
    for i=2:n
       Q(:,i) = q(:,i) - g;
       Mp = M + m(i);
       g = 1/Mp*(m(i)*q(:,i) + M*g);
       P(:,i) = M/Mp* p(:,i) - m(i)/Mp*G;
       G = p(:,i) + G;
       M = Mp;
    end
    Q(:,1) = g;
    P(:,1) = G;   


elseif(CASE == 5 || CASE == 6 || CASE == 7 || CASE == 8 || CASE == 9)
%     Simply Cartesian coordinates
    Q = q;
    P = p;
end
    
end



function [Q,V] = convertcartv(m,q,v)
global CASE
n = size(m,2);
M = sum(m);
Q = zeros(3,n);
V = zeros(3,n);

if(CASE == 1 || CASE == 2 || CASE == 3)
    if(CASE == 1 || CASE == 2)
        a = m/M;
    else
       a(1) = 1;
       a(2:n) = 0;
    end
    for i=1:n
        Q(:,1) = Q(:,1) + a(i)*q(:,i); 
        V(:,1) = V(:,1) + m(i)/M*v(:,i);
    end
    for i=2:n
       V(:,i) = v(:,i) - V(:,1); 
       Q(:,i) = q(:,i) - q(:,1); 
    end
    V
    v
    blergh
elseif(CASE == 4) %WHJ case
    g = q(:,1);
    G = p(:,1);
    M = m(1);
    for i=2:n
       Q(:,i) = q(:,i) - g;
       Mp = M + m(i);
       g = 1/Mp*(m(i)*q(:,i) + M*g);
       V(:,i) = v(:,i) - 1/M*G;
       G = p(:,i) + G;
       M = Mp;
    end
    Q(:,1) = g;
    V(:,1) = G/M; %center of mass vel


elseif(CASE == 5 || CASE == 6 || CASE == 7 || CASE == 8 || CASE == 9)
%     Simply Cartesian coordinates
    Q = q;
    P = p;
end
    
end





function [Q,P] = convertcartp(m,q,p)
global CASE
n = size(m,2);
M = sum(m);
Q = zeros(3,n);
P = zeros(3,n);

if(CASE == 1 || CASE == 2 || CASE == 3)
    if(CASE == 1 || CASE == 2)
        a = m/M;
    else
       a(1) = 1;
       a(2:n) = 0;
    end
    for i=1:n
        Q(:,1) = Q(:,1) + a(i)*q(:,i); 
        P(:,1) = P(:,1) + p(:,i);
    end
    for i=2:n
       P(:,i) = p(:,i) - a(i)*P(:,1); 
       Q(:,i) = q(:,i) - q(:,1); 
    end
elseif(CASE == 4) %WHJ case
    g = q(:,1);
    G = p(:,1);
    M = m(1);
    for i=2:n
       Q(:,i) = q(:,i) - g;
       Mp = M + m(i);
       g = 1/Mp*(m(i)*q(:,i) + M*g);
       P(:,i) = M/Mp* p(:,i) - m(i)/Mp*G;
       G = p(:,i) + G;
       M = Mp;
    end
    Q(:,1) = g;
    P(:,1) = G;   


elseif(CASE == 5 || CASE == 6 || CASE == 7 || CASE == 8 || CASE == 9)
%     Simply Cartesian coordinates
    Q = q;
    P = p;
end
    
end


function [q,p] = convert2cart(m,Q,P)
global CASE
n = size(m,2);
M = sum(m);
p = zeros(3,n);
q = zeros(3,n);

if(CASE == 1 || CASE == 2 || CASE == 3)
    if(CASE == 1 || CASE == 2)
        a = m/M;
    else
       a(1) = 1;
       a(2:n) = 0;
    end
    s1 = zeros(3,1);
    s2 = zeros(3,1);
    for i=2:n
        s1 = s1 + a(i)*Q(:,i);
        s2 = s2 + P(:,i);
    end
    p(:,1) = a(1)*P(:,1) - s2;
    q(:,1) = Q(:,1) - s1;
    for i=2:n
       p(:,i) =  a(i)*P(:,1) + P(:,i); 
       q(:,i) = Q(:,1) + Q(:,i) - s1; 
    end
elseif(CASE == 4) %WHJ case
    g = Q(:,1);
    G = P(:,1);
    M = sum(m);
    for i=n:-1:2
       Mp = M - m(i);
       q(:,i) = Mp/M*Q(:,i) + g;
       g = -m(i)/M*Q(:,i) + g;
       p(:,i) = P(:,i) + m(i)/M*G;
       G = -P(:,i) + Mp/M*G;
       M = Mp;
    end
    q(:,1) = g;
    p(:,1) = G;   
elseif(CASE == 5 || CASE == 6 || CASE == 7 || CASE == 8 || CASE == 9)
%     Simply Cartesian coordinates
    q = Q;
    p = P;
end
    
end





function [q,v] = convert2cartv(m,Q,V)
global CASE
n = size(m,2);
M = sum(m);
v = zeros(3,n);
q = zeros(3,n);


if(CASE == 4) %WHJ case
    g = Q(:,1);
    M = sum(m);
    G = V(:,1)*M; %total mom
    for i=n:-1:2
       Mp = M - m(i);
       q(:,i) = Mp/M*Q(:,i) + g;
       g = -m(i)/M*Q(:,i) + g;
       v(:,i) = Mp/M*V(:,i) + 1/M*G;
       G = -V(:,i)*m(i)*Mp/M + Mp/M*G;
       M = Mp;
    end
    q(:,1) = g;
    v(:,1) = G/m(1);   
end
    
end




function [q] = pos2cart(Q,m)
n = size(m,2);
q = zeros(3,n);
g = Q(:,1);

M = sum(m);
for i=n:-1:2
   Mp = M - m(i);
   q(:,i) = Mp/M*Q(:,i) + g;
   g = -m(i)/M*Q(:,i) + g;
   M = Mp;
end
q(:,1) = g; 



end

function [q,p] = cart2jac(Q,P,m)
n = size(m,2);
g = Q(:,1);
G = P(:,1);
M(1) = m(1);
for i=2:n
   M(i) = M(i-1)+m(i); 
end
for i=n:-1:2
   q(:,i) = M(i-1)/M(i)*Q(:,i) + g;
   p(:,i) = P(:,i) + m(i)/M(i)*G;
   g = -m(i)/M(i)*Q(:,i) + g;
   G = -P(:,i) + M(i-1)/M(i)*G;
end
q(:,1) = g;
p(:,1) = G;

end

function [Pdot] = mom2jac(pdot,m)
n = size(m,2);
G(:,1) = pdot(:,1);
M(1) = m(1);
for i=2:n
   M(i) = M(i-1)+m(i); 
end
for i=2:n
   Pdot(:,i) = M(i-1)/M(i)*pdot(:,i) - m(i)/M(i)*G;
   G = pdot(:,i) + G;
end
Pdot(:,1) = G; %this should be 0, change in total mom


end

function [Vdot] = mom2jacv(vdot,m)
n = size(m,2);
G(:,1) = m(1)*vdot(:,1);
M(1) = m(1);
for i=2:n
   M(i) = M(i-1)+m(i); 
end
for i=2:n
   Vdot(:,i) = vdot(:,i) - 1/M(i-1)*G;
   G = m(i)*vdot(:,i) + G;
end
Vdot(:,1) = G/M(n); %CM velocity


end


function [m,x0,p0] = inplummer(n)
%Standard units, check Chapter 10 Vol. 11
%Art of Computational Science
% E = -1/4;
% rv = 1;
% rh ~ 4/5
% vrms ~.71
% => tcross ~ 1.4

m = ones(1,n)./n;
fidx = fopen('xplummer.txt','r');
fidv = fopen('vplummer.txt','r');
x0 = fscanf(fidx,'%g %g %g\n',[3 n]);
v0 = fscanf(fidv,'%g %g %g\n',[3 n]);
xcm = zeros(3,1);
vcm = zeros(3,1);
for i=1:n
     xcm(:) = xcm(:) + m(i)*x0(:,i);
     vcm(:) = vcm(:) + m(i)*v0(:,i);
end
for i=1:n
    x0(1:3,i) = x0(1:3,i) - xcm(1:3);
    v0(1:3,i) = v0(1:3,i) - vcm(1:3);
end
scalefac = 16/(3*pi); %for nbody units
x0 = x0/scalefac;
v0 = v0*sqrt(scalefac);
for i=1:n
   p0(1:3,i) = v0(1:3,i)*m(i); 
end

end


function [m,x0,p0] = inplummerbin(n,nbin)
%Standard units, check Chapter 10 Vol. 11
%Art of Computational Science
% E = -1/4;
% rv = 1;
% rh ~ 4/5
% vrms ~.71
% => tcross ~ 1.4

m = ones(1,n)./n;
fidx = fopen('xplummer.txt','r');
fidv = fopen('vplummer.txt','r');
x0 = fscanf(fidx,'%g %g %g\n',[3 n]);
v0 = fscanf(fidv,'%g %g %g\n',[3 n]);
xcm = zeros(3,1);
vcm = zeros(3,1);
for i=1:n
     xcm(:) = xcm(:) + m(i)*x0(:,i);
     vcm(:) = vcm(:) + m(i)*v0(:,i);
end
for i=1:n
    x0(1:3,i) = x0(1:3,i) - xcm(1:3);
    v0(1:3,i) = v0(1:3,i) - vcm(1:3);
end
scalefac = 16/(3*pi); %for nbody units
x0 = x0/scalefac;
v0 = v0*sqrt(scalefac);
p0 = vtop(v0,m);

a = 1e-2;
e = .9;
[x0,p0,m] = binadd(x0,p0,m,nbin,a,e);

end











function [m,x0,p0] = in3bodsuvakov(i)
m = [1 1 1];
x0 = [-1 0 0; 1 0 0; 0 0 0];
x0 = x0';
p1v = [0.30689 0.39295 0.18428 ...
    0.46444 0.43917 0.40592 0.38344 ...
    0.08330 0.350112 0.08058 0.55906 ...
    0.51394 0.28270 0.41682 0.41734];
p2v = [0.12551 0.09758 0.58719 0.39606...
    0.45297 0.23016 0.37736 0.12789 0.07934 ...
     0.58884 0.34919 0.30474 0.32721 ...
     0.33033 0.31310];
 Tv = [6.2356 7.0039 63.5345 14.8939 28.6703 13.8658...
    25.8406 10.4668 79.4759 21.2710 55.5018 ...
    17.3284 10.9626 55.7898 54.2076];
p1 = p1v(i);
p2 = p2v(i);
T = Tv(i);
v0 = [p1 p2 0; p1 p2 0; -2*p1 -2*p2 0];
v0 = v0';
n = size(m,2);
p0 = zeros(size(v0,1),size(v0,2));
for i=1:n
    p0(:,i) = v0(:,i)*m(i);
end
end





function [m,x0,v0,pair] = in3bodpyth()
% Need units where GNEWT is 1
m = [4 5 3];
x0 = [-2 -1 0; 1 -1 0; 1 3 0];
x0 = x0';
v0 = [0 0 0; 0 0 0; 0 0 0];
v0 = v0';
n = size(m,2);
for i=1:n
    for j=1:n
%         Kepler solver group
        pair(i,j) = 0;
        pair(j,i) = 0;
    end
end

end




function [m,xp,vp] = inwisd(ein)
global GNEWT
% step from 1e-4 to 1e0
% time span 3000 years

fac = 1;
n = 3;
% m = [1.0 fac*1/1047.355 fac*1/3498.5];
m = [1.0 fac*1/1047.355 0];

a1 = 5.2;
a2 = 9.58;
% a2 = a1;
e1 = 0.05;
e2 = ein;

x1 = a1*(1-e1);
x2 = a2*(1-e2);
v1 = sqrt(GNEWT*(m(1)+m(2))/a1*(1+e1)/(1-e1));
v2 = sqrt(GNEWT*(m(1)+m(3))/a2*(1+e2)/(1-e2));

% th = 0;
% R = [cos(th) 0 sin(th); 0 0 0; -sin(th) 0 cos(th)];
% c1 = [x2 0 0];
% c2 = [0 0 v2];
% c1r = c1*R;
% c2r = c2*R;
% 
% xp = [0 0 0; x1 0 0; c1r];
% xp = xp';
% vp = [0 0 0; 0 v1 0; c2r];
% vp = vp';


xp = [0 0 0; x1 0 0; x2 0 0];
xp = xp';
vp = [0 0 0; 0 v1 0; 0 0 v2];
vp = vp';


vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vp(:,j);
    xcm(:) = xcm(:)+m(j)*xp(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary; barycentric velocities
for j=1:n
    vp(:,j) = vp(:,j)-vcm(:); 
    xp(:,j) = xp(:,j)-xcm(:);
end





end







function [mt,xt,vt,om,order] = inwisdtest(ein)
global GNEWT
% step from 1e-4 to 1e0
% time span 3000 years

fac = 1;
n = 3;
m = [1.0 fac*1/1047.355 fac*1/3498.5];
% m = [1.0 fac*1/1047.355 0];
% m = [1.0 fac*1/100 0];

a1 = 5.2;
a2 = 9.58;
% a2 = a1;
e1 = 0;
e2 = ein;

x1 = a1*(1-e1);
x2 = a2*(1+e2);
v1 = sqrt(GNEWT*(m(1)+m(2))/a1*(1+e1)/(1-e1));
v2 = sqrt(GNEWT*(m(1)+m(3))/a2*(1-e2)/(1+e2)); %apocenter start

xp = [0 0 0; x1 0 0; x2 0 0];
xp = xp';
vp = [0 0 0; 0 v1 0; 0 0 v2];
vp = vp';


vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vp(:,j);
    xcm(:) = xcm(:)+m(j)*xp(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary; barycentric velocities
for j=1:n
    vp(:,j) = vp(:,j)-vcm(:); 
    xp(:,j) = xp(:,j)-xcm(:);
end

om = sqrt(GNEWT*(m(1)+m(2))/a1^3);
if(norm(xp(:,2)) < norm(xp(:,3)) )
    order = 1;
    xt = xp;
    vt = vp;
    mt = m;
else

%     order = 1;
%     xt = xp;
%     vt = vp;
%     mt = m;




    order = -1;
    xt = xp;
    vt = vp;
    mt = m;
    xt(:,2) = xp(:,3);
    xt(:,3) = xp(:,2);
    vt(:,2) = vp(:,3);
    vt(:,3) = vp(:,2);
    mt(2) = m(3);
    mt(3) = m(2);
end





end














function [m,xp,vp,om,order] = inwisdm(ein)
global GNEWT
% step from 1e-4 to 1e0
% time span 3000 years

fac = 1;
n = 3;
m = [1.0 fac*1/1047.355 fac*1/3498.5 fac*0.00004365];

a1 = 5.2;
a2 = 9.58;
a3 = 19.2;
e1 = 0.05;
e2 = 0.047;
e3 = ein;

x1 = a1*(1-e1);
x2 = a2*(1+e2);
x3 = a3*(1+e3);
v1 = sqrt(GNEWT*(m(1)+m(2))/a1*(1+e1)/(1-e1));
v2 = sqrt(GNEWT*(m(1)+m(3))/a2*(1-e2)/(1+e2)); %apocenter start
v3 = sqrt(GNEWT*(m(1)+m(4))/a2*(1-e3)/(1+e3)); %apocenter start

xp = [0 0 0; x1 0 0; x2 0 0; x3 0 0];
xp = xp';
vp = [0 0 0; 0 v1 0; 0 v2 0; 0 0 v3];
vp = vp';


vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*vp(:,j);
    xcm(:) = xcm(:)+m(j)*xp(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary; barycentric velocities
for j=1:n
    vp(:,j) = vp(:,j)-vcm(:); 
    xp(:,j) = xp(:,j)-xcm(:);
end

om = sqrt(GNEWT*(m(1)+m(2))/a1^3);
order = [1 2 3]; %number of planets size






end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [m,xout,vout] = in3boddunctwo()
global GNEWT
a0 = 0.0125;
a1 = 1.0;
e0 = 0.6;
e1 = 0;
a2 = 0.0130;
% a2 = 0.0125;
e2 = 0.2;
a1p = 3.0;
e1p = e1;
m0 = 1e-3;
m = [m0 m0 1 m0 m0];
mu0 = (m(1)*m(2))/(m(1)+m(2));
mt0 = m(1)+m(2);
mu1 = (m(3)*mt0)/(m(3)+mt0);
mt1 = m(3)+mt0;
x0 = a0*(1+e0);
x1 = a1*(1+e1);
x1p = a1p*(1+e1);
x2 = a2*(1+e2);
v0 = sqrt(GNEWT*mt0/a0*(1-e0)/(1+e0));
v1 = sqrt(GNEWT*mt1/a1*(1-e1)/(1+e1));
v1p = sqrt(GNEWT*mt1/a1p*(1-e1p)/(1+e1p));
v2 = sqrt(GNEWT*mt0/a2*(1-e2)/(1+e2));
fa0 = m(1)/mt0;
fb0 = m(2)/mt0;
fa1 = mt0/mt1;
fb1 = m(3)/mt1;

% switch bodies around
xout = [fa1*x1 0 0; -fb1*x1-fb0*x0 0 0; -fb1*x1+fa0*x0 0 0; -fb1*x1p-fb0*x2 0 0; -fb1*x1p+fa0*x2 0 0;];
xout = xout';
vout = [0 fa1*v1 0; 0 -fb1*v1-fb0*v0 0; 0 -fb1*v1+fa0*v0 0; 0 -fb1*v1p-fb0*v2 0; 0 -fb1*v1p+fa0*v2 0];
vout = vout';
m = [1 m0 m0 m0 m0];


[xout,vout] = adjustcm(xout,vout,m);

% P = 2*pi*sqrt(a0^3/(GNEWT*(2*m0)));
% P/365.25
% P = 2*pi*sqrt(a1p^3/(GNEWT*(m0+1)));
% P/365.25
% blergh
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [m,x0,v0] = in3bodfig8()
% The period is approximately T = 6.32591398
m = [1 1 1];
xa = .97000436;
xb = -.24308753;
va = -.93240737;
vb = -.86473146;
x0 = [xa xb 0; -xa -xb 0; 0 0 0];
x0 = x0';
v0 = [-va/2 -vb/2 0; -va/2 -vb/2 0; va vb 0];
v0 = v0';
end



function [consist] = consistency(x1,p1,x2,p2,h)
abp = abs((p2-p1)/h);
abp = max(max(abp));
abx = abs((x2-x1)/h);
abx = max(max(abx));
consist = max([abx abp]);

end


function [x,p] = hw(x,p,h,m)
n = size(m,2);
dim = size(p,1);
tau = h;
[v] = ptov(p,m);
dmr = zeros(dim,n);
dmv = zeros(dim,n);
for i=1:n
    for j=1:n
        if(j ~= i)
            mij = m(i) + m(j);
            mu = m(i)*m(j)/mij;
            rr0 = x(:,i)-x(:,j);
            vv0 = v(:,i)-v(:,j);
            r0 = rr0-vv0*tau/2;
            v0 = vv0;
            [delx,delv] = keplerm(mij,r0,v0,tau);
            r1 = r0+delx;
            v1 = v0 + delv;
            rr1 = r1 - v1*tau/2;
            vv1 = v1;
            dmr(:,i) = dmr(:,i) + mu*(rr1-rr0);
            dmv(:,i) = dmv(:,i) + mu*(vv1-vv0);
        end  
    end
end
for i=1:n
    x(:,i) = x(:,i) + dmr(:,i)/m(i);
    v(:,i) = v(:,i) + dmv(:,i)/m(i);
end
[p] = vtop(v,m);

end








function [x,p] = goncal(x,p,h,m)
n = size(m,2);
dim = size(p,1);
[x] = drift(x,p,h/2,m);
tau = h;
[v] = ptov(p,m);
dmr = zeros(dim,n);
dmv = zeros(dim,n);
for i=1:n
    for j=1:n
        if(j ~= i)
            mij = m(i) + m(j);
            mu = m(i)*m(j)/mij;
            rr0 = x(:,i)-x(:,j);
            vv0 = v(:,i)-v(:,j);
            r0 = rr0-vv0*tau/2;
            v0 = vv0;
            [delx,delv] = keplerm(mij,r0,v0,tau);
            r1 = r0+delx;
            v1 = v0 + delv;
            rr1 = r1 - v1*tau/2;
            vv1 = v1;
            dmr(:,i) = dmr(:,i) + mu*(rr1-rr0);
            dmv(:,i) = dmv(:,i) + mu*(vv1-vv0);
        end  
    end
end
for i=1:n
    x(:,i) = x(:,i) + dmr(:,i)/m(i);
    v(:,i) = v(:,i) + dmv(:,i)/m(i);
end
[p] = vtop(v,m);
[x] = drift(x,p,h/2,m);

end





function [x,p] = phi1(x,p,h,m)
n = size(m,2);
[x] = drift(x,p,h,m);
for i=1:(n-1)
    for j=(i+1):n
        [x] = driftij(x,p,i,j,-h/2,m);
        [x,p] = keplerij(m,x,p,i,j,h/2);
    end
end
end



function [x,p] = phi2(x,p,h,m)
n = size(m,2);
[x] = drift(x,p,h/2,m);
for i=1:(n-1)
    for j=(i+1):n
        [x] = driftij(x,p,i,j,-h/2,m);
        [x,p] = keplerij(m,x,p,i,j,h/2);
    end
end
for i=(n-1):-1:1
    for j=n:-1:(i+1)
        [x,p] = keplerij(m,x,p,i,j,h/2);
        [x] = driftij(x,p,i,j,-h/2,m);
    end
end
[x] = drift(x,p,h/2,m);

end









function [x,v] = hb15(x,v,h,m,pair)
n = size(m,2);

[x] = drift(x,v,h/2);
[v] = kickpair(x,v,h/2,m,pair);
for i=1:(n-1)
    for j=(i+1):n
        if(pair(i,j) == 0)
            [x] = driftij(x,v,i,j,-h/2,m);
            [x,v] = keplerij(m,x,v,i,j,h/2);
        end
    end
end
for i=(n-1):-1:1
    for j=n:-1:(i+1)
        if(pair(i,j) == 0)
            [x,v] = keplerij(m,x,v,i,j,h/2);
            [x] = driftij(x,v,i,j,-h/2,m);
        end
    end
end
[v] = kickpair(x,v,h/2,m,pair);
[x] = drift(x,v,h/2);

end




function [x,v] = whd(x,v,h,m)
% n = 1 is the sun here.  Using heliocentric positions and barycentric
% velocities.
global GNEWT
n = size(m,2);
gm = GNEWT*m(1);
[x] = driftwhd(x,v,h/2,m);
[v] = kickwhd(x,v,h/2,m);
for i=2:n
    xo = x(:,i);
    vo = v(:,i);
    [xn,vn] = kepler_stepxv(gm,xo,vo,h);
    x(:,i) = xn(:);
    v(:,i) = vn(:);
end
[v] = kickwhd(x,v,h/2,m);
[x] = driftwhd(x,v,h/2,m);

end




function [x,v] = whdw(x,v,h,m)
% n = 1 is the sun here.  Using barycentric
% velocities.
global GNEWT
n = size(m,2);
gm = GNEWT*m(1);
% Update x_0 only
[x] = driftwhdw(x,v,h/2,m);
% Update v's of planets
[v] = kickwhd(x,v,h/2,m);

% Begin B step
psum0 = zeros(3,1);
% total initial momentum
for i=1:n
    psum0(:) = psum0(:) + m(i)*v(:,i);
end
psum = zeros(3,1);
% update x,v of planets
for i=2:n
    xo = x(:,i)-x(:,1);
    vo = v(:,i);
    [xn,vn] = kepler_stepxv(gm,xo,vo,h);
    x(:,i) = xn(:)+x(:,1);
    v(:,i) = vn(:);
    psum(:) = psum(:) + m(i)*v(:,i);
end
% update v_0 only
psum = zeros(3,1);
% final momentum without Sun
for i=2:n
    psum(:) = psum(:) + m(i)*v(:,i);
end
% conserve total momentum
v(:,1) = (psum0(:) -psum(:))/m(1);
% End B step

[v] = kickwhd(x,v,h/2,m);
[x] = driftwhdw(x,v,h/2,m);

end


function [x,v] = h2pert(xp,vp,m,isun)
global year
npart = size(xp,2);
sumx = zeros(3,1);
sumv = zeros(3,1);
for i=1:npart
    if(i == isun)
        continue
    end
    sumx(1:3) = sumx(1:3) - xp(1:3,i)*m(i);
    sumv(1:3) = sumv(1:3) - vp(1:3,i)*m(i);
end
x(1:3,isun) = sumx/sum(m);
v(1:3,isun) = sumv/sum(m);
for i=1:npart
    if(i == isun)
        continue
    end
    x(1:3,i) = x(1:3,isun) + xp(1:3,i);
    v(1:3,i) = v(1:3,isun) + vp(1:3,i);
end

end










function [x] = drift(x,v,h)
n = size(x,2);
for i=1:n
    x(:,i) = x(:,i) + h*v(:,i);
end

end

function [x] = driftwhd(x,v,h,m)
% Assuming sun has index 1.
n = size(x,2);
sump = zeros(3,1);
msun = m(1);
for i=2:n
    sump(:) = sump(:) + m(i)*v(:,i);
end
for i=2:n
    x(:,i) = x(:,i) + h/msun*sump(:);
end

end

function [x] = driftwhdw(x,v,h,m)
% Assuming sun has index 1.

x(:,1) = x(:,1) + h*v(:,1);


end


function [x,v] = keplerij(m,x,v,i,j,h)
global GNEWT
    xr = x(:,i)-x(:,j);
    vr = v(:,i)-v(:,j);
    gm = GNEWT*(m(i) + m(j));
    if(gm == 0)
        x(:,i) = x(:,i) + h*v(:,i);
        x(:,j) = x(:,j) + h*v(:,j);
    else
    %   xn,vn,xr,vr are 3*1 matrices.
        [xn,vn] = kepler_stepxv(gm,xr,vr,h);
        delx = xn - xr;
        delv = vn - vr;
        %go back absolute coord, add cm motion
        [x,v] = centerm(m,x,v,delx,delv,i,j,h);
    end
end










function [r,th,pr,pth] = phasespace(x,v)
r = sqrt(sum(x.^2));
pr = sum((x.*v))/r;
th = atan(x(2)/x(1));
pth = r^2/(1+(x(2)/x(1))^2)*(v(2)/x(1) - x(2)/x(1)^2*v(1));



end





% 
% function [omtildc,omc,ic,em,Lc,a] = calcorbinc(Q,P,m)
% global GNEWT
% [q,p] = convert2cart(m,Q,P);
% for i=1:6
%     v(:,i) = p(:,i)/m(i);
% end
% 
% qr = q(:,6)-q(:,1);
% vr = v(:,6)-v(:,1);
% qrm = norm(qr);
% vrm = norm(vr);
% M = m(6) + m(1);
% mu = m(6)*m(1)/M;
% pr = vr*mu;
% mudan = GNEWT*M;
% 
% h = cross(qr,pr);
% e = (cross(pr,h)-GNEWT*mu^2*M*qr/qrm)/(GNEWT*mu^2*M);
% om = atan2(h(1),-h(2));
% i = atan2(sqrt(h(1)^2 + h(2)^2),h(3));
% ecosw = e(1)*cos(om) + e(2)*sin(om);
% esinw = e(3)/sin(i);
% omtild = atan2(esinw,ecosw) + om;
% em = norm(e);
% a = (2/qrm - vrm^2/mudan)^(-1);
% nu = acos(dot(e,qr)/(em*qrm));
% L = omtild + nu;
% 
% 
% % elements omtild, om, i, em, a,L
% convert = 360/(2*pi);
% omtildc = omtild*convert;
% omc = om*convert;
% ic = i*convert;
% Lc = L*convert;
% 
% end
% 
% 
% 





function [theta] = calcinc(Q,P,m)
[q,p] = convert2cart(m,Q,P);

Lpl = cross(q(:,6),p(:,6));
Lsolar = zeros(3,1);
for i=1:5
    Lsolar = Lsolar + cross(q(:,i),p(:,i));
end
theta = acos(dot(Lsolar,Lpl)/(norm(Lsolar)*norm(Lpl)));



end

function [a,em] = calcorb(Q,v,m,ik,jk)
global GNEWT


qr = Q(:,ik)-Q(:,jk);
vr = v(:,ik)-v(:,jk);
qrm = norm(qr);
vrm = norm(vr);
M = m(ik) + m(jk);
mu = m(ik)*m(jk)/M;
pr = vr*mu;
mudan = GNEWT*M;

h = mu*cross(qr,vr);
e = (cross(pr,h)-GNEWT*mu^2*M*qr/qrm)/(GNEWT*mu^2*M);
em = norm(e);
a = (2/qrm - vrm^2/mudan)^(-1);

end


function [E,L] = consq(m,Q,P)
% Q
% P
global GNEWT CASE
n = size(m,2);

if(CASE == 1 || CASE == 2 || CASE == 3)
%     Democratic Heliocentric or Canonical Heliocentric
    M = sum(m);
    if(CASE == 1 || CASE == 2)
        a = m/M;
    else
        a(1) = 1;
        a(2:n) = 0;
    end
    T1 = a(1)^2/(2*m(1));
    T2 = 0;
    T3 = zeros(3,1);
    T4 = zeros(3,1);
    V1 = 0;
    V2 = 0;
    v1 = Q(:,1);
    v2 = P(:,1);
    L = cross(v1,v2);
    for i=2:n
        L(:) = L(:) + cross(Q(:,i),P(:,i));
        T1 = T1 + a(i)^2/(2*m(i));
        T2 = T2 + sum(P(:,i).^2)/(2*m(i));
        T3 = T3 + P(:,i);
        T4 = T4 + (a(i)/m(i) - a(1)/m(1))*P(:,i);
        V1 = V1 - GNEWT*m(i)*m(1)/(sqrt(sum(Q(:,i).^2)));
        for j=i+1:n
            dist = sqrt(sum((Q(:,i)-Q(:,j)).^2));
            V2 = V2 - GNEWT*m(i)*m(j)/dist;
        end
    end
    T1 = T1*sum(P(:,1).^2);
    T3 = sum(T3(:).^2)/(2*m(1));
    T4 = dot(P(:,1),T4(:));
    E = T1+T2+T3+T4+V1+V2;
elseif(CASE == 4)
    [q] = pos2cart(Q,m);
    L = zeros(3,1);
    E = 0;
    for i = 1:n
 %Kinetic Energy part
            L(:) = L(:) + cross(Q(:,i),P(:,i));
            if(i == 1)
                M(i) = m(1);
                mp(i) = sum(m);
            else
                M(i) = M(i-1) + m(i);
                mp(i) = m(i)*M(i-1)/M(i);
            end
            E = E + 1/(2*mp(i))*sum(P(:,i).^2); 
        for j=i+1:n
           xij = q(:,i) - q(:,j);
           rij = norm(xij);
           E = E - GNEWT*m(i)*m(j)/rij;
        end
    end
     
elseif(CASE == 5 || CASE == 6 || CASE == 7 || CASE == 8 || CASE == 9)
    L = zeros(3,1);
    E = 0;
    for i = 1:n
        L(:) = L(:) + cross(Q(:,i),P(:,i));
        E = E + 1/(2*m(i))*sum(P(:,i).^2); 
        for j=i+1:n
           xij = Q(:,i) - Q(:,j);
           rij = sqrt(sum(xij.^2));
           E = E - GNEWT*m(i)*m(j)/rij;
        end
    end
    
    
end

    





end

  



function [E,L] = consqv(m,Q,v)

global GNEWT CASE
n = size(m,2);
% %     Democratic Heliocentric or Canonical Heliocentric
[x,junk] = convert2cart(m,Q,v);
if(CASE == 1)
    T = 0;
    V = 0;
    L = 0;
    for i=1:n
        L = L + m(i)*cross(x(:,i),v(:,i));
        T = T + 1/2*m(i)*sum(v(:,i).^2);
        for j=i+1:n
            Qij = sqrt(sum((x(:,i)-x(:,j)).^2));
            V = V - GNEWT*m(i)*m(j)/Qij;
        end
    end
end
E = T + V;
L = L;


% 
% % 
% if(CASE == 1)
%     n = size(m,2);
% %     Democratic Heliocentric or Canonical Heliocentric
%     M = sum(m);
%     a = m/M;
%     [junk,P] = convertcart(m,Q,v); %only need P from Cartesian
%     T1 = a(1)^2/(2*m(1));
%     T2 = 0;
%     T3 = zeros(3,1);
%     T4 = zeros(3,1);
%     V1 = 0;
%     V2 = 0;
%     v1 = Q(:,1);
%     v2 = P(:,1);
%     L = cross(v1,v2);
%     for i=2:n
%         L(:) = L(:) + cross(Q(:,i),P(:,i));
%         T1 = T1 + a(i)^2/(2*m(i));
%         T2 = T2 + sum(P(:,i).^2)/(2*m(i));
%         T3 = T3 + P(:,i);
%         T4 = T4 + (a(i)/m(i) - a(1)/m(1))*P(:,i);
%         V1 = V1 - GNEWT*m(i)*m(1)/(sqrt(sum(Q(:,i).^2)));
%         for j=i+1:n
%             dist = sqrt(sum((Q(:,i)-Q(:,j)).^2));
%             V2 = V2 - GNEWT*m(i)*m(j)/dist;
%         end
%     end
%     T1 = T1*sum(P(:,1).^2);
%     T3 = sum(T3(:).^2)/(2*m(1));
%     T4 = dot(P(:,1),T4(:));
%     E = T1+T2+T3+T4+V1+V2;
% 
% 
% 
% end
% 


if(CASE == 4)
    [q] = pos2cart(Q,m);
    L = zeros(3,1);
    E = 0;
    for i = 1:n
 %Kinetic Energy part
            if(i == 1)
                M(i) = m(1);
                mp(i) = sum(m);
            else
                M(i) = M(i-1) + m(i);
                mp(i) = m(i)*M(i-1)/M(i);
            end
            L(:) = L(:) + cross(Q(:,i),mp(i)*V(:,i));
            E = E + 1/2*mp(i)*sum(V(:,i).^2); 
        for j=i+1:n
           xij = q(:,i) - q(:,j);
           rij = norm(xij);
           E = E - GNEWT*m(i)*m(j)/rij;
        end
    end
     
    
end

    





end




function [m,x0,v0,om,a2bod,radius] = intestpart()
global GNEWT
G = GNEWT;
%Distance is in AU, mass in solar masses, t in day
msun = 1;
mu = 0.01;
mjup = mu*msun/(msun-mu);

% %take care of massless m3 IC
n = 3; %2 or 3 particles
m = zeros(1,n);
x0 = zeros(3,n);
v0 = zeros(3,n);

%massive particles IC
a2bod = 5.2;
ejup = 0;
om = sqrt(G*(mjup+msun)/a2bod^3);
m(1) = msun;
m(3) = mjup;
mu0 = (m(1)*m(3))/(m(1)+m(3));
mt0 = m(1)+m(3);
fa0 = m(1)/mt0;
fb0 = m(3)/mt0;
v2bod = sqrt(G*(msun+mjup)/(a2bod)*(1-ejup)/(1+ejup));
% Start jupiter apoapse
x2bod = a2bod*(1+ejup);
x0(1,3) = fa0*x2bod;
v0(2,3) = fa0*v2bod;
% Sun body
x0(1,1) = -fb0*x2bod;
v0(2,1) = -fb0*v2bod;

% %take care of massless m3 IC
xhel = 4.42;
% xhel = 1.0;
vhel = 0.0072;
x0(1,2) = xhel + x0(1,1);
v0(2,2) = vhel + v0(2,1);
% % switch mass indices

vcm = zeros(3,1);
xcm = zeros(3,1);
for j=1:n
    vcm(:) = vcm(:)+m(j)*v0(:,j);
    xcm(:) = xcm(:)+m(j)*x0(:,j);
end
vcm = vcm/sum(m);
xcm = xcm/sum(m);
% Adjust so CoM is stationary; barycentric velocities
for j=1:n
    v0(:,j) = v0(:,j)-vcm(:); 
    x0(:,j) = x0(:,j)-xcm(:);
end

radius = (m(3)/(3*m(1)))^(1/3)*a2bod; %Jupiter Hill radius
% a2bod

end


function [x,v] = perturbic(x,v,i)
delta = 1e-3;
xscale = 10;
vscale = 0.01;

if(i == 1)
    x(1,3) = x(1,3) + xscale*delta;
elseif(i == 2)
    x(1,3) = x(1,3) - xscale*delta;
elseif(i == 3)
    x(2,3) = x(2,3) + xscale*delta;
elseif(i == 4)
    x(2,3) = x(2,3) - xscale*delta;
elseif(i == 5)
    v(1,3) = v(1,3) + vscale*delta;
elseif(i == 6)
    v(1,3) = v(1,3) - vscale*delta;
elseif(i == 7)
    v(2,3) = v(2,3) + vscale*delta;
elseif(i == 8)
    v(2,3) = v(2,3) - vscale*delta;
end
    
    



end




function [a,e,T,th] = cartes2el(gm,x,v,t)
%x is 2D, v is 2D

r = sqrt(sum(x.^2));
l = x(1)*v(2) - x(2)*v(1);
Ax = v(2)*l/gm-x(1)/r;
Ay = -v(1)*l/gm - x(2)/r;
th = atan(Ay/Ax);
vs = sum(v.^2);
u = sum(x.*v);
alpha=2.0*gm/r-vs;
a=gm/alpha;
ec = 1-r/a;
es=u/sqrt(gm*a);
e=sqrt(ec^2+es^2);
E = atan2(es,ec);
en=sqrt(gm/a^3);
M = (E-es)*en;
P = 2*pi/en;
T = t-(E-es)/en;
while(abs(T) > P)
   T = T + (-1)*sign(T)*P; 
end

end



function [jac] = jacobi(m,Q,v,om,asemi,t)
global GNEWT

n = size(Q,2);
M = sum(m);
a = m/M;
s1 = zeros(3,1);
for i=2:n
    s1 = s1 + a(i)*Q(:,i);
end
x(:,1) = Q(:,1) - s1;
for i=2:n
   x(:,i) = Q(:,1) + Q(:,i) - s1; 
end



i = 3;  
cmt = cos(om*t);
smt = sin(om*t);
xrot = x(1,i)*cmt + x(2,i)*smt;
yrot = -x(1,i)*smt + x(2,i)*cmt;
xdrot = v(1,i)*cmt+v(2,i)*smt + ...
    om*(-x(1,i)*smt + x(2,i)*cmt);
ydrot = -v(1,i)*smt + v(2,i)*cmt + ...
    om*(-x(1,i)*cmt - x(2,i)*smt);
mu1 = m(2)/(m(1)+m(2));
mu2 = m(1)/(m(1)+m(2));
r1 = sqrt((xrot + mu1*asemi)^2 + yrot^2);
r2 = sqrt((xrot - mu2*asemi)^2 + yrot^2);
U = -1/2*om^2*(xrot^2 + yrot^2) - GNEWT*m(1)/r1 - ...
    GNEWT*m(2)/r2;
jac = 1/2*(xdrot^2 + ydrot^2) + U; 

end






function [xrot2,yrot2] = rotatecoord(Q,t,om,m)
n = size(Q,2);
M = sum(m);
a = m/M;
s1 = zeros(3,1);
for i=2:n
    s1 = s1 + a(i)*Q(:,i);
end
x(:,1) = Q(:,1) - s1;
for i=2:n
   x(:,i) = Q(:,1) + Q(:,i) - s1; 
end

i = 3;
cmt = cos(om*t);
smt = sin(om*t);
xrot1 = x(1,i)*cmt + x(2,i)*smt;
yrot1 = -x(1,i)*smt + x(2,i)*cmt;

i = 2;
xrot2 = x(1,i)*cmt + x(2,i)*smt;
yrot2 = -x(1,i)*smt + x(2,i)*cmt;

i = 1;
xrot3 = x(1,i)*cmt + x(2,i)*smt;
yrot3 = -x(1,i)*smt + x(2,i)*cmt;

end


function [qrot] = rotatecoordjac(q,t,om,m)
n = size(q,2);

cmt = cos(om*t);
smt = sin(om*t);

rotmat = [cmt smt 0; -smt cmt 0; 0 0 1;];
for i=1:3
    qrot(:,i) = q(:,i)'*rotmat';
end



end




function [jac] = jacobinonrot(m,Q,v,om)
global GNEWT

G = GNEWT;

% Check eq. (31) of Dehnen & Hernandez paper
n = size(Q,2);
M = sum(m);
a = m/M;
s1 = zeros(3,1);
for i=2:n
    s1 = s1 + a(i)*Q(:,i);
end
x(:,1) = Q(:,1) - s1;
for i=2:n
   x(:,i) = Q(:,1) + Q(:,i) - s1; 
end


% x = z;

rsq = sum(x(:,3).^2);
vsq = sum(v(:,3).^2);
ang = x(1,3)*v(2,3) - x(2,3)*v(1,3);
r1 = sqrt(sum((x(:,3)-x(:,1)).^2));
r2 = sqrt(sum((x(:,3)-x(:,2)).^2));


U = - G*m(1)/r1 - G*m(2)/r2;
jac = 1/2*vsq-om*ang + U;


end



function [jac] = jacobinonrotv(m,Q,V,om)
global GNEWT

G = GNEWT;
[x,v] = convert2cartv(m,Q,V);
prim = 1;
sec = 2;
ind = 3;

rsq = sum(x(:,ind).^2);
vsq = sum(v(:,ind).^2);
ang = x(1,ind)*v(2,ind) - x(2,ind)*v(1,ind);
r1 = sqrt(sum((x(:,ind)-x(:,prim)).^2));
r2 = sqrt(sum((x(:,ind)-x(:,sec)).^2));


U = - G*m(prim)/r1 - G*m(sec)/r2;
jac = 1/2*vsq-om*ang + U;


end



% From MERCURY ANALYSIS PROGRAMS
% function [jac] = jacobitrem(m,x,v,om,a,t)
% global G
% pretty sure these are barycentric
% n = size(x,2);
% 
%         rsq = sum(x(:,3).^2);
%         vsq = sum(v(:,3).^2);
%         ang = x(1,3)*v(2,3) - x(2,3)*v(1,3);
%         r1 = sqrt(sum((x(:,3)-x(:,1)).^2));
%         r2 = sqrt(sum((x(:,3)-x(:,2)).^2));
%         
%      
%         U = - G*m(1)/r1 - G*m(2)/r2;
%         jac = 1/2*vsq-om*ang + U; 
% 
% 
% end






function [Q,P] = map(Q,P,h,m)
global CASE GNEWT
n = size(m,2);

if(CASE == 1)
%     Same for WHD and Canonical Heliocentric maps
%     map B for h/2
    [Q,P] = mapbd(Q,P,h/2,m);
%     map A for time h
    [Q,P] = mapad(Q,P,h,m);
    %     map B for h/2
    [Q,P] = mapbd(Q,P,h/2,m);
    
%     [Q,P] = mapad(Q,P,h/2,m);
elseif(CASE == 2)
    %     map B for h/2
    [Q,P] = mapbdt(Q,P,h/2,m);
%     map A for time h
    [Q,P] = mapadt(Q,P,h,m);
    %     map B for h/2
    [Q,P] = mapbdttrans(Q,P,h/2,m);
elseif(CASE == 3)
    [Q,P] = mapbhel(Q,P,h/2,m);
%     map A for time h
    [Q,P] = mapad(Q,P,h,m);
    %     map B for h/2
    [Q,P] = mapbhel(Q,P,h/2,m);
    
elseif(CASE == 4)
%             map B for h/2
    [Q,P] = mapbj(Q,P,h/2,m);
    
%     map A for time h
    [Q,P] = mapaj(Q,P,h,m);
    %     map B for h/2
    [Q,P] = mapbj(Q,P,h/2,m);
    
%     [Q,P] = mapaj(Q,P,h/2,m);
%     
%     [Q,P] = mapbj(Q,P,h,m);
%     
%     [Q,P] = mapaj(Q,P,h/2,m);
elseif(CASE == 5)
    [Q,P] = mapbc(Q,P,h/2,m);
    
    [Q,P] = mapac(Q,P,h,m);

    [Q,P] = mapbc(Q,P,h/2,m);
    
elseif(CASE == 6)
    [Q,P] = mapbct(Q,P,h/2,m);
    
    [Q,P] = mapact(Q,P,h,m);

    [Q,P] = mapbcttrans(Q,P,h/2,m);
    
elseif(CASE == 7 || CASE == 8)
    for i=1:n
       v(:,i) = P(:,i)/m(i); 
       if(CASE == 7)
           for j=1:n
                if(i == 1 || j == 1) 
                    pair(i,j) = 0;
                    pair(j,i) = 0;
                else
                    pair(i,j) = 1;
                    pair(j,i) = 1;
                end
           end
       elseif(CASE == 8)
           for j=1:n 
                pair(i,j) = 0;
                pair(j,i) = 0;
           end
       end
    end
    [Q,v] = hb15(Q,v,h,m,pair);
    for i=1:n
        P(:,i) = v(:,i)*m(i);
    end
    
elseif(CASE == 9)
    [Q,P] = mapbl(Q,P,h/2,m);
    
    [Q,P] = mapal(Q,P,h,m);
    
    [Q,P] = mapbl(Q,P,h/2,m);
    
    
end



end










function [Q,v] = mapv(Q,v,h,m)
global CASE 


if(CASE == 1)
%     Same for WHD and Canonical Heliocentric maps
%     map B for h/2
    [Q,v] = mapbdv(Q,v,h/2,m);
%     map A for time h
    [Q,v] = mapadv(Q,v,h,m);
    %     map B for h/2
    [Q,v] = mapbdv(Q,v,h/2,m);
    
%     [Q,P] = mapad(Q,P,h/2,m);
end

    
if(CASE == 4)
%             map B for h/2
    [Q,V] = mapbjv(Q,V,h/2,m);
%     [E0,L0] = consqv(m,Q,V);
%     E0
%     blergh

%     map A for time h
    [Q,V] = mapajv(Q,V,h,m);
%         [E0,L0] = consqv(m,Q,V);
%     E0
%     blergh
    %     map B for h/2
    [Q,V] = mapbjv(Q,V,h/2,m);
    
end



end





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% BEGIN MERCURIUS CODE K(X) V, barycentric velocities
function [Q,v] = mapmercuriusv(Q,v,h,m)
  
    [Q] = map1mercuriusv(Q,v,h/2,m);
    [v] = map2mercuriusv(Q,v,h/2,m);
    [Q,v] = map3mercuriusv(Q,v,h,m);
    [v] = map2mercuriusv(Q,v,h/2,m);
    [Q] = map1mercuriusv(Q,v,h/2,m);
    
end

function [Q] = map1mercuriusv(Q,v,h,m)
% What I put in for Q below doesn't affect P calculation
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end

for i=2:n
    spk = sp - P(:,i);
    Q(:,i) = Q(:,i) + h/(m(1))*spk;
end


end

function [v] = map2mercuriusv(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 RADIUS
n = size(m,2);
R = RADIUS; %constant for now

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        x = (qmod - HTRANS1*R)/(HTRANS2*R);
        [Kij,Kijp] = Kcalc(x);
        Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
        faci = GNEWT*m(j)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        facj = GNEWT*m(i)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end


end

function [Q,v] = map3mercuriusv(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);
    
% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end
    

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
        if(qmod < HTRANS3*R)
            a = Q(:,i);
            b = Q(:,j);
            c = v(:,i);
            d = v(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverv(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            v(:,i) = z(7:9,2);
            v(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
        end
    end
end
for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
        mu = m(i)*m(1)/(m(i) + m(1));
        gm = GNEWT*(m(1) + m(i));
        a = Q(:,i);
        b = v(:,i)*(m(i) + m(1))/m(1);
       
        [a,b] = kepler_stepxv(gm,a,b,h);
        fact = m(1)/(m(i) + m(1));
        v(:,i) = b(:) * fact;
        Q(:,i) = a(:);
       
    %     check did not miss close encounter with anyone (should be fast calc); K did not become less than 1
        for j=2:n
            if(i == j)
                continue
            end
            qvec = Q(:,i) - Q(:,j);
            qmod = sqrt(sum(qvec.^2));
            x = (qmod - HTRANS1*R)/(HTRANS2*R);
            [Kij,Kijp] = Kcalc(x);
            if(Kij < 1)
                errormissedenc
            end
        end
    
    
    end
end
% Calculate solar velocity, (assume CoM frame?)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    

% Q

end

function xprime = odesolverv(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
x = (qmod - HTRANS1*R)/(HTRANS2*R);
[Kij,Kijp] = Kcalc(x);
Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
faci = -GNEWT*m(1)/qmagi^3*Qi;
facj = -GNEWT*m(1)/qmagj^3*Qj;
fac = GNEWT/qmod^3*qvec* ...
        (1 - Kij + qmod*Kijp);
facip = (m(i) + m(1))/m(1);
facjp = (m(j) + m(1))/m(1);
xprime = [vi(1)*facip; vi(2)*facip; vi(3)*facip;
    vj(1)*facjp; vj(2)*facjp; vj(3)*facjp; ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end







function xprime = odesolvervjac(t,z,m,M,Ms,mp,i,j)
    global GNEWT 
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
Vi(:) = z(7:9);
Vj(:) = z(10:12);


eps1 = M(i-1)/M(i);
distv = Qj - eps1*Qi;
dist3 = norm(distv).^3;
% 
%j terms
T1 = -GNEWT*m(i)/dist3*distv*M(i+1)/M(i);
T2 = -GNEWT*Ms(j)/(norm(Qj))^3*Qj;

%i terms
T3 = -GNEWT*m(j)/dist3*(Qi - Qj);
T4 = -GNEWT*Ms(i)/(norm(Qi))^3*Qi;

% %j terms
% T1 = -GNEWT*m(i)*m(j)/dist3*distv/mp(j);
% T2 = -GNEWT*Ms(j)/(norm(Qj))^3*Qj;
% 
% %i terms
% T3 = -GNEWT*m(j)/dist3*(Qi - Qj);
% T4 = -GNEWT*Ms(i)/(norm(Qi))^3*Qi;

xprime = [Vi(1); Vi(2); Vi(3);
    Vj(1); Vj(2); Vj(3); ...
    T3(1) + T4(1); T3(2) + T4(2); T3(3) + T4(3); ...
    T1(1) + T2(1); T1(2) + T2(2); T1(3) + T2(3);];
for i=1:12
    if(isnan(xprime(i)) == 1)
        xprime
        T1
        T2
        dist3
        distv
        mp
        m
        blergh
    end
end
% xprime

    
end




function xprime = odesolvervjacp(t,z,m,M,Ms,mp,i,j)
    global GNEWT 

n = size(m,2);
vdot = zeros(3,n);
Q(:,1) = z(1:3); %j > i
Q(:,2) = z(4:6);
Q(:,3) = z(7:9);
V(:,1) = z(10:12);
V(:,2) = z(13:15);
V(:,3) = z(16:18);


[q] = pos2cart(Q,m);
for i=1:n
    for j=i+1:n
        xij = q(:,i)-q(:,j);
        rij = sqrt(sum(xij.^2));
        fac = GNEWT/rij^3*xij; 
        vdot(:,i) = vdot(:,i) - fac*m(j);
        vdot(:,j) = vdot(:,j) + fac*m(i);
    end
end
[Vdot] = mom2jacv(vdot,m);




xprime = [V(1,1); V(2,1); V(3,1);...
    V(1,2); V(2,2); V(3,2); ...
    V(1,3); V(2,3); V(3,3); ...
    Vdot(1,1); Vdot(2,1); Vdot(3,1);...
    Vdot(1,2); Vdot(2,2); Vdot(3,2);...
    Vdot(1,3); Vdot(2,3); Vdot(3,3);];


    
end






function xprime = odesolvervtot(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));

faci = -GNEWT*m(1)/qmagi^3*Qi;
facj = -GNEWT*m(1)/qmagj^3*Qj;
fac = GNEWT/qmod^3*qvec;
facip = (m(i) + m(1))/m(1);
facjp = (m(j) + m(1))/m(1);
xprime = [vi(1)*facip; vi(2)*facip; vi(3)*facip;
    vj(1)*facjp; vj(2)*facjp; vj(3)*facjp; ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end














% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% END MERCURIUS CODE V




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% BEGIN MERCURIUS CODE K(X)
function [Q,P] = mapmercurius(Q,P,h,m)
 
    [Q] = map1mercurius(Q,P,h/2,m);
    [P] = map2mercurius(Q,P,h/2,m);
    [Q,P] = map3mercurius(Q,P,h,m);
    [P] = map2mercurius(Q,P,h/2,m);
    [Q] = map1mercurius(Q,P,h/2,m);
    
end

function [Q] = map1mercurius(Q,P,h,m)
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end

for i=2:n
    spk = sp - P(:,i);
    Q(:,i) = Q(:,i) + h/(m(1))*spk;
end


end

function [P] = map2mercurius(Q,P,h,m)
global GNEWT HTRANS1 HTRANS2 RADIUS
n = size(m,2);
R = RADIUS; %constant for now

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        x = (qmod - HTRANS1*R)/(HTRANS2*R);
        [Kij,Kijp] = Kcalc(x);
        Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end

end

function [Q,P] = map3mercurius(Q,P,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
        if(qmod < HTRANS3*R)
            a = Q(:,i);
            b = Q(:,j);
            c = P(:,i);
            d = P(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolver(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            P(:,i) = z(7:9,2);
            P(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
        end
    end
end
for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
        mu = m(i)*m(1)/(m(i) + m(1));
        gm = GNEWT*(m(1) + m(i));
        a = Q(:,i);
        b = P(:,i)/mu;
        [a,b] = kepler_stepxv(gm,a,b,h);
        P(:,i) = b(:)*mu;
        Q(:,i) = a(:);
    end
end
    
% P
% Q

end

function xprime = odesolver(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
Pi(:) = z(7:9);
Pj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
mui = m(i)*m(1)/(m(i) + m(1));
muj = m(j)*m(1)/(m(j) + m(1));
x = (qmod - HTRANS1*R)/(HTRANS2*R);
[Kij,Kijp] = Kcalc(x);
Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
gmi = GNEWT*(m(1) + m(i));
gmj = GNEWT*(m(1) + m(j));
faci = -gmi*mui/qmagi^3*Qi;
facj = -gmj*muj/qmagj^3*Qj;
facij = GNEWT*m(i)*m(j)/qmod^3*qvec* ...
        (1 - Kij + qmod*Kijp);
xprime = [Pi(1)/mui; Pi(2)/mui; Pi(3)/mui;
    Pj(1)/muj; Pj(2)/muj; Pj(3)/muj; ...
    faci(1)- facij(1); faci(2)- facij(2); faci(3)- facij(3); ...
    facj(1)+ facij(1); facj(2)+ facij(2); facj(3)+ facij(3);];
    
end





function [K,Kp] = closeenc(Q)
global HTRANS1 HTRANS2 HTRANS3 RADIUS
R = RADIUS;

% 1e4 AU is large enough
n = size(Q,2);
rmin = 1e4; 
Kk = 0;
Kpk = 0;
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        if(qmod < rmin)
            rmin = qmod;
            x = (qmod - HTRANS1*R)/(HTRANS2*R);
            [K,Kp] = Kcalc(x);
            Kk = K;
            Kpk = Kp;
        end
    end
end


end









function [K,Kp] = Kcalc(x)
global CASEPOLY
% Wisdom switching function
if(x < 0)
    K = 0;
    Kp = 0;
elseif(x > 1)
    K = 1;
    Kp = 0;
else

    if(CASEPOLY == 0)
        
% Begin polynomial function with no quotient
% 0 derivatives
K = x;
Kp = 1;
    elseif(CASEPOLY == 1)
% % 
% 1 derivative
K = -2*x^3 + 3*x^2;
Kp = -6*x^2 + 6*x;

    elseif(CASEPOLY == 2)
% 
% 2 derivatives defined
K = 6*x^5 - 15*x^4 + 10*x^3;
Kp = 30*x^4 - 60*x^3 + 30*x^2;

    elseif(CASEPOLY == 3)
% 
% 3 derivatives defined
K = -20*x^7 + 70*x^6 - 84*x^5 + 35*x^4;
Kp = -140*x^6 + 420*x^5 - 420*x^4 + 140*x^3;

    elseif(CASEPOLY == 4)

% 4 derivatives defined
K = 70*x^9 - 315*x^8 + 540*x^7 - 420*x^6 + 126*x^5;
Kp = 630*x^8 - 2520*x^7 + 3780*x^6 - 2520*x^5 + 630*x^4;

    elseif(CASEPOLY == 5)
% 5 derivatives defined
K = -252*x^(11) + 1386*x^(10) - 3080*x^9 + 3465*x^8 - 1980*x^7 + 462*x^6;
Kp = -2772*x^(10) + 13860*x^9 - 27720*x^8 + 27720*x^7 - 13860*x^6 + 2772*x^5;
    end
end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% END MERCURIUS CODE








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% BEGIN MERCURIUS CODE WHD 

function [Q,P] = mapmercuriusWHD(Q,P,h,m)
    [Q,P] = map1mercuriusWHD(Q,P,h/2,m);
    [Q,P] = map2mercuriusWHD(Q,P,h,m);
    [Q,P] = map1mercuriusWHD(Q,P,h/2,m);
    
end



function [Q,P] = map1mercuriusWHD(Q,P,h,m)
global GNEWT HTRANS1 HTRANS2 HATRANS3 RADIUS
R = RADIUS;
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        x = (qmod - HTRANS1*R)/(HTRANS2*R);
        [Kij,Kijp] = Kcalc(x);
        Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end

end


function [Q,P] = map2mercuriusWHD(Q,P,h,m)
global GNEWT HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
        if(qmod < HTRANS3*R)
            a = Q(:,i);
            b = Q(:,j);
            c = P(:,i);
            d = P(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverWHD(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            P(:,i) = z(7:9,2);
            P(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
        end
    end
end
for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
%         flagi(i)
%         i
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = P(:,i)/m(i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    P(:,i) = b(:)*m(i);
    Q(:,i) = a(:);
    end
end


end


function xprime = odesolverWHD(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
Pi(:) = z(7:9);
Pj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
x = (qmod - HTRANS1*R)/(HTRANS2*R);
[Kij,Kijp] = Kcalc(x);
% Kij
% Kijp
Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
% Only gravitating, reduced masses, change for gmi,gmj
gmi = GNEWT*m(1);
gmj = GNEWT*m(1);
faci = -gmi*m(i)/qmagi^3*Qi;
facj = -gmj*m(j)/qmagj^3*Qj;
% facij stays invariant for WHD, WHDS
facij = GNEWT*m(i)*m(j)/qmod^3*qvec* ...
        (1 - Kij + qmod*Kijp);
xprime = [Pi(1)/m(i); Pi(2)/m(i); Pi(3)/m(i);
    Pj(1)/m(j); Pj(2)/m(j); Pj(3)/m(j); ...
    faci(1)- facij(1); faci(2)- facij(2); faci(3)- facij(3); ...
    facj(1)+ facij(1); facj(2)+ facij(2); facj(3)+ facij(3);];
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% END MERCURIUS CODE WHD





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% BEGIN MERCURIUS CODE WHD V

function F0 = selector(Q,v,h,m) %return integrator selector, <0 or >0 to switch integrator
global HTRANS3 RADIUS GNEWT

n = size(Q,2);
R = RADIUS;
lgn = 1000;
F0 = lgn;
hfac = 30;

% F_{i j} =  sqrt((3 V_{ij}^2+ G (m_i + m_j)/Q_{ij}) / Q_{ij}



for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        vvec = v(:,i) - v(:,j);
        qmod = sqrt(sum(qvec.^2));
        vmod2 = (sum(vvec.^2));
        F1 = sqrt(3*vmod2 + GNEWT*(m(i) + m(j))/qmod)/qmod;
        Ft = 1/F1 - hfac*h;
        if(Ft < F0)
            F0 = Ft;
        end
%         if(F2 < F0)
%             F0 = F2;
% %             F0
%         end
    end
end

end

function [F0,ik,jk] = selectorpl(Q,V,m) %return integrator selector, <0 or >0 to switch integrator
global HTRANS3 RADIUS

R = RADIUS;
n = size(Q,2);
[q,v] = convert2cartv(m,Q,V);
ik = 0;
jk = 0;
F0 = 0;

for i=2:n
    for j=i+1:n
        qm = norm(q(:,i) - q(:,j));
%         qm
%         HTRANS3*RADIUS
%         blergh
        if(qm < HTRANS3*RADIUS)
            if(F0 == 1) %error is two close encounters
                blergh
            end
            F0 = 1;
            ik = i;
            jk = j;
        end

    end

end

% F0 = 0;

end


function F0 = selectorJac(Q,V,m,order) %return integrator selector, <0 or >0 to switch integrator

n = size(Q,2);
[q,v] = convert2cartv(m,Q,V);
delta = 9; %in au
Qmag1 = norm(q(:,2));
Qmag2 = norm(q(:,3));

if(order == 1) %usual order
    delta = Qmag1;
    if(Qmag2 - delta > 0)

        F0 = 1; %usual Jacobi
    else
        F0 = -1; %reversed Jacobi
    end
elseif(order == -1) %reversed order
    delta = Qmag2;
    if(Qmag1 - delta > 0)

        F0 = 1; %usual Jacobi
    else
        F0 = -1; %reversed Jacobi
    end
end


end



function F0 = selectorJacm(Q,V,m,order) %return integrator selector, <0 or >0 to switch integrator

n = size(Q,2);
[q,v] = convert2cartv(m,Q,V);
for i=2:n
    qmag(i) = norm(q(:,i)); %of pairs is n-2
    if(i > 2)
        del = qmag(i) - qmag(i-1);
        if(del > 0)
            F0(i-2) = 1; %this is pair number, order is +1 of this
        else
            F0(i-2) = -1;
        end
    end
end




end




function [Q1,V1,m1,order1] = switchjac(m,Q0,V0,order,F0)

[qt,vt] = convert2cartv(m,Q0,V0);
n = size(qt,2);
q1 = qt;
v1 = vt;
m1 = m;
if(F0 > 0)
    if(order > 0)
        Q1 = Q0;
        V1 = V0;
        m1 = m;
        order1 = order;
        return
    elseif(order < 0)
        q1(:,2) = qt(:,3);
        q1(:,3) = qt(:,2);
        v1(:,2) = vt(:,3);
        v1(:,3) = vt(:,2);
        m1(2) = m(3);
        m1(3) = m(2);
        [Q1,V1] = convertcartv(m1,q1,v1);
        order1 = 1;
        return
    end
elseif(F0 < 0)
    if(order < 0)
        Q1 = Q0;
        V1 = V0;
        m1 = m;
        order1 = order;
        return
    elseif(order > 0)
        q1(:,2) = qt(:,3);
        q1(:,3) = qt(:,2);
        v1(:,2) = vt(:,3);
        v1(:,3) = vt(:,2);
        m1(2) = m(3);
        m1(3) = m(2);
        [Q1,V1] = convertcartv(m1,q1,v1);
        order1 = -1;
        return
    end
end

end











function [Q1,V1,m1,order1] = switchjacm(m,Q0,V0,order,F0)

[qt,vt] = convert2cartv(m,Q0,V0);
n = size(F0,2);

for i=1:n
    if(F0(i) < 0) % switching time (does not work if same particle needs to switch w/ two pairs)
        [q1,v1,m1,order1] = switchind(qt,vt,m,order,i+1,i+2);
    end
end

[Q1,V1] = convertcartv(m1,q1,v1);


end



function [Q1,V1,m1,order1] = switchjacmind(m,Q0,V0,order,F0,ik)

[qt,vt] = convert2cartv(m,Q0,V0);
[q1,v1,m1,order1] = switchind(qt,vt,m,order,ik+1,ik+2);

[Q1,V1] = convertcartv(m1,q1,v1);


end



function [q1,v1,m1,order1] = switchind(qt,vt,m,order,i0,j0)

% order
q1 = qt;
v1 = vt;
m1 = m;
order1 = order;

q1(:,i0) = qt(:,j0);
q1(:,j0) = qt(:,i0);
v1(:,i0) = vt(:,j0);
v1(:,j0) = vt(:,i0);
m1(i0) = m(j0);
m1(j0) = m(i0);
order1(i0) = order(j0);
order1(j0) = order(i0);
% order1


end




function [Q1,P1,m1,order1] = switchjacp(m,Q0,P0,order,F0)

[qt,pt] = convert2cartv(m,Q0,P0);
n = size(qt,2);
q1 = qt;
p1 = pt;
m1 = m;
if(F0 > 0)
    if(order > 0)
        Q1 = Q0;
        P1 = P0;
        m1 = m;
        order1 = order;
        return
    elseif(order < 0)
        q1(:,2) = qt(:,3);
        q1(:,3) = qt(:,2);
        p1(:,2) = pt(:,3);
        p1(:,3) = pt(:,2);
        m1(2) = m(3);
        m1(3) = m(2);
        [Q1,P1] = convertcartv(m1,q1,p1);
        order1 = 1;
        return
    end
elseif(F0 < 0)
    if(order < 0)
        Q1 = Q0;
        P1 = P0;
        m1 = m;
        order1 = order;
        return
    elseif(order > 0)
        q1(:,2) = qt(:,3);
        q1(:,3) = qt(:,2);
        p1(:,2) = pt(:,3);
        p1(:,3) = pt(:,2);
        m1(2) = m(3);
        m1(3) = m(2);
        [Q1,P1] = convertcartv(m1,q1,p1);
        order1 = -1;
        return
    end
end

end



function [Q,P,m,F,casen,order] = mapAnd(Q0,P0,h,m,F0,order0)

if(F0 > 0) %usual order (M1)
    [Q1,P1,m1,order1] = switchjacp(m,Q0,P0,order0,F0); %switch indices if needed, map usual
    [Q1,P1] = mapv(Q1,P1,h,m1);
    F1 = selectorJac(Q1,P1,m1,order1);
    if(F1 > 0) %map usual order still correct
        Q = Q1;
        P = P1;
        m = m1;
        F = F1;
        order = order1;
        casen = 1;
        return
    end
end
%map with reverse order (M2)
[Q1,P1,m1,order1] = switchjacp(m,Q0,P0,order0,-1); %switch indices if needed, map reversed now
[Q1,P1] = mapv(Q1,P1,h,m1);
F1 = selectorJac(Q1,P1,m1,order1);
Q = Q1;
P = P1;
m = m1;
F = F1;
order = order1;
casen = 2;
return
    
end








function [Q,V,m,F,casen,order,flgr] = mapAndm(Q0,V0,h,m0,F0,order0)
[flg0,ik] = tstf(F0); %test if require order switching, ik gives M2 map
if(flg0 == 0) %M1 map, no switching
    [Q1,V1] = mapv(Q0,V0,h,m0);
    F1 = selectorJacm(Q1,V1,m0,order0);
    [flg1,ik] = tstf(F1); %this sets M2 map, does not need to equal ik above
    if(flg1 == 0) %M1 map, no switching, was fine
        Q = Q1;
        V = V1;
        m = m0;
        F = F1;
        order = order0;
        casen = 1;
        flgr=0;
        return
    end
end
%use M2 with a switch, either ik from before step or test step
[Q,V,m,order] = switchjacmind(m0,Q0,V0,order0,F0,ik);
[Q,V] = mapv(Q,V,h,m);
F = selectorJacm(Q,V,m,order);
casen = 2;
flgr = 1;
return

end




function [x,v,F,intyosh] = mapwalter(Q0,v0,h,m,F0)
global KTP1 KTP2 KTP3
dM = false;
intyosh = true;
C0 = F0 < 0;
% Define F0 = qmod - HTRANS3*R: if small, close encounter

if(C0)
    [xtry,vtry] = mapmercuriusWHDva(Q0,v0,h,m);
    intyosh = true;
else
    [xtry,vtry] = mapmercuriusWHDvb(Q0,v0,h,m);
    intyosh = false;
end
Ftry = selector(xtry,vtry,h,m);
Ctry = (1/2*(F0 + Ftry) < 0);
if(Ctry == C0)
    x = xtry;
    v = vtry;
    F = Ftry;
    KTP1 = KTP1+1;
% 
%     if(C0)
%         [xtry2,vtry2] = mapmercuriusWHDvb(Q0,v0,h,m);
%     else
%         [xtry2,vtry2] = mapmercuriusWHDva(Q0,v0,h,m);
%     end
%     Ftry2 = selector(xtry2,vtry2,h,m);
%     Ctry2 = (1/2*(F0 + Ftry2) < 0);
% %     Ctry2
% %     Ctry
%     if(Ctry2 ~= C0)
%         t
%         Ftry2
%         Ftry
%         Ctry2
%         C0
%         blergh
%     end
    return
end %this should cover most cases and save cost, just one integration

if(C0)
    [xalt,valt] = mapmercuriusWHDvb(Q0,v0,h,m);
    intyosh = false;
else
    [xalt,valt] = mapmercuriusWHDva(Q0,v0,h,m);
    intyosh = true;
end
Falt = selector(xalt,valt,h,m);
Calt = (1/2*(F0 + Falt) < 0);
if(Calt ~= C0)
    x = xalt;
    v = valt;
    F = Falt;
    KTP2 = KTP2 + 1;
    return
end
%final inconsistent cases
KTP3 = KTP3 + 1;
if(dM == C0)
    x = xalt;
    v = valt;
    F = Falt;
    return
else
    x = xtry;
    v = vtry;
    F = Ftry;
    return
end

end






function [Q,V,m,F,order] = mapbasic(Q0,V0,h,m,F0,order0)
if(F0 > 0) %normal order
    [Q,V,m,order] = switchjacp(m,Q0,V0,order0,F0); %switch indices if needed
    [Q,V] = mapv(Q,V,h,m);
    F = selectorJac(Q,V,m,order);
elseif(F0 < 0) %reversed order
    [Q,V,m,order] = switchjacp(m,Q0,V0,order0,F0); %switch indices if needed
    [Q,V] = mapv(Q,V,h,m);
    F = selectorJac(Q,V,m,order);
end




end

function [Q,V,m,F,order,flg,Fpl,ik,jk] = mapbasicm(Q0,V0,h,m,F0,order0,Fpl,ik,jk)

% Test if need switch Jacobi indices, do if needed
[flg,ik2] = tstf(F0);
if(flg == 1)
    [Q,V,m,order] = switchjacm(m,Q0,V0,order0,F0);
else

    Q = Q0;
    V = V0;
    m = m;
    order = order0;
end

%test if need close encounter or usual map
if(Fpl == 0)
    [Q,V] = mapjaca(Q,V,h,m);
elseif(Fpl == 1)
    [Q,V] = mapjacb(Q,V,h,m,ik,jk);
end

F = selectorJacm(Q,V,m,order);
[Fpl,ik,jk] = selectorpl(Q,V,m);


end





function [Q,V,m,F,order,flg,Fpl,ik,jk,casen] = mapbasicmrevers(Q0,V0,h,m,F0,order0,Fpl,ik,jk)

% Test if need switch Jacobi indices, do if needed
[flg,ik2] = tstf(F0);
if(flg == 1)
    [Q,V,m,order] = switchjacm(m,Q0,V0,order0,F0);
else

    Q = Q0;
    V = V0;
    m = m;
    order = order0;
end
% Fpl = 1;

%test if need close encounter or usual map
if(Fpl == 0)
%     Fpl
    [Q1,V1] = mapjaca(Q,V,h,m);
    [Fpl1,ik1,jk1] = selectorpl(Q,V,m);
    if(Fpl == 0)
        Q = Q1;
        V = V1;
        Fpl = Fpl1;
        ik = ik1;
        jk = jk1;
        F = selectorJacm(Q,V,m,order);
        casen = 1;
        return
    end
end
% all else use close encounter map
[Q,V] = mapjacb(Q,V,h,m,ik,jk);
[Fpl,ik,jk] = selectorpl(Q,V,m);
F = selectorJacm(Q,V,m,order);
casen = 2;

return


end






function [flg,ik] = tstf(F0)
flg = 0;
n = size(F0,2);
ik = 0;
for i=1:n
    if(F0(i) < 0) %assume only one switch happens at once
        if(flg == 0)
            flg = 1;
            ik = i;
        elseif(flg == 1)
            blergh;
        end
    end
end



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solar close encounters

% Not close with Sun (Kepler)

function [Q,v] = mapWHa(Q,P,h,m)
    [Q,P] = map1WHa(Q,P,h/2,m);
    [Q,P] = map2WHa(Q,P,h,m);
    [Q,P] = map1WHa(Q,P,h/2,m);

end

% Close with Sun

function [Q,v] = mapWHb(Q,P,h,m)
    [Q,P] = map1WHb(Q,P,h/2,m);
    [Q,P] = map2WHb(Q,P,h,m);
    [Q,P] = map1WHb(Q,P,h/2,m);

end




function [Q,P] = map1WHa(Q,P,h,m)
global GNEWT 
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        faci = GNEWT*m(j)/qmod^3*qvec;
        facj = GNEWT*m(i)/qmod^3*qvec;
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end







function [Q,P] = map2WHa(Q,P,h,m)
global GNEWT

n = size(m,2);
for i=2:n
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = P(:,i)/m(i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    P(:,i) = b(:)*m(i);
    Q(:,i) = a(:);
end




end





function [Q,P] = map1WHb(Q,P,h,m)
global GNEWT 
n = size(m,2);


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        faci = GNEWT*m(j)/qmod^3*qvec;
        facj = GNEWT*m(i)/qmod^3*qvec;
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end












function [Q,v] = map2WHb(Q,v,h,m)
global GNEWT 

n = size(m,2);
sp = zeros(3,1);
for i = 2:n
    sp = sp + P(:,i);
end
j = 0;
for i=2:n
        a = Q(:,i);
        b = P(:,i);
        z0 = [a; b];
%             check z0 is a column vector
        tspan = [0 h];
%             Only two outputs...
        [z, info] = BulirschStoer(@(t,z) odesolverWHDclose(t,z,m,i,sp), tspan,z0,TOL);
        Q(:,i) = z(1:3,2);
        P(:,i) = z(4:6,2);  
    end
end










function xprime = odesolverWHDclose(t,z,m,i,sp)
    global GNEWT 

Qi(:) = z(1:3);
Pi(:) = z(7:9);
qmagi = sqrt(sum(Qi.^2));
gm = GNEWT*m(1)*m(i);
faci = -gm/qmagi^3*Qi;
xprime = [Pi(1)/m(i) + sp(1)/m(1); Pi(2)/m(i)+sp(2)/m(1); Pi(3)/m(i)+sp(3)/m(1);
    faci(1); faci(2); faci(3);];
    
end








% No close encounter case
function [Q,v] = mapWHDv(Q,v,h,m)
    [Q,v] = map1WHDv(Q,v,h/2,m);
    [Q,v] = map2WHDv(Q,v,h,m);
    [Q,v] = map1WHDv(Q,v,h/2,m);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Q,v] = mapmercuriusWHDv(Q,v,h,m)
    [Q,v] = map1mercuriusWHDv(Q,v,h/2,m);
    [Q,v] = map2mercuriusWHDv(Q,v,h,m);
    [Q,v] = map1mercuriusWHDv(Q,v,h/2,m);

end

% Close encounter case
function [Q,v] = mapmercuriusWHDva(Q,v,h,m)
    [Q,v] = map1mercuriusWHDva(Q,v,h/2,m);
    [Q,v] = map2mercuriusWHDva(Q,v,h,m);
    [Q,v] = map1mercuriusWHDva(Q,v,h/2,m);
end

% No close encounter case
function [Q,v] = mapmercuriusWHDvb(Q,v,h,m)
    [Q,v] = map1mercuriusWHDvb(Q,v,h/2,m);
    [Q,v] = map2mercuriusWHDvb(Q,v,h,m);
    [Q,v] = map1mercuriusWHDvb(Q,v,h/2,m);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Close encounter case
function [Q,v] = mapmercuriusWHDsva(Q,v,h,m)
    [Q,v] = map1mercuriusWHDsva(Q,v,h/2,m);
    [Q,v] = map2mercuriusWHDsva(Q,v,h,m);
    [Q,v] = map1mercuriusWHDsva(Q,v,h/2,m);
end

% No close encounter case
function [Q,v] = mapmercuriusWHDsvb(Q,v,h,m)
    [Q,v] = map1mercuriusWHDsvb(Q,v,h/2,m);
    [Q,v] = map2mercuriusWHDsvb(Q,v,h,m);
    [Q,v] = map1mercuriusWHDsvb(Q,v,h/2,m);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v] = map1mercuriusWHDv(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 RADIUS
R = RADIUS;
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        x = (qmod - HTRANS1*R)/(HTRANS2*R);
        [Kij,Kijp] = Kcalc(x);
        Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
        faci = GNEWT*m(j)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        facj = GNEWT*m(i)/qmod^3*qvec* ...
        (Kij - Kijp*qmod);
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q,v] = map1mercuriusWHDva(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 RADIUS
R = RADIUS;
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


end

function [Q,v] = map1mercuriusWHDsva(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 RADIUS
R = RADIUS;
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v] = map1WHDv(Q,v,h,m)
global GNEWT
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
%         can simplify this...
        faci = GNEWT*m(j)/qmod^3*qvec;
        facj = GNEWT*m(i)/qmod^3*qvec;
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end





function [Q,v] = map1mercuriusWHDvb(Q,v,h,m)
global GNEWT
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
%         can simplify this...
        faci = GNEWT*m(j)/qmod^3*qvec;
        facj = GNEWT*m(i)/qmod^3*qvec;
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end

function [Q,v] = map1mercuriusWHDsvb(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 RADIUS
R = RADIUS;
[junk,P] = convertcart(m,Q,v);
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
%         can simplify this...
        faci = GNEWT*m(j)/qmod^3*qvec;
        facj = GNEWT*m(i)/qmod^3*qvec;
        v(:,i) = v(:,i) - h * faci;
        v(:,j) = v(:,j) + h * facj;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v] = map2mercuriusWHDv(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
        if(qmod < HTRANS3*R)
            a = Q(:,i);
            b = Q(:,j);
            c = v(:,i);
            d = v(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverWHDv(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            v(:,i) = z(7:9,2);
            v(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
        end
    end
end
for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
        gm = GNEWT*m(1);
        a = Q(:,i);
        b = v(:,i);
        [a,b] = kepler_stepxv(gm,a,b,h);
        v(:,i) = b(:);
        Q(:,i) = a(:);



%         check did not miss close encounter; K did not become less than 1
        for j=2:n
            if(i == j)
                continue
            end
            qvec = Q(:,i) - Q(:,j);
            qmod = sqrt(sum(qvec.^2));
            x = (qmod - HTRANS1*R)/(HTRANS2*R);
            [Kij,Kijp] = Kcalc(x);
            if(Kij < 1)
                errormissedenc
            end
        end
    end
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end



function [Q,v] = map2mercuriusWHDsv(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
        if(qmod < HTRANS3*R)
            a = Q(:,i);
            b = Q(:,j);
            c = v(:,i);
            d = v(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverWHDv(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            v(:,i) = z(7:9,2);
            v(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
        end
    end
end
for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
        gm = GNEWT*m(1);
        a = Q(:,i);
        b = v(:,i);
        [a,b] = kepler_stepxv(gm,a,b,h);
        v(:,i) = b(:);
        Q(:,i) = a(:);



%         check did not miss close encounter; K did not become less than 1
        for j=2:n
            if(i == j)
                continue
            end
            qvec = Q(:,i) - Q(:,j);
            qmod = sqrt(sum(qvec.^2));
            x = (qmod - HTRANS1*R)/(HTRANS2*R);
            [Kij,Kijp] = Kcalc(x);
            if(Kij < 1)
                errormissedenc
            end
        end
    end
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Q,v] = map2mercuriusWHDva(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
            a = Q(:,i);
            b = Q(:,j);
            c = v(:,i);
            d = v(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverWHDva(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            v(:,i) = z(7:9,2);
            v(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
    end
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end



function [Q,v] = map2mercuriusWHDsva(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum(qvec.^2));
            a = Q(:,i);
            b = Q(:,j);
            c = v(:,i);
            d = v(:,j);
            z0 = [a; b; c; d];
%             check z0 is a column vector
            tspan = [0 h];
%             Only two outputs...
            [z, info] = BulirschStoer(@(t,z) odesolverWHDva(t,z,m,i,j), tspan,z0,TOL);
            Q(:,i) = z(1:3,2);
            Q(:,j) = z(4:6,2);
            v(:,i) = z(7:9,2);
            v(:,j) = z(10:12,2);    
            flagi(i) = 1;
            flagi(j) = 1;
    end
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,v] = map2mercuriusWHDvb(Q,v,h,m)
global GNEWT

n = size(m,2);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end


for i=2:n
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = v(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    v(:,i) = b(:);
    Q(:,i) = a(:);
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end







function [Q,v] = map2WHDv(Q,v,h,m)
global GNEWT

n = size(m,2);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end


for i=2:n
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = v(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    v(:,i) = b(:);
    Q(:,i) = a(:);
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end







function [Q,v] = map2mercuriusWHDsvb(Q,v,h,m)
global GNEWT HTRANS1 HTRANS2 HTRANS3 TOL RADIUS

n = size(m,2);
R = RADIUS;
n = size(m,2);
flagi = zeros(n,1);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end


for i=2:n
%     Assume only one close encounter at a time
    if(flagi(i) ~= 1) %Kepler case when no close encounters at all with i
        gm = GNEWT*m(1);
        a = Q(:,i);
        b = v(:,i);
        [a,b] = kepler_stepxv(gm,a,b,h);
        v(:,i) = b(:);
        Q(:,i) = a(:);
    end
end


% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));
    




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










function xprime = odesolverWHDv(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
x = (qmod - HTRANS1*R)/(HTRANS2*R);
[Kij,Kijp] = Kcalc(x);
Kijp = Kijp/(HTRANS2*R); %derivative w.r.t. x
% Only gravitating, reduced masses, change for gmi,gmj
gm = GNEWT*m(1);
faci = -gm/qmagi^3*Qi;
facj = -gm/qmagj^3*Qj;
% facij stays invariant for WHD, WHDS
fac = GNEWT/qmod^3*qvec* ...
        (1 - Kij + qmod*Kijp);
xprime = [vi(1); vi(2); vi(3);
    vj(1); vj(2); vj(3); ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end



function xprime = odesolverWHDvtot(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
% Only gravitating, reduced masses, change for gmi,gmj
gm = GNEWT*m(1);
faci = -gm/qmagi^3*Qi;
facj = -gm/qmagj^3*Qj;
del = m(2)/m(1)*vi;
% facij stays invariant for WHD, WHDS
fac = GNEWT/qmod^3*qvec;
xprime = [vi(1)+del(1); vi(2)+del(2); vi(3)+del(3);
    vj(1)+del(1); vj(2)+del(2); vj(3)+del(3); ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function xprime = odesolverWHDva(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
% Only gravitating, reduced masses, change for gmi,gmj
gm = GNEWT*m(1);
faci = -gm/qmagi^3*Qi;
facj = -gm/qmagj^3*Qj;
% facij stays invariant for WHD, WHDS
fac = GNEWT/qmod^3*qvec;
xprime = [vi(1); vi(2); vi(3);
    vj(1); vj(2); vj(3); ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end


function xprime = odesolverWHDsva(t,z,m,i,j)
    global GNEWT HTRANS1 HTRANS2 RADIUS
    R = RADIUS;
%     z is 12 X 1 in size; Calculate Pdot Qdot
Qi(:) = z(1:3);
Qj(:) = z(4:6);
vi(:) = z(7:9);
vj(:) = z(10:12);
qvec = Qi - Qj;
qmod = sqrt(sum(qvec.^2));
qmagi = sqrt(sum(Qi.^2));
qmagj = sqrt(sum(Qj.^2));
% Only gravitating, reduced masses, change for gmi,gmj
gm = GNEWT*m(1);
faci = -gm/qmagi^3*Qi;
facj = -gm/qmagj^3*Qj;
% facij stays invariant for WHD, WHDS
fac = GNEWT/qmod^3*qvec;
xprime = [vi(1); vi(2); vi(3);
    vj(1); vj(2); vj(3); ...
    faci(1)- fac(1)*m(j); faci(2)- fac(2)*m(j); faci(3)- fac(3)*m(j); ...
    facj(1)+ fac(1)*m(i); facj(2)+ fac(2)*m(i); facj(3)+ fac(3)*m(i);];
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % 
% END MERCURIUS CODE WHD V





function [q,p] = mapact(q,p,h,m)
global GNEWT
n = size(m,2);

ptot = zeros(3,1);
ppart = zeros(3,1);
for i=1:n
   ptot = ptot + p(:,i); 
end

for i=2:n
    gm = GNEWT*(m(1) + m(i));
    mu = m(i)*m(1)/(m(i) + m(1));
    [a,b] = kepler_stepxv(gm,q(:,i)-q(:,1),p(:,i)/mu,h);
    p(:,i) = b * mu;
    q(:,i) = a + q(:,1);
    ppart = ppart + p(:,i);
end
p(:,1) = ptot - ppart;


end


function [q,p] = mapbct(q,p,h,m)
global GNEWT
n = size(m,2);

for i=2:n
    for j=i+1:n
        qvec = q(:,i) - q(:,j);
        qmod = sqrt(sum(qvec.^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        p(:,i) = p(:,i) - h * fac;
        p(:,j) = p(:,j) + h * fac;
    end
end

q(:,1) = q(:,1) + h*p(:,1)/m(1);

for i=2:n
    q(:,i) = q(:,i) - h/m(1)*p(:,i);
end


end


function [q,p] = mapbcttrans(q,p,h,m)
global GNEWT
n = size(m,2);

q(:,1) = q(:,1) + h*p(:,1)/m(1);

for i=2:n
    q(:,i) = q(:,i) - h/m(1)*p(:,i);
end


for i=2:n
    for j=i+1:n
        qvec = q(:,i) - q(:,j);
        qmod = sqrt(sum(qvec.^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        p(:,i) = p(:,i) - h * fac;
        p(:,j) = p(:,j) + h * fac;
    end
end



end



function [q,p] = mapbc(q,p,h,m)
global GNEWT
n = size(m,2);


for i=2:n
    for j=i+1:n
        xij = q(:,i)-q(:,j);
        rij = sqrt(sum(xij.^2));
        fac = GNEWT*m(i)*m(j)/rij^3*xij; 
        p(:,i) = p(:,i) - h*fac;
        p(:,j) = p(:,j) + h*fac;
    end
end
q(:,1) = q(:,1) + h*p(:,1)/m(1);


end



function [q,p] = mapac(q,p,h,m)
global GNEWT
n = size(m,2);

gm = GNEWT*m(1);
ptot = zeros(3,1);
ppart = zeros(3,1);
for i=1:n
   ptot = ptot + p(:,i); 
end


for i=2:n
    [a,b] = kepler_stepxv(gm,q(:,i)-q(:,1),p(:,i)/m(i),h);
    q(:,i) = a + q(:,1);
    p(:,i) = m(i)*b;
    ppart = ppart + p(:,i);
end
p(:,1) = ptot - ppart;




end










function [Q,P] = mapbj(Q,P,h,m)
global GNEWT MSCRIPT
n = size(m,2);
[q] = pos2cart(Q,m);
% Calculate \dot{p} (in Cartesian)
pdot = zeros(3,n);
for i=1:n
    for j=i+1:n
        xij = q(:,i)-q(:,j);
        rij = sqrt(sum(xij.^2));
        fac = GNEWT*m(i)*m(j)/rij^3*xij; 
        pdot(:,i) = pdot(:,i) - fac;
        pdot(:,j) = pdot(:,j) + fac;
    end
end
[Pdot] = mom2jac(pdot,m);

M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
end
for i=2:n
    Qmag = sqrt(sum(Q(:,i).^2));
    Pdot2 = GNEWT*mp(i)*Ms(i)/Qmag^3*Q(:,i);
    P(:,i) = P(:,i) + h*(Pdot(:,i) + Pdot2);
end




end



function [Q,V] = mapbjv(Q,V,h,m)
global GNEWT MSCRIPT
n = size(m,2);
[q] = pos2cart(Q,m);
% Calculate \dot{p} (in Cartesian)
vdot = zeros(3,n);
for i=1:n
    for j=i+1:n
        xij = q(:,i)-q(:,j);
        rij = sqrt(sum(xij.^2));
        fac = GNEWT/rij^3*xij; 
        vdot(:,i) = vdot(:,i) - fac*m(j);
        vdot(:,j) = vdot(:,j) + fac*m(i);
    end
end
[Vdot] = mom2jacv(vdot,m);

M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
end
for i=2:n
    Qmag = sqrt(sum(Q(:,i).^2));
    Vdot2 = GNEWT*Ms(i)/Qmag^3*Q(:,i);
    V(:,i) = V(:,i) + h*(Vdot(:,i) + Vdot2);
end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Q,V] = mapbjvb(Q,V,h,m,ik,jk)
global GNEWT MSCRIPT
n = size(m,2);
[q] = pos2cart(Q,m);
% Calculate \dot{p} (in Cartesian)
vdot = zeros(3,n);
for i=1:n
    for j=i+1:n
        if(i == ik && j == jk) %skip close encounter pair
            continue
        end
        xij = q(:,i)-q(:,j);
        rij = sqrt(sum(xij.^2));
        fac = GNEWT/rij^3*xij; 
        vdot(:,i) = vdot(:,i) - fac*m(j);
        vdot(:,j) = vdot(:,j) + fac*m(i);
    end
end
[Vdot] = mom2jacv(vdot,m);

M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
end
for i=2:n
    Qmag = sqrt(sum(Q(:,i).^2));
    Vdot2 = GNEWT*Ms(i)/Qmag^3*Q(:,i);
    V(:,i) = V(:,i) + h*(Vdot(:,i) + Vdot2);
%     test = (Vdot(:,i) + Vdot2);
%     test
end
% test = (Vdot + Vdot2);
% test




end







% function [Q,V] = mapbjvbp(Q,V,h,m,ik,jk)
% global GNEWT MSCRIPT
% n = size(m,2);
% [q] = pos2cart(Q,m);
% % Calculate \dot{p} (in Cartesian)
% vdot = zeros(3,n);
% for i=1:n
%     for j=i+1:n
%         if(i == ik && j == jk) %skip close encounter pair
%             continue
%         end
%         xij = q(:,i)-q(:,j);
%         rij = sqrt(sum(xij.^2));
%         fac = GNEWT/rij^3*xij; 
%         vdot(:,i) = vdot(:,i) - fac*m(j);
%         vdot(:,j) = vdot(:,j) + fac*m(i);
%     end
% end
% [Vdot] = mom2jacv(vdot,m);
% 
% M(1) = m(1);
% for i=2:n
%     M(i) = M(i-1) + m(i);
%     mp(i) = m(i)*M(i-1)/M(i);
%     if(MSCRIPT == 1)
%         Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
%     else
%         Ms(i) = M(i);
%     end
% end
% for i=2:n
%     Qmag = sqrt(sum(Q(:,i).^2));
%     Vdot2 = GNEWT*Ms(i)/Qmag^3*Q(:,i);
%     V(:,i) = V(:,i) + h*(Vdot(:,i) + Vdot2);
% %     test = (Vdot(:,i) + Vdot2);
% %     test
% end
% % test = (Vdot + Vdot2);
% % test
% 
% 
% 
% 
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin switching method in Jacobi coordinates

function [Q,V] = mapjaca(Q,V,h,m) %do not have a planetary close encounter
    [Q,V] = mapbjv(Q,V,h/2,m);
 
    [Q,V] = mapajv(Q,V,h,m);

    [Q,V] = mapbjv(Q,V,h/2,m);
%            Q
%     blergh
    

end







function [Q,V] = mapjacb(Q,V,h,m,ik,jk) %does have a planetary close encounter
    
%     [Q,V] = mapbjvb(Q,V,h/2,m,ik,jk);
%     [Q,V] = mapajvb(Q,V,h,m,ik,jk);
%     [Q,V] = mapbjvb(Q,V,h/2,m,ik,jk);

    [Q,V] = mapajvbp(Q,V,h,m,ik,jk);

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,P] = mapaj(Q,P,h,m)
global GNEWT MSCRIPT
n = size(m,2);
M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
    gm = GNEWT*Ms(i);
    a = Q(:,i);
    b = P(:,i)/mp(i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    P(:,i) = b(:)*mp(i);
    Q(:,i) = a(:);
end




end


function [Q,V] = mapajv(Q,V,h,m)
global GNEWT MSCRIPT
n = size(m,2);
M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
    gm = GNEWT*Ms(i);
    a = Q(:,i);
    b = V(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    V(:,i) = b(:);
    Q(:,i) = a(:);
end




end

function [Q,V] = mapajvb(Q,V,h,m,ik,jk)
global GNEWT MSCRIPT TOL
n = size(m,2);


flg = 0;
M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
end

for i=2:n
    if((i == ik || i == jk) && (flg == 0)) %%adjacent, j > i
%         ik
        z0 = [Q(:,ik); Q(:,jk); V(:,ik); V(:,jk)];
        tspan = [0 h];
        [z, info] = BulirschStoer(@(t,z) odesolvervjac(t,z,m,M,Ms,mp,ik,jk),tspan,z0,TOL);
         Q(:,ik) = z(1:3,2);
         Q(:,jk) = z(4:6,2);
         V(:,ik) = z(7:9,2);
         V(:,jk) = z(10:12,2);
        flg = 1;
        continue;
    end
    if((i == ik || i == jk))
        continue
    end
    gm = GNEWT*Ms(i); %kepler solver for non pairs...
    a = Q(:,i);
    b = V(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    V(:,i) = b(:);
    Q(:,i) = a(:);
end




end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Q,V] = mapajvbp(Q,V,h,m,ik,jk)
global GNEWT MSCRIPT TOL

n = size(m,2);
M(1) = m(1);
for i=2:n
    M(i) = M(i-1) + m(i);
    mp(i) = m(i)*M(i-1)/M(i);
    if(MSCRIPT == 1)
        Ms(i) = m(1)*M(i)/M(i-1); %Wisdom-Holman choice
    else
        Ms(i) = M(i);
    end
end
% ik
% jk

z0 = [Q(:, 1); Q(:,2); Q(:,3); V(:,1); V(:,2); V(:,3)];
tspan = [0 h];
[z, info] = BulirschStoer(@(t,z) odesolvervjacp(t,z,m,M,Ms,mp,ik,jk),tspan,z0,TOL);
Q(:,1) = z(1:3,2);
Q(:,2) = z(4:6,2);
Q(:,3) = z(7:9,2);
V(:,1) = z(10:12,2);
V(:,2) = z(13:15,2);
V(:,3) = z(16:18,2);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,P] = corrector(Q,P,h,m)
% third-order corrector try
c = 2;
alpha = 10^(1/2);
a = alpha*[3/10 1/5];
b = alpha*[-1/72 1/24];

for i=c:-1:1
    [Q,P] = zmap(Q,P,a(i)*h,b(i)*h,m);
end


end

function [Q,P] = invcorrector(Q,P,h,m)
% third-order corrector try
c = 2;
alpha = 10^(1/2);
a = alpha*[3/10 1/5];
b = alpha *[-1/72 1/24];


for i=1:c
    [Q,P] = zmap(Q,P,a(i)*h,-b(i)*h,m);
end


end

function [Q,P] = zmap(Q,P,h1,h2,m)
[Q,P] = xmap(Q,P,-h1,-h2,m);
[Q,P] = xmap(Q,P,h1,h2,m);


end


function [Q,P] = xmap(Q,P,h1,h2,m)
[Q,P] = mapbd(Q,P,-h1,m);
[Q,P] = mapad(Q,P,h2,m);
[Q,P] = mapbd(Q,P,h1,m);


end


function [Q,P] = mapad(Q,P,h,m)
global GNEWT


n = size(m,2);
for i=2:n
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = P(:,i)/m(i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    P(:,i) = b(:)*m(i);
    Q(:,i) = a(:);
end



end





function [Q,v] = mapadv(Q,v,h,m)
global GNEWT

n = size(m,2);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end



n = size(m,2);
for i=2:n
    gm = GNEWT*m(1);
    a = Q(:,i);
    b = v(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    v(:,i) = b(:);
    Q(:,i) = a(:);
end

% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));

end



function [Psum0] = calcm(v,m)
n = size(m,2);

% total momentum initially-- conserved
Psum0 = zeros(3,1);
for i = 1:n
    Psum0(:) = Psum0(:) + m(i)*v(:,i); 
end

end


function [xo,vo] = adjustcm(x,v,m)
n = size(m,2);

ms = sum(m);
% total momentum initially-- conserved
vcm = zeros(3,1);
xcm = zeros(3,1);
for i = 1:n
    vcm(:) = vcm(:) + m(i)*v(:,i); 
    xcm(:) = xcm(:) + m(i)*x(:,i);
end
vcm = vcm/ms;
xcm = xcm/ms;

for i=1:n
    xo(:,i) = x(:,i) - xcm(:);
    vo(:,i) = v(:,i) - vcm(:);
end


end


function [v] = adjustSun(Q,v,m,Psum0)
n = size(m,2);

% Calculate solar velocity, (does not assume CoM stationary)
Psum = zeros(3,1);
for i = 2:n
    Psum(:) = Psum(:) + m(i)*v(:,i); 
end
% Psum0
v(:,1) = 1/m(1)*(Psum0(:)-Psum(:));

end



function [Q,v] = Keppair(Q,v,h,m,i,j)
global GNEWT


gm = GNEWT*m(1);
a = Q(:,i);
b = v(:,i);
[a,b] = kepler_stepxv(gm,a,b,h);
v(:,i) = b(:);
Q(:,i) = a(:);

a = Q(:,j);
b = v(:,j);
[a,b] = kepler_stepxv(gm,a,b,h);
v(:,j) = b(:);
Q(:,j) = a(:);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q,v,kepdone] = Kepfull(Q,v,h,m,lev,levc,kepdone,levmax)
global GNEWT

n = size(m,2);
intflg = zeros(1,n);
gm = GNEWT*m(1);
for i=1:n
    for j=i+1:n
        if(levc(i,j) == lev)
            if(lev == levmax) %yes, always do Kepler solvers, even if they were done before at same level
                if(intflg(i) == 0)
%                     display('Kep')
                    a = Q(:,i);
                    b = v(:,i);
                    [a,b] = kepler_stepxv(gm,a,b,h);
                    v(:,i) = b(:);
                    Q(:,i) = a(:);
                    kepdone(lev,i) = 1;
                    intflg(i) = 1; %don't loop over same particle twice
%                     i
%                     lev
%                     display('max')
                end

                if(intflg(j) == 0)
%                     display('Kep')
                    a = Q(:,j);
                    b = v(:,j);
                    [a,b] = kepler_stepxv(gm,a,b,h);
                    v(:,j) = b(:);
                    Q(:,j) = a(:);
                    kepdone(lev,j) = 1;
                    intflg(j) = 1;
%                     j
%                     lev
%                     display('max')
                end
            else %Kepler solvers only if higher levels have not done a Kepler solver, even if current level already done
                flagi = 0;
                flagj = 0;
                for k=lev+1:levmax
                    if(kepdone(k,i) == 1)
                        flagi = 1;
                    end
                    if(kepdone(k,j) == 1)
                        flagj = 1;
                    end
                end
                if(flagi == 0)
                    if(intflg(i) == 0)
%                         display('Kep')
                        a = Q(:,i);
                        b = v(:,i);
                        [a,b] = kepler_stepxv(gm,a,b,h);
                        v(:,i) = b(:);
                        Q(:,i) = a(:);
                        kepdone(lev,i) = 1;
                        intflg(i) = 1;
%                         i
%                         lev
%                         display('mid')
                    end
                end
                if(flagj == 0)
                    if(intflg(j) == 0)
%                         display('Kep')
                        a = Q(:,j);
                        b = v(:,j);
                        [a,b] = kepler_stepxv(gm,a,b,h);
                        v(:,j) = b(:);
                        Q(:,j) = a(:);
                        kepdone(lev,j) = 1;
                        intflg(j) = 1;
%                         j
%                         lev
%                         display('mid')
                    end
                end
            end
        end
    end
end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v] = Keppairsmooth(Q,v,h,m,i,j)
global GNEWT


gm = GNEWT*m(1);
a = Q(:,i);
b = v(:,i);
[a,b] = kepler_stepxv(gm,a,b,h);
v(:,i) = b(:);
Q(:,i) = a(:);

a = Q(:,j);
b = v(:,j);
[a,b] = kepler_stepxv(gm,a,b,h);
v(:,j) = b(:);
Q(:,j) = a(:);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,v] = Keppairsmoothfull(Q,v,h,m)
global GNEWT

n = size(m,2);
gm = GNEWT*m(1);

for i=2:n
    a = Q(:,i);
    b = v(:,i);
    [a,b] = kepler_stepxv(gm,a,b,h);
    v(:,i) = b(:);
    Q(:,i) = a(:);
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,P] = mapbl(Q,P,h,m)
global GNEWT
sp = zeros(3,1);
n = size(m,2);

for i=1:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = norm(qvec);
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end

end


function [Q,P] = mapal(Q,P,h,m)
global GNEWT
n = size(m,2);

for i=1:n
    Q(:,i) = Q(:,i) + h*P(:,i)/m(i);
end

end



function [Q,v] = mapbdv(Q,v,h,m)
global GNEWT
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + m(i)*v(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        fac = GNEWT/qmod^3*qvec;
        v(:,i) = v(:,i) - h * fac*m(j);
        v(:,j) = v(:,j) + h * fac*m(i);
    end
end



end




function [Q] = mapSun(Q,v,h,m)
global GNEWT
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + m(i)*v(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end


end


function [v] = interpair(Q,v,h,m,ik,jk)
global GNEWT
n = size(m,2);

% pair interaction
qvec = Q(:,ik) - Q(:,jk);
qmod = sqrt(sum((Q(:,ik)-Q(:,jk)).^2));
fac = GNEWT/qmod^3*qvec;
v(:,ik) = v(:,ik) - h * fac*m(jk);
v(:,jk) = v(:,jk) + h * fac*m(ik);

end



function [v] = interfull(Q,v,h,m,lev,levc)
global GNEWT
n = size(m,2);
for i=2:n
    for j=i+1:n
        if(levc(i,j) ==  lev) %check if pair at this level
%             display('interf')
%             i
%             j
%             lev
            qvec = Q(:,i) - Q(:,j);
            qmod = sqrt(sum(qvec.^2));
            fac = GNEWT/qmod^3*qvec;
            v(:,i) = v(:,i) - h * fac*m(j);
            v(:,j) = v(:,j) + h * fac*m(i);
        end
    end
end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [v] = interpairsmooth(Q,v,h,m,rlev,levc,ik,jk)
global GNEWT
n = size(m,2);
if(levc == 1)
    r1 = rlev(1);
    r2 = rlev(2);
    rvec = Q(:,ik)-Q(:,jk);
    r = sqrt(sum(rvec.^2));
    if(r > r1)
        fac = GNEWT/r^3*rvec;
    elseif(r2 <= r && r < r1)
        x = (r1-r)/(r1-r2);
        fct = 2*x^3 - 3*x^2 + 1;
        fac = GNEWT/r^3*rvec*fct;
    else
        return %first level no contribution
    end
else
    r1 = rlev(levc-1);
    r2 = rlev(levc);
    r3 = rlev(levc+1);
    rvec = Q(:,ik)-Q(:,jk);
    r = sqrt(sum(rvec.^2));
    if(r2 <= r && r < r1)
        x = (r1-r)/(r1-r2);
        fct = 2*x^3 - 3*x^2 + 1;
        fac = GNEWT/r^3*rvec*(1-fct);
    elseif(r3 <= r && r < r2)
        x = (r2-r)/(r2-r3);
        fct = 2*x^3 - 3*x^2 + 1;
        fac = GNEWT/r^3*rvec*fct;
    else
        return
    end
end

v(:,ik) = v(:,ik) - h * fac*m(jk);
v(:,jk) = v(:,jk) + h * fac*m(ik);



end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [v] = interpairsmoothfull(Q,v,h,m,rlev,levc)
global GNEWT
n = size(m,2);

for i=2:n
    for j=i+1:n
        if(levc == 1)
            r1 = rlev(1);
            r2 = rlev(2);
            rvec = Q(:,i)-Q(:,j);
            r = sqrt(sum(rvec.^2));
            if(r > r1)
                fac = GNEWT/r^3*rvec;
            elseif(r2 <= r && r < r1)
                x = (r1-r)/(r1-r2);
                fct = 2*x^3 - 3*x^2 + 1;
                fac = GNEWT/r^3*rvec*fct;
            else
                return %first level no contribution
            end
        else
            ind = levc;
            r1 = rlev(ind-1);
            r2 = rlev(ind);
            r3 = rlev(ind+1);
            rvec = Q(:,i)-Q(:,j);
            r = sqrt(sum(rvec.^2));
            if(r2 <= r && r < r1)
                x = (r1-r)/(r1-r2);
                fct = 2*x^3 - 3*x^2 + 1;
                fac = GNEWT/r^3*rvec*(1-fct);
            elseif(r3 <= r && r < r2)
                x = (r2-r)/(r2-r3);
                fct = 2*x^3 - 3*x^2 + 1;
                fac = GNEWT/r^3*rvec*fct;
            else
                return
            end
        end
        v(:,i) = v(:,i) - h * fac*m(j);
        v(:,j) = v(:,j) + h * fac*m(i);
    end
end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,v] = driftop(Q,v,m,lev,rlev,hlev,hsub,levmax,jk,kk)
if(lev < levmax)
    for i=1:hsub
        [v] = interpairsmooth(Q,v,hlev(lev)/2,m,rlev,lev,jk,kk);
        [Q,v] = driftop(Q,v,m,lev+1,rlev,hlev,hsub,levmax,jk,kk);
        [v] = interpairsmooth(Q,v,hlev(lev)/2,m,rlev,lev,jk,kk);
    end
elseif(lev == levmax)
    [Q,v] = Keppairsmooth(Q,v,hlev(lev-1),m,jk,kk); %decrease by one because added one before
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,v,kepdone] = driftopfull(Q,v,m,lev,rlev,hlev,hsub,levmax,levc,kepdone)
if(lev < levmax)
    for i=1:hsub
        [v] = interpairsmoothfull(Q,v,hlev(lev)/2,m,rlev,lev);
        [Q,v,kepdone] = driftopfull(Q,v,m,lev+1,rlev,hlev,hsub,levmax,levc,kepdone);
        [v] = interpairsmoothfull(Q,v,hlev(lev)/2,m,rlev,lev);
    end
elseif(lev == levmax)
    [Q,v] = Keppairsmoothfull(Q,v,hlev(lev-1),m); %decrease by one because added one before
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,v,kepdone] = driftopdiscrete(Q,v,m,lev,hlev,hsub,levmax,levc,kepdone)
if(lev <= levmax)
    for i=1:hsub
        [v] = interfull(Q,v,hlev(lev)/2,m,lev,levc);%one interaction at each level
        if(lev < levmax)
            [Q,v,kepdone] = driftopdiscrete(Q,v,m,lev+1,hlev,hsub,levmax,levc,kepdone);
        end
        [Q,v,kepdone] = Kepfull(Q,v,hlev(lev),m,lev,levc,kepdone,levmax); %one Kepler solver at each level; Kepler operator commutes with H_I and other H_K, because different particles
        [v] = interfull(Q,v,hlev(lev)/2,m,lev,levc);
    end
elseif(lev == levmax) 
    return
end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Q,v,kepdone,levsub] = driftopdiscrete_m(Q,v,m,lev,hlev,hsub,levmax,levc,kepdone,levsub,rlev)
if(lev <= levmax)
    for i=1:hsub
        [v] = interfull(Q,v,hlev(lev)/2,m,lev,levc);%one interaction at each level
        if(lev < levmax)
            [Q,v,kepdone,levsub] = driftopdiscrete_m(Q,v,m,lev+1,hlev,hsub,levmax,levc,kepdone,levsub,rlev);
        end
        [Q,v,kepdone] = Kepfull(Q,v,hlev(lev),m,lev,levc,kepdone,levmax); %one Kepler solver at each level; Kepler operator commutes with H_I and other H_K, because different particles
        [v] = interfull(Q,v,hlev(lev)/2,m,lev,levc); %looping over levels in each subfunction
%         return what pairs were treated at this level, calculate new
%         levels for those pairs.  Array levp(i,j,k).  But actually, just
%         need to store max level, not even that k value..?
% i
% hsub
        [levsub] = calclevvsub(Q,v,m,rlev,hlev(1),lev,levc,levsub);
    end
elseif(lev == levmax) 
    return
end




end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [Q,P] = mapbd(Q,P,h,m)
global GNEWT
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
for i=2:n
    Q(:,i) = Q(:,i) + h/(m(1))*sp;
end
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end



end






function [Q,P] = mapbhel(Q,P,h,m)
global GNEWT
sp = zeros(3,1);
n = size(m,2);
for i = 2:n
    sp = sp + P(:,i);
end
Q(:,1) = Q(:,1) + h/m(1)*(P(:,1)-sp);

for i=2:n
    Q(:,i) = Q(:,i) + h/m(1)*(sp - P(:,1));
end
for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end




end

function [Q,P] = mapadt(Q,P,h,m)
global GNEWT

n = size(m,2);
for i=2:n
    gm = GNEWT*(m(1) + m(i));
    a = Q(:,i);
    mu = m(i)*m(1)/(m(i) + m(1));
    b = P(:,i)/mu;
    [a,b] = kepler_stepxv(gm,a,b,h);
    P(:,i) = b(:)*mu;
    Q(:,i) = a(:);
end

end



function [Q,P] = mapbdt(Q,P,h,m)
global GNEWT
n = size(m,2);

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end

sp = zeros(3,1);
for i = 2:n
    sp = sp + P(:,i);
end

for i=2:n
    spk = sp - P(:,i);
    Q(:,i) = Q(:,i) + h/(m(1))*spk;
end


end


function [Q,P] = mapbdttrans(Q,P,h,m)
global GNEWT
n = size(m,2);

sp = zeros(3,1);
for i = 2:n
    sp = sp + P(:,i);
end

for i=2:n
    spk = sp - P(:,i);
    Q(:,i) = Q(:,i) + h/(m(1))*spk;
end

for i=2:n
    for j=i+1:n
        qvec = Q(:,i) - Q(:,j);
        qmod = sqrt(sum((Q(:,i)-Q(:,j)).^2));
        fac = GNEWT*m(i)*m(j)/qmod^3*qvec;
        P(:,i) = P(:,i) - h * fac;
        P(:,j) = P(:,j) + h * fac;
    end
end




end





function [l1] = l1norm(x1,p1,x2,p2,x0,p0)
delx = sum(sum(abs(x1-x2),2),1);
xn = sum(sum(abs(x0),2),1);
delp = sum(sum(abs(p1-p2),2),1);
pn = sum(sum(abs(p0),2),1);
l1 = delx/xn + delp/pn;

end








function [x] = driftij(x,v,i,j,h,m)
    x(:,i) = x(:,i) + h*v(:,i);
    x(:,j) = x(:,j) + h*v(:,j);
end



  
function [x,v] = centerm(m,x,v,delx,delv,i,j,h)
    mij = m(i)+m(j);
    vcm = (m(i)*v(:,i) + m(j)*v(:,j))/mij;
    x(:,i) = x(:,i) + m(j)/mij*delx + h*vcm;
    x(:,j) = x(:,j) - m(i)/mij*delx + h*vcm;
    v(:,i) = v(:,i) + m(j)/mij*delv;
    v(:,j) = v(:,j) - m(i)/mij*delv;
end









function [e] = finde(x,p,m,i,j)
    gm = m(i)+m(j);
    xr = x(:,i)-x(:,j);
    vr = p(:,i)/m(i)-p(:,j)/m(j);
    r0 = sqrt(sum(xr.^2));
    v0s = sum(vr.^2);
    u = sum(xr.*vr);
    alpha=2.0*gm/r0-v0s;
    zetap = gm - alpha*r0;
    %determine initial guesses*/
    if(alpha <= 0)
        %hyperbolic motion
        a=gm/alpha;
        en=sqrt(-gm/(a^3));
        ch=1-r0/a;
        sh=u/sqrt(-a*gm);
        e=sqrt(ch*ch-sh*sh);
    else
        %elliptic motion
        a=gm/alpha;
        en=sqrt(gm/a^3);
        ec = 1-r0/a;
        es=u/(en*a^2);
        e=sqrt(ec^2+es^2);
    end
end








   





function [fij] = ftwo(x,m,i,j)
rij = x(:,i)-x(:,j);
r2 = sum(rij.^2);
r3 = r2.^1.5;
fij = - rij*(m(i)*m(j)./r3)';
%}\medskip
end





%Normal DKD kick step
function [p] = kick(x,p,h,m)
n = size(m,2);
for i=1:n
    j = [1:i-1,i+1:n];
    rij = x(:,i)*ones(1,n-1) - x(:,j);
    r2 = sum(rij.^2);
    r3 = r2.^1.5;
    force = - rij*(m(i)*m(j)./r3)';
    p(:,i) = p(:,i) + h*force;
end
end

function [a] = accel(x,m)
global GNEWT
n = size(m,2);
d = size(x,1);
a = zeros(d,n);
for i=1:n
    j = [1:i-1,i+1:n];
    rij = x(:,i)*ones(1,n-1) - x(:,j);
    r2 = sum(rij.^2);
    r3 = r2.^1.5;
    a(:,i) = - rij*(GNEWT*m(j)./r3)';
end
end


function [v] = kickpair(x,v,h,m,pair)
global GNEWT
n = size(m,2);
for i=1:n
    for j=i+1:n
        if(pair(i,j) == 1)
            rij = x(:,i) - x(:,j);
            r2 = sum(rij.^2);
            r3 = r2.^(1.5);
            fac = h*GNEWT*rij./r3;
            v(:,i) = v(:,i) - m(j)*fac;
            v(:,j) = v(:,j) + m(i)*fac;
        end
    end
end
end


function [v] = kickwhd(x,v,h,m)
global GNEWT
n = size(m,2);
% Kick only non-Solar particles.  Assume Sun has i=0.
% Particles have heliocentric coordinates
for i=2:n
    for j=i+1:n
        rij = x(:,i) - x(:,j);
        r2 = sum(rij.^2);
        r3 = r2.^(1.5);
        fac = h*GNEWT*rij./r3;
        v(:,i) = v(:,i) - m(j)*fac;
        v(:,j) = v(:,j) + m(i)*fac;
    end
end
end




function [v] = kicksun(x,v,h,m)
global GNEWT
n = size(m,2);
% Kick only the Sun.  Assume Sun has i=1.
for i=2:n
    rij = x(:,i) - x(:,1);
    r2 = sum(rij.^2);
    r3 = r2.^(1.5);
    fac = h*GNEWT*rij./r3;
    v(:,1) = v(:,1) + m(i)*fac;
end
end





























% % % Begin Kepler solver
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 


function [X,S2,C2,SUCCESS] = solve_universal_newton(kc,r0,beta,b,eta,zeta,h,X)
  xnew = X;
%   Initialize in case of no success
  S2 = 0;
  C2 = 0;
  X = 0;

  count = 0;
  flag = 0;
  while(flag == 0)
    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    cc = eta*g1 + zeta*g2;
    xnew = (h + (x*cc - (eta*g2 + zeta*g3)))/(r0 + cc);
    count = count + 1;
    if(count > 10)
        SUCCESS = 0;
        return
    end
    if(abs((x-xnew)/xnew) <= 1.e-8)
        flag = 1;
    end
  end
%    compute the output 
  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);
  
  X = x;
  S2 = s2;
  C2 = c2;
  SUCCESS = 1;
end

function [X,S2,C2,SUCCESS] = solve_universal_laguerre(kc,r0,beta,b,eta,zeta,h,X)
    xnew = X;
    X = 0;
    S2 = 0;
    C2 = 0;

  count = 0;
  flag = 0;
  while(flag == 0)
    c5=5.0;
    c16 = 16.0;
    c20 = 20.0;

    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    fp = r0 + eta*g1 + zeta*g2;
    g0 = 1.0 - beta*g2;
    fpp = eta*g0 + zeta*g1;
    dx = -c5*f/(fp + sqrt(abs(c16*fp*fp - c20*f*fpp)));
    xnew = x + dx;
    count = count + 1;
    if(count > 10) 
        SUCCESS = 0;
        return
    end
    if(abs(dx) <= 2.e-7*abs(xnew))
        flag = 1;
    end
  end

%   /* refine with a last newton */
  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);
  g1 = 2.0*s2*c2/b;
  g2 = 2.0*s2*s2/beta;
  g3 = (x - g1)/beta;
  cc = eta*g1 + zeta*g2;
  xnew = (h + (x*cc - (eta*g2 + zeta*g3)))/(r0 + cc);

  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);

  X = x;
  S2 = s2;
  C2 = c2;

  SUCCESS = 1;
end

function [X,S2,C2,SUCCESS] = solve_universal_bisection(kc,r0,beta,b,eta,zeta,h,X)

  xnew = X;
  err = 1.e-9*abs(xnew);

%   /* bisection limits due to Rein et al. */
  invperiod = b*beta/(2.*UNIVERSAL_PI*kc);
  X_per_period = 2.*UNIVERSAL_PI/b;
  X_min = X_per_period * floor(h*invperiod);
  X_max = X_min + X_per_period;
  xnew = (X_max + X_min)/2.;

  count = 0;
  flag = 0;
  while(flag == 0)
    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    if (f>=0.)
      X_max = x;
    else
      X_min = x;
    end
    xnew = (X_max + X_min)/2.;
    count = count + 1;
    if(count > 10)
        SUCCESS = 0;
        return
    end
    if(abs(x - xnew) <= err)
        flag = 1;
    end

  end
  x = xnew;
  
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);

  X = x;
  S2 = s2;
  C2 = c2;
  SUCCESS = 1;
end

function [sgn] = sign(x)
if (x > 0)
    sgn = 1.0;
elseif(x < 0)
    sgn = -1.0;
else
    sgn = 0;
end
end


function [x1] = cubic1(a,b,c)
  
  Q = (a*a - 3.0*b)/9.0;
  R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  if(R*R < Q*Q*Q) 
    theta = acos(R/sqrt(Q*Q*Q));
    x1 = 0;
    return;
  else
      sgn = sign(R);
    A = -sgn*(abs(R) + sqrt(R*R - Q*Q*Q))^(1/3);
    if(A == 0.0) 
      B = 0.0;
    else 
      B = Q/A;
    end
    x1 = (A + B) - a/3.0;
  end
end


function [X,S2,C2,SUCCESS] = solve_universal_parabolic(kc,r0,beta,b,eta,zeta,h,X)
  b = 0.0;
  s2 = 0.0;
  c2 = 1.0;
%   
%     g1 = x;
%     g2 = x*x/2.0;
%     g3 = x*x*x/6.0;
%     g = eta*g1 + zeta*g2;
%   
%      f = r0*x + eta*g2 + zeta*g3 - h; 
%      this is just a cubic equation 
  x = cubic1(3.0*eta/zeta, 6.0*r0/zeta, -6.0*h/zeta);

  s2 = 0.0;
  c2 = 1.0;
%   /* g1 = 2.0*s2*c2/b; */
%   /* g2 = 2.0*s2*s2/beta; */

  X = x;
  S2 = s2;
  C2 = c2;
  
  SUCCESS = 1;
end

function [X,S2,C2,SUCCESS] = solve_universal_hyperbolic_newton(kc,r0,minus_beta,b,eta,zeta,h,X)

  xnew = X;
  count = 0;
  flag = 0;
  while(flag == 0)
    x = xnew;
    arg = b*x/2.0;
    if(abs(arg)>200.0) 
        SUCCESS = 0;
        return
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    g = eta*g1 + zeta*g2;
    xnew = (x*g - eta*g2 - zeta*g3 + h)/(r0 + g);
    count = count + 1;
    if(count > 10) 
        SUCCESS = 0;
        return;
    end
    if(abs(x - xnew) <= 1.e-9*abs(xnew))
        flag = 1;
    end
  end
  
  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  X = x;
  S2 = s2;
  C2 = c2;

  SUCCESS = 1;
end

function [X,S2,C2,SUCCESS] = solve_universal_hyperbolic_laguerre(kc,r0,beta,b,eta,zeta,h,X)

  xnew = X;
  count = 0;
  flag = 0;
  while(flag == 0)
    c5=5.0; c16 = 16.0; c20 = 20.0;

    x = xnew;
    arg = b*x/2.0;
    if(abs(arg)>50.0) 
        SUCCESS = 0;
        return;
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    fp = r0 + eta*g1 + zeta*g2;
    g0 = 1.0 + minus_beta*g2;
    fpp = eta*g0 + zeta*g1;
    den = (fp + sqrt(abs(c16*fp*fp - c20*f*fpp)));
    if(den == 0.0) 
        SUCCESS = 0;
        return;
    end
    dx = -c5*f/den;
    xnew = x + dx;
    count = count + 1;
    if(count > 20) 
        SUCCESS = 0;
        return;
    end
    if(abs(x - xnew) <= 1.e-9*abs(xnew))
        flag = 1;
    end
  end
  
%   /* refine with a step of Newton */

    x = xnew;
    arg = b*x/2.0;
    if(abs(arg)>200.0) 
        SUCCESS = 0;
        return;
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    g = eta*g1 + zeta*g2;
    xnew = (x*g - eta*g2 - zeta*g3 + h)/(r0 + g);

  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  X = x;
  S2 = s2;
  C2 = c2;

  SUCCESS = 1;
end

function [X,S2,C2,SUCCESS] = solve_universal_hyperbolic_bisection(kc,r0,beta,b,eta,zeta,h,X)
  xnew = X;

  err = 1.e-10*abs(xnew);

  X_min = 0.5*xnew;
  X_max = 10.0*xnew;

    x = X_min;
    arg = b*x/2.0;
    if(abs(arg)>200.0) 
        SUCCESS = 0;
        return
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    fmin = r0*x + eta*g2 + zeta*g3 - h;

    x = X_max;
    arg = b*x/2.0;
    if(abs(arg)>200.0) 
      x = 200.0/(b/2.0);
      arg = 200.0;
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    fmax = r0*x + eta*g2 + zeta*g3 - h;

    if(fmin*fmax > 0.0) 
        SUCCESS = 0;
        return
    end

  count = 0;
  flag = 0;
  
  while(flag == 0)
    x = xnew;
    arg = b*x/2.0;
    if(abs(arg)>200.0) 
        SUCCESS = 0;
        return;
    end
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    if (f>=0.)
      X_max = x;
    else
      X_min = x;
    end
    xnew = (X_max + X_min)/2.;
    count = count + 1;
    if(count > 100) 
        SUCCESS = 0;
        return;
    end
    if(abs(x - xnew) <= err)
        return
    end

  end
  
  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  X = x;
  S2 = s2;
  C2 = c2;

  SUCCESS = 1;
end

function [s] = kepler_step(kc,dt,s0)
  x = s0(1);
  y = s0(2);
  z = s0(3);
  xd = s0(4);
  yd = s0(5);
  zd = s0(6);
  
  r0 = sqrt(x*x + y*y + z*z);
  v2 = xd*xd + yd*yd + zd*zd;
  eta = x*xd + y*yd + z*zd;
  beta = 2*kc/r0 - v2;
  zeta = kc - beta*r0;
  b = sqrt(abs(beta));
  [s] = kepler_step_depth(kc, dt, beta, b, s0, 0, r0, v2, eta, zeta);
end

function [s] = kepler_step_depth(kc,dt,beta,b,s0,depth,r0,v2,eta,zeta)

  if(depth > 30)
    printf('kepler depth exceeded\n');
     error('Program exit')
  end

  [s,SUCCESS] = kepler_step_internal(kc, dt, beta, b, s0, r0, v2, eta, zeta);
  if(SUCCESS == 0) 
    [ss] = kepler_step_depth(kc, dt/4.0, beta, b, s0, depth+1, r0, v2, eta, zeta);

      x = ss(1);
      y = ss(2);
      z = ss(3);
      xd = ss(4);
      yd = ss(5);
      zd = ss(6);
  
      r0 = sqrt(x*x + y*y + z*z);
      v2 = xd*xd + yd*yd + zd*zd;
      eta = x*xd + y*yd + z*zd;
      beta = 2*kc/r0 - v2;
      zeta = kc - beta*r0;
      b = sqrt(abs(beta));
      
      [s] = kepler_step_depth(kc, dt/4.0, beta, b, ss, depth+1, r0, v2, eta, zeta);

      x = s(1);
      y = s(2);
      z = s(3);
      xd = s(4);
      yd = s(5);
      zd = s(6);
  
      r0 = sqrt(x*x + y*y + z*z);
      v2 = xd*xd + yd*yd + zd*zd;
      eta = x*xd + y*yd + z*zd;
      beta = 2*kc/r0 - v2;
      zeta = kc - beta*r0;
      b = sqrt(abs(beta));
      
      [ss] = kepler_step_depth(kc, dt/4.0, beta, b, s, depth+1, r0, v2, eta, zeta);

      x = ss(1);
      y = ss(2);
      z = ss(3);
      xd = ss(4);
      yd = ss(5);
      zd = ss(6);
  
      r0 = sqrt(x*x + y*y + z*z);
      v2 = xd*xd + yd*yd + zd*zd;
      eta = x*xd + y*yd + z*zd;
      beta = 2*kc/r0 - v2;
      zeta = kc - beta*r0;
      b = sqrt(abs(beta));
      
      [s] = kepler_step_depth(kc, dt/4.0, beta, b, ss, depth+1, r0, v2, eta, zeta);

  end
end

function [s] =  new_guess(r0,eta,zeta,dt)

  if(zeta ~= 0.0) 
    s = cubic1(3.0*eta/zeta, 6.0*r0/zeta, -6.0*dt/zeta);
  elseif(eta ~= 0.0) 
    reta = r0/eta;
    disc = reta*reta + 8.0*dt/eta;
    if(disc >= 0.0) 
      s = -reta + sqrt(disc);
    else
      s = dt/r0;
    end
  else
    s = dt/r0;
  end        
end

function [s,SUCCESS] = kepler_step_internal(kc,dt,beta,b,s0,r0,v2,eta,zeta)


  if(beta < 0.0) 
    x0 = new_guess(r0, eta, zeta, dt);
    x = x0;
   [x,s2,c2,SUCCESS] = solve_universal_hyperbolic_newton(kc, r0, -beta, b, eta, zeta, dt, x);
    if(SUCCESS == 0) 
      x = x0;
      [x,s2,c2,SUCCESS] = solve_universal_hyperbolic_laguerre(kc, r0, -beta, b, eta, zeta, dt, x);
    end
% 
%     /*
%     if(flag == FAILURE) {
%       x = x0;
%       flag = solve_universal_hyperbolic_bisection(kc, r0, -beta, b, eta, zeta, dt, &x, &s2, &c2);
%     }
%       */
        if(SUCCESS == 0)
            return;
        end


    a = kc/(-beta);
    G1 = 2.0*s2*c2/b;
    c = 2.0*s2*s2;
    G2 = c/(-beta);
    ca = c*a;
    r = r0 + eta*G1 + zeta*G2;
    bsa = (a/r)*(b/r0)*2.0*s2*c2;

    elseif(beta > 0.0)
%     /* x0 = dt/r0; */

    x0 = dt/r0;
    ff = zeta*x0*x0*x0 + 3.0*eta*x0*x0;
    fp = 3.0*zeta*x0*x0 + 6.0*eta*x0 + 6.0*r0;
    x0 = x0 - ff/fp;

    x = x0;
    [x,s2,c2,SUCCESS] = solve_universal_newton(kc, r0, beta, b, eta, zeta, dt, x);
    if(SUCCESS == 0)
      x = x0;
      [x,s2,c2,SUCCESS] = solve_universal_laguerre(kc, r0, beta, b, eta, zeta, dt, x);
%     /*
%     if(flag == FAILURE) {
%       x = x0;
%       flag = solve_universal_bisection(kc, r0, beta, b, eta, zeta, dt, &x, &s2, &c2);
%     }
%       */
    end
    if(SUCCESS == 0) 
        s = 0;
      return
    end

    a = kc/beta;
    G1 = 2.0*s2*c2/b;
    c = 2.0*s2*s2;
    G2 = c/beta;
    ca = c*a;
    r = r0 + eta*G1 + zeta*G2;
    bsa = (a/r)*(b/r0)*2.0*s2*c2;

  else
    x = dt/r0;
    [x,s2,c2,SUCCESS] = solve_universal_parabolic(kc, r0, beta, b, eta, zeta, dt, x);
    if(SUCCESS == 0)
        exit
    end
    
    G1 = x;
    G2 = x*x/2.0;
    ca = kc*G2;
    r = r0 + eta*G1 + zeta*G2;
    bsa = kc*x/(r*r0);
  end

    

%     /* f = 1.0 - (ca/r0); */
    fhat = -(ca/r0);
    g = eta*G2 + r0*G1;
    fdot = -bsa;
%     /* gdot = 1.0 - (ca/r); */
    gdothat = -(ca/r);
      x = s0(1);
      y = s0(2);
      z = s0(3);
      xd = s0(4);
      yd = s0(5);
      zd = s0(6);
    
    s(1) = x + fhat*x + g*xd;
    s(2) = y + fhat*y + g*yd;
    s(3) = z + fhat*z + g*zd;
    s(4) = xd + fdot*x + gdothat*xd;
    s(5) = yd + fdot*y + gdothat*yd;
    s(6) = zd + fdot*z + gdothat*zd;

  SUCCESS = 1;
  
end

% % % End Kepler solver
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 




function dydt = kepode(t,y)
r = sqrt(y(1)^2 + y(2)^2);
dydt = [y(3); y(4); -y(1)/r; -y(2)/r];

end

function [x,v] = kepler_stepxv(GM,x,v,h)
[s] = converxv(x,v);
[s] = kepler_step(GM,h,s);
[x,v] = converst(s);

end




function [s] = converxv(x,v)
x = x';
v = v';
s = [x v];

end


function [x,v] = converst(s)
x = [s(1) s(2) s(3)];
v = [s(4) s(5) s(6)];
x = x';
v = v';

end







function [m,x0,v0,pair] = in3bodfig8bin(nbin)
% T = 6.32591398 approximately
m = [1 1 1];
xa = .97000436;
xb = -.24308753;
va = -.93240737;
vb = -.86473146;
x0 = [xa xb 0; -xa -xb 0; 0 0 0];
x0 = x0';
v0 = [-va/2 -vb/2 0; -va/2 -vb/2 0; va vb 0];
v0 = v0';
a = 1e-2;
e = .9;
[x0,v0,m] = binadd(x0,v0,m,nbin,a,e);
n = size(m,2);
for i=1:n
    for j=1:n
%         Kepler solver group
        pair(i,j) = 0;
        pair(j,i) = 0;
    end
end
end













function [x0,v0,m] = binadd(x0,v0,m,nbin,a,e)
% This version will only work for special pair ordering problem
n = size(m,2);
if(nbin > 0)
    mb = [m(1)/2 m(1)/2];
    mu = (mb(1)*mb(2))/(mb(1)+mb(2));
    mt = mb(1)+mb(2);
    for i=1:nbin        
        if(i == 1)
            xdisp = a*(1+e);
            vdisp = sqrt(mt/a*(1-e)/(1+e));
        else
            xdisp = a*(1-e);
            vdisp = sqrt(mt/a*(1+e)/(1-e));
        end
        %add particle at the end
        x0(1:3,n+i) = x0(1:3,i);
        v0(1:3,n+i) = v0(1:3,i);
        x0(1,i) = x0(1,i) - xdisp/2;
        x0(1,n+i) = x0(1,n+i) + xdisp/2;
        v0(2,i) = v0(2,i) - vdisp/2;
        v0(2,n+i) = v0(2,n+i) + vdisp/2; 
        m(i) = mb(1);
        m(n+i) = mb(2); 
    end
end
end





