function symbatwo(hin,tmax)
% tmax in period units
tmaxp = tmax*2*pi;
hp = hin*2*pi;
h = hp;

a = 1;
e = 0.999;
P = 2*pi;
q = [a*(1+e) 0];
p = [0 sqrt((1-e)/(1+e))]; %apocenter, counterclock
t = 0;
str = sprintf('data/fcons.txt');
fcons = fopen(str,'w');
[E0] = consq(q,p);
nout = 1e4;
tout = linspace(0,tmaxp,nout);
ind = 1;
rmax = sqrt(2)*a;
rmin = 0;
delr = rmax-rmin;
levtot = 30;
levv = zeros(levtot,1);

rsub = sqrt(2);
% rsub = 3;
hsub = 2; %this gives h~r^(3/2), Levison et al. (1998); must also be ratio for dividing global substep

rlev(1) = rmax;
hlev(1) = h;
steps(1) = 1;
for i=2:levtot
    rlev(i) = delr/(rsub^(i-1));
    hlev(i) = h/(hsub^(i-1)); 
    steps(i) = steps(i-1)*hsub;
end
casen = 1;
[levc] = calclev(q,rlev);
[ap,ep,Tp,thp] = cartes2el(1,q,p,t);
fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
    t,E0,levc,ap,ep,Tp,thp,0,0);
t1 = cputime;
loop = 0;
loopt = 0;
stepn = 0;
levu = levc;
rd = 0;

while(t < tmaxp)
    stepn = stepn + 1;
%     t/tmaxp
%     [q,p,levc,casen,levu] = mapAnd(q,p,levc,rlev,hlev,steps);
%     [q,p,levu] = mapsymba(q,p,rlev,hlev,hsub);
%     [q,p,levc] = mapbasic(q,p,levc,rlev,hlev,steps);
levuo = levu;
%     [q,p,levv,levc,levu,rd] = mapadstep(q,p,levv,levc,rlev,hlev,hsub);

    [q,p,levc,loop,steps,levu] = globalstep(q,p,h,levc,rlev,hlev,t);
%     loopt = loopt + rd;

%     t = t + hlev(levu);
    t = t + h;
    if(abs(levuo - levu) > 0)
        levuo
        levu
    end
    
    if(t > tout(ind))
%         t/tmaxp
        ind = ind + 1;
        [E] = consq(q,p);
        dE = (E - E0)/E0;
        [ap,ep,Tp,thp] = cartes2el(1,q,p,t);
        fprintf(fcons,['%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n'], ...
            t,E,levc,ap,ep,Tp,thp,loop,stepn);
    end

end
fclose('all');
t2 = cputime;
dt = t2-t1;
[E] = consq(q,p);
dE = (E - E0)/E0;
dt
dE
loopt







end



function [q,p] = Tdrift(q,p,h)
q = q+ h*p;

end


function [q] = Tdriftr(q,p,h)
% global DRIFTC
q = q+ h*p;

% DRIFTC = DRIFTC + 1;

end


function [q,p] = Pkick(q,p,h)

r = sqrt(q(1)^2 + q(2)^2);
p = p - h*q/r^3;


end

function [p] = Pkickr(q,p,h,rlev,lev)
if(lev == 1)
    r1 = rlev(1);
    r2 = rlev(2);
    r = sqrt(q(1)^2 + q(2)^2);
    if(r > r1)
        force = -q/r^3;
    elseif(r2 <= r && r < r1)
        x = (r1-r)/(r1-r2);
        fct = 2*x^3 - 3*x^2 + 1;
        force = -q/r^3*fct;
    else
        return %first level makes no contribution
    end
else
    r1 = rlev(lev-1);
    r2 = rlev(lev);
    r3 = rlev(lev+1);
    r = sqrt(q(1)^2 + q(2)^2);
    if(r2 <= r && r < r1)
        x = (r1-r)/(r1-r2);
        fct = 2*x^3 - 3*x^2 + 1;
        force = -q/r^3*(1-fct);
    elseif(r3 <= r && r < r2)
        x = (r2-r)/(r2-r3);
        fct = 2*x^3 - 3*x^2 + 1;
        force = -q/r^3*fct;
    else
        return
    end
end

p = p + h*force;


end

function [q,p] = leap(q,p,stepsc,hlevc)
for i=1:stepsc
    [q,p] = Pkick(q,p,hlevc/2);
    [q,p] = Tdrift(q,p,hlevc);
    [q,p] = Pkick(q,p,hlevc/2);
end


end


function [q,p,levels,steps] = basic_global_step(q,p,h0,levc,collect_levels,rlev,hlev)
levels = 0;
steps = round(h0/hlev(levc));
% steps
hlevc = hlev(levc);

    [levc] = calclev(q,rlev);
    levels(1) = levc;

for i=1:steps


    [q,p] = Pkick(q,p,hlevc/2);
    [q,p] = Tdrift(q,p,hlevc);
    [q,p] = Pkick(q,p,hlevc/2);

    [levc] = calclev(q,rlev);
    levels(i+1) = levc;



end

%     [levc] = calclev(q,rlev);
%     levels(2) = levc;







end






function [q,p] = leaprecursive(q,p,hlev,levc,rlev,hsub)

for i=1:hsub
    [q,p] = Pkick(q,p,hlevc/2);
    [q,p] = Tdrift(q,p,hlevc);
    [q,p] = Pkick(q,p,hlevc/2);
    [levcp] = calclev(q1,rlev);
    if(levcp > levc) %going higher as much as needed; if level is smaller, do nothing, will need to repeat all step
        [q,p] = leaprecursive(q,p,hlev,levc,rlev,hsub);
    end
end

end


function [q,p,levmax] = driftop(q,p,lev,rlev,hlev,hsub,levmax)
% lev is the next one we're attempting
[flg] = rmincheck(q,p,hlev,rlev,lev-1); %if flg = 1, we're inside next radius shell during the next step, can't get back out...

if(flg == 1)
    for i=1:hsub
        if(lev > levmax)
            levmax = lev;
        end
        [p] = Pkickr(q,p,hlev(lev)/2,rlev,lev);
        [q,p,levmax] = driftop(q,p,lev+1,rlev,hlev,hsub,levmax);
        [p] = Pkickr(q,p,hlev(lev)/2,rlev,lev);
    end
elseif(flg == 0)
    [q] = Tdriftr(q,p,hlev(lev-1)); %decrease by one because added one before
end


end





function [levc] = calclev(q,rlev)

si = size(rlev,2);
r = sqrt(q(1)^2 + q(2)^2);
if(r > rlev(1)) %if greater than first curoff, then done
    levc = 1;
    return
end
for i=1:si
    if((r <= rlev(i)) && (r > rlev(i+1)) )
        levc = i+1;
        break
    end


end

end


function [flag] = levcheck(q,rlev,lev)

si = size(rlev,2);
r = sqrt(q(1)^2 + q(2)^2);
if(r > rlev(2)) %if greater than first curoff, then done
    levc = 1;
    return
end
for i=2:si
    if((r <= rlev(i-1)) && (r > rlev(i+1)) )
        levc = i;
        break
    end
end

if(levc == lev)
    flag = 1;
else
    flag = 0;
end

end



function [q,p,levc,casen,levu] = mapAnd(q,p,levc,rlev,hlev,steps)

% [q1,p1] = leapr(q,p,1,rlev,hlev,hsub,levc,levmax);
% levc
[q1,p1] = leap(q,p,steps(levc),hlev(levc));
[levcp] = calclev(q1,rlev);

if(levcp == levc) %step is the same after the fact
    levc = levcp;
    q = q1;
    p = p1;
    casen = 1;
    levu = levc;
    return
elseif(levcp > levc ) %step is now bigger; keep result
%     levcp - levc
    levc = levcp;
    q = q1;
    p = p1;
    casen=2;
    levu = levc;
    return
elseif(levcp < levc) %step is smaller; and means always smaller
%     [q1,p1] = leapr(q,p,1,rlev,hlev,hsub,levcp,levmax);
%     levcp - levc
% levcp
    [q1,p1] = leap(q,p,steps(levcp),hlev(levcp));
    levc = levcp;
    q = q1;
    p = p1;
    casen = 3;
    levu = levcp;
    return

end




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q,p,levv,levc,levu,rd] = mapadstep(q,p,levv,levc,rlev,hlev,hsub)

levu = levc;
[q1,p1] = leap(q,p,1,hlev(levc));
[levcp] = calclev(q1,rlev);
rd = 0;
if(levcp > levc)
    rd = 1;
    [q1,p1] = leap(q,p,1,hlev(levcp));
    levu = levcp;
    levv(levcp) = levv(levcp) + 1; %add step after attempt
    if(levv(levcp) == hsub)
        levv(levcp) = 0;
    end
    q = q1;
    p = p1;
    levc = levcp;
else %case when lower or equal level
    levv(levc) = levv(levc) + 1;
    if(levv(levc) == hsub)
        levv(levc) = 0;
    end
    if(levcp < levc) 
%         levcp
        levt = levc;
        while(levt > levcp)
            if(levv(levt) == 0)
                levt = levt - 1;
            else
                levt = levt;
                break;
            end
        end
%         levt
%         blergh
        levcp = levt; %go as low as possible in level
    end %endif
    q = q1;
    p = p1;
    levc = levcp; %may be levt or levcp from first lines (same as initial)
end





end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [q,p,levc1,loop,steps,levc] = globalstep(q,p,h0,levc,rlev,hlev,t)
loop = 0;
%sum of steps will still be h0
[q1,p1,levels,steps] = basic_global_step(q,p,h0,levc,1,rlev,hlev);
iN = max(levels); %if iN > levc, redo
check_again = iN > levc + 1; %or statement when two checks...
loop = 0;
if(iN > levc)
    while(iN > levc)
        loop = loop+1;
        levc = iN;
        [q1,p1,levels,steps] = basic_global_step(q,p,h0,levc,1,rlev,hlev);
        if(check_again)
            iN = max(levels);
        end
        check_again = iN > levc+1;
    end
    if(loop > 1)
        loop
    end
    levc1 = levels(end);
    q = q1;
    p = p1;
else
    levc1 = levels(end);
    q = q1;
    p = p1;
end




end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [q,p,levc,levu,flg] = mapAndsimple(q,p,rlev,levc,hlev)
% levc next step
% levu used this step
flg = 0;
[q1,p1] = leap(q,p,1,hlev(levc));
[levcp] = calclev(q1,rlev);
levu = levc;


if(levcp >= levc ) %step is now smaller; keep result
    levc = levcp;
    q = q1;
    p = p1;
    return
elseif(levcp < levc) %step is bigger;
    [q1,p1] = leap(q,p,1,hlev(levcp));
    [levcpp] = calclev(q1,rlev);
    if(levcpp ~= levcp)
        flg = 1;
    end
    levc = levcp;
    levu = levcp;
    q = q1;
    p = p1;
    return
end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [q,p,levmax] = mapAndsub(q,p,levc,rlev,hlev,hsub)


[q,p] = leaprecursive(q,p,hlev,levc,rlev,hsub);
%M^i steps



end


function [q,p,levmax] = mapsymba(q,p,rlev,hlev,hsub)

[p] = Pkickr(q,p,hlev(1)/2,rlev,1);
[q,p,levmax] = driftop(q,p,2,rlev,hlev,hsub,1);
% if(levmax > 5)
%     levmax
%     norm(q)
%     rlev
%     blergh
%     
% end
[p] = Pkickr(q,p,hlev(1)/2,rlev,1);

end







function [flg] = rmincheck(q,p,hlev,rlev,lev)

flg = 0;
qm = norm(q);

rcrit = rlev(lev); %check rcrit and dt during current, larger shell size
dt = hlev(lev); 

if(qm < rcrit) %already too small, this is also the case if they're moving away
    flg = 1;
%     display('case 1')
%     qm
%     rcrit
    return
end

a = dot(q,p);
if(a < 0) %particles moving away from each other, can keep
    vdotr = dot(q,p);
    v2 = sum(p.^2);
    q2 = sum(q.^2);
    tmin = -vdotr/v2;
    if(tmin < dt)
        r2min = q2 - vdotr^2/v2; %seems a very rare case, and doesn't even need to reduce level, when it does happen
%         display('case 2')
%         lev
%         sqrt(r2min)
%         rcrit
%         tmin
%         dt
%         blergh
    else
        r2min = q2 + 2*vdotr*dt + v2*dt*dt;
%         display('case 3')
%         sqrt(r2min)
    end
    rmin = sqrt(r2min);
    rmin = min(qm,rmin); %Symba check
    if(rmin < rcrit)
        flg = 1;
%         display('case 4')
%         rcrit
%         rmin
    else
        flg = 0;
%         display('case 5')
%         rcrit 
%         rmin
    end
    return
    
end

% display('case 6')
% qm
% rcrit
% a


end






function [q,p,levc] = mapbasic(q,p,levc,rlev,hlev,steps)

[q1,p1] = leap(q,p,steps(levc),hlev(levc));
[levcp] = calclev(q1,rlev);
levc = levcp;
q = q1;
p = p1;
return



end



function E = consq(q,p)

E = 1/2*(p(1)^2 + p(2)^2) - 1/sqrt(q(1)^2 + q(2)^2);

end



function [a,e,T,th] = cartes2el(gm,x,v,t)
%x is 2D, v is 2D

del = 1e-10;
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
if(T < 0)
    T = T + P;
end
% blergh


end



