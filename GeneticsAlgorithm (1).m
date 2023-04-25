%data input
clear
clc
A = [1907290000 1622435000 1758670000 1832980000 1849080000 19856350000 1744080000 1895867000];
B = [235315000 200637000 193206000 220453000 221615000 186426000 180453000 210453000];
C1 = [117900 30300 40200 38100 66300 85800 66000 61500; 91200 45900 111900 63000 39600 75600 55800 49800; 62700 85500 117300 68400 47100 84000 92400 91500];
L1 = [39.3 10.1 13.4 12.7 22.1 28.6 22 20.5; 30.4 15.3 37.3 21 13.2 25.2 18.6 16.6; 20.9 28.5 39.1 22.8 15.7 28 30.8 30.5];
C2 = [68700 64500 73500 67800 72300 73200 66300 7890000000 88200 61800 66600 69600; 57900 59400 50700 55800 53700 51300 57300 68700 67800 36600 56400 59400; 26400 26100 25500 27000 30000 31500 24000 37500 93600 64800 38700 46800;62700 65100 83100 61200 59400 52800 66300 77100 75000 44100 77400 101100;67500 65700 62700 60900 58500 53400 65700 75900 47400 10500 85200 100500; 21600 26100 24600 27000 29400 33000 21900 10800 106800 75900 30300 23400; 30600 35100 33600 36000 38100 42000 30900 19500 115800 84900 44100 5100; 36000 24600 36600 27300 31200 35100 28200 39300 108300 95400 16800 19800];
L2 = [22.9 21.5 24.5 22.6 24.1 24.4 22.1 26.3 29.4 20.6 22.2 23.2; 19.3 19.8 16.9 18.6 17.9 17.1 19.1 22.9 22.6 12.2 18.8 19.8; 8.8 8.7 8.5 9 10 10.5 8 12.5 31.2 21.6 12.9 15.6; 20.9 21.7 27.7 20.4 19.8 17.6 22.1 25.7 25 15.7 25.8 33.7; 22.5 21.9 20.9 20.3 19.5 17.8 21.9 25.3 15.8 3.5 28.4 33.5; 7.2 8.7 8.2 9 9.8 11 7.3 3.6 35.6 25.3 10.1 7.8; 10.2 11.7 11.2 12 12.7 14 10.3 6.5 38.6 28.3 14.7 1.7; 12 8.2 12.2 9.1 10.4 11.7 9.4 13.1 36.1 31.8 5.6 6.6];
demand = [1358 1267 1035 392 273 5100 5050 1590 1570 4400 2300 1350];
f = [173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000;173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000 173000];
s = [1 2 3]; %supplier
v = 45;
t = L2/v;
T = [0.42 0.33 0.33 0.42 0.38 0.33 0.25 0.38 0.58 0.42 0.42 0.25];
Tmax = [0.58 0.5 0.5 0.58 0.55 0.5 0.42 0.55 0.75 0.58 0.58 0.42];
Tmin = [0.25 0.17 0.17 0.25 0.22 0.17 0.08 0.22 0.42 0.25 0.25 0.08];
fi = [0.05 0.04 0.08]; %contaminate rate of freshfood
alpha1 = 0.7;
alpha2 = 0.3;
r = [126000 195000 215000]; %price of freshfood 

%% initialization
pmax = 14;
population = initialization(pmax);
for a = 1:pmax
    ar1 = population(a,(1:23));
    [d4,consumer] = output3(ar1);
    [k1,k2,k3] = splitout3(d4);
    [dbn01,dbn02,dbn03] = intonumarr1(k1,k2,k3);
    dcp = [dbn01,dbn02,dbn03];
    % fitnessfunction
    op1 = myf1(A,B,dcp);
    op2 = myf2(C1,C2,s,L1,L2,dcp,demand,consumer);
    op3 = myf3(f,dcp,t,T,consumer);
    [op4,corruptp] = myf4(t,dcp,alpha1,alpha2,L1,L2,v,fi,consumer,s,r);

    totalPop = op1 + op2 + op3 + op4;
    freshnessPop = sum(corruptp);
    DCservicetimePop = myf5(consumer,t,Tmin,Tmax,dcp);
    DCtimePop = sum(DCservicetimePop);

    finalobjPop = 0.6*totalPop + 0.2*freshnessPop + 0.2*DCtimePop;
    cstPop(a,1) = totalPop; %,freshnessPop,DCtimePop,finalobjPop];    
end
%% sorting for population 
[x5,I1] = sort(cstPop,"descend");
sortpop = population(I1,:);
%% solution for each array
iter = 0;
while iter <= 5000
if(iter >= 0 && iter <= 5000)       
%% Main loop
for c = 1:10
    out1 = 0;
    check = 0;
    feasible = 0;
    %crossover DC
        while check == 0
            [r1,r2] = randm();
            check = choosing(r1,r2);
        end
        [parent1,parent2] = choosing2(population,r1,r2); 
        [x,y] = splitintodcarray(parent1,parent2);
        [childrendc1,childrendc2] = crs(x,y);
    %muitation DC
        while out1 == 0 || feasible == 0 
            ran = randomize(10);
        while ran > 0
            m1 = round(random('Uniform',1,10)-0.5);
            m2 = round(random('Uniform',1,3)-0.5);
            if(m2 == 1)
                if (childrendc1(m1) == 0)
                    childrendc1(m1) = 1;
                else 
                    childrendc1(m1) = 0;
                end
            end
            if(m2 == 2)
                if (childrendc2(m1) == 0)
                    childrendc2(m1) = 1;
                else 
                    childrendc2(m1) = 0;
                end
            end
            ran = ran - 1;
        end
        [dc1,dc2,dc3] = splitdc(childrendc1,childrendc2);
        [dcn1, dcn2, dcn3] = dctonum(dc1,dc2,dc3);
        feasible = myfunc4(dcn1,dcn2,dcn3);
        out1 = myfunc3(dcn1,dcn2,dcn3);
        dcnchil = [dc1,dc2,dc3];
        dcn = [dcn1,dcn2,dcn3];
        
        end  

%crossover customer 
[p1,p2] = splitcus(parent1,parent2);
ran2 = randomized(15);
v1 = cross(p1,p2,ran2);
%muitation customer
r3 = randperm(14);
v5 = transfering(v1,r3);
%fitnessfunction
o1 = cost1(A,B,dcn);
o2 = cost2(C1,C2,s,L1,L2,dcn,demand,v5);
o3 = cost3(f,dcn,t,T,v5);
[o4,corrupt] = cost4(t,dcn,alpha1,alpha2,L1,L2,v,fi,v5,s,r);

totalcost = o1 + o2 + o3 + o4;

freshness = sum(corrupt);
DCservicetime = sv(v5,t,Tmin,Tmax,dcn);
DCtime = sum(DCservicetime);

seq = array(dcnchil,v5);

finalobj = 0.6*totalcost + 0.2*freshness + 0.2*DCtime;
%building population
chil(c,(1:23)) = seq(1,(1:23));
cst(c,1) = totalcost; %,freshness,DCtime,finalobj];
end
%% sorting children
[y5,I] = sort(cst,"ascend");
sortchil = chil(I,:);
%% replace and perform new population
j = 1;
while j <= 7
      if(x5(1) > y5(j))
         x5(1) = y5(j);
         sortpop(1,(1:23)) = sortchil(j,(1:23));
      end
      [temp,I2] = sort(x5,"descend");
      sortpop1 = sortpop(I2,:);
      x5 = temp;
      sortpop = sortpop1;
      j = j + 1;
end
population = sortpop;
bestcost = x5(pmax);
bestsolution = population(pmax,(1:23));
iter = iter + 1;
disp(x5);
end
end
%% final array of DC - customer
function seq = array(dcn,v5)
    seq = [dcn,v5];
end
%% forming right form sequence from population
function [d4,consumer] = output3(ar1) 
    d4 = ar1(1:9);
    consumer = ar1(10:23);
end
function [k1,k2,k3] = splitout3(d4)
        k1((1:3)) = d4((1:3));
        k2((1:3)) = d4((4:6));
        k3((1:3)) = d4((7:9));
end
function [dbn01,dbn02,dbn03] = intonumarr1(k1,k2,k3)
dbn01 = 0;
dbn02 = 0;
dbn03 = 0;
    for i = 0:2
            dbn01 = dbn01 + (2^i)*k1(i+1);
            dbn02 = dbn02 + (2^i)*k2(i+1);
            dbn03 = dbn03 + (2^i)*k3(i+1);
    end
end

%% population
function population = initialization(pmax)
sas = 0;
feasiblepop = 0;
p = 1;
    while p <= pmax || sas == 0 || feasiblepop == 0
            population(p,(1:23)) = [randi([0 1],1,9),randperm(14,14)];
            for j = 10:23
                if(population(p,j) == 13 || population(p,j) == 14)
                    population(p,j) = 0;
                end
            end
            [d,cus] = splitarr(population,p);
            [db1,db2,db3] = splitout2(d,p);
            [dbn1,dbn2,dbn3] = intonumarr(db1,db2,db3,p);
            feasiblepop = checkingout2(dbn1,dbn2,dbn3);
            sas = checkingout(dbn1,dbn2,dbn3);
            if(feasiblepop == 1 && sas == 1)
                p = p + 1;
            else
                p = p;
            end
    end
end
%% split into cus and dc array
function [d,cus] = splitarr(population,p)
    d(p,(1:9)) = population(p,(1:9));
    cus(p,(1:14)) = population(p,(10:23));
end
%% split DC for population for population
function [db1,db2,db3] = splitout2(d,p)
        db1(p,(1:3)) = d(p,(1:3));
        db2(p,(1:3)) = d(p,(4:6));
        db3(p,(1:3)) = d(p,(7:9));
end
%% changing into dc-number for population
function [dbn1,dbn2,dbn3] = intonumarr(db1,db2,db3,p)
dbn1 = 0;
dbn2 = 0;
dbn3 = 0;
    for i = 0:2
            dbn1 = dbn1 + (2^i)*db1(p,i+1);
            dbn2 = dbn2 + (2^i)*db2(p,i+1);
            dbn3 = dbn3 + (2^i)*db3(p,i+1);
    end
end
%% checking for population
function sas = checkingout(dbn1,dbn2,dbn3)
            if(dbn1 == dbn2)
                sas = 0;
            elseif(dbn2 == dbn3)
                sas = 0;
            elseif(dbn1 == dbn3)
                sas = 0;
            else
                sas = 1;
            end
end
function feasiblepop = checkingout2(dbn1,dbn2,dbn3) 
        if(dbn1 == 0)   
            feasiblepop = 0;
        elseif(dbn2 == 0)
            feasiblepop = 0;
        elseif(dbn3 == 0)
            feasiblepop = 0;
        else
            feasiblepop = 1;
        end
       
end
%% fitness function for population
%% operation cost
function op1 = myf1(A,B,dcp) 
k = 0;
k1 = 0;
     for i = 1:3 
         k = k + A(1,dcp(i));
         k1 = k1 + B(1,dcp(i));
         op1 = k + k1; %total fixed and operation cost
     end
end
%% transportation cost
function op2 = myf2(C1,C2,s,L1,L2,dcp,demand,consumer)
trans1 = 0;
for h = 1:3
    for i = 1:3
       trans1 = trans1 + C1(s(h),dcp(i))*L1(s(h),dcp(i));
    end
end
trans2 = zeros(3,length(consumer));
pos = 1;
i = 1;
while pos <= length(consumer) 
    while consumer(pos) == 0
        if(pos < length(consumer))
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(consumer) && consumer(pos) == 0)
        i = i + 1;
        pos = pos;
        trans2(i,pos) = 0;
        break
    end
    trans2(i,pos) = C2(dcp(i),consumer(pos))*L2(dcp(i),consumer(pos))*demand(consumer(pos));
    pos = pos + 1;
end
trans2 = sum(sum(trans2));
op2 = trans1 + trans2; %total of transportation cost
end
%% Delayed penalty cost
function op3 = myf3(f,dcp,t,T,consumer)
i = 1;
pos = 1;
p1 = zeros(3, length(consumer));
while pos <= length(consumer)
    while consumer(pos) == 0
        if(pos < length(consumer))
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(consumer) && consumer(pos) == 0)
        i = i + 1;
        pos = pos;
        p1(i,pos) = 0;
        break
    end  
    if t(dcp(i),consumer(pos)) > T(consumer(pos))
        p1(i,pos) = f(dcp(i),consumer(pos)) * (t(dcp(i),consumer(pos)) - T(consumer(pos)));
    else
        p1(i,pos) = 0;
    end
     pos = pos + 1;
end
op3 = sum(sum(p1));
end
%% freshness cost
function [op4,corruptp] = myf4(t,dcp,alpha1,alpha2,L1,L2,v,fi,consumer,s,r)
lnew = zeros(3,length(consumer));
pos = 1;
i = 1;
h = 1;
while pos <= length(consumer)
    while consumer(pos) == 0
        if(pos < length(consumer))
            h = h + 1;
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(consumer) && consumer(pos) == 0)
        i = i + 1;
        h = h + 1;
        pos = pos;
        lnew(i,pos) = 0;
        break
    end 
        lnew(i,pos) = (alpha1*((L1(s(h),dcp(i))/v) + t(dcp(i),consumer(pos))) + alpha2*(L1(s(h),dcp(i)) + L2(dcp(i),consumer(pos))));
        pos = pos + 1; 
end
lnew = sum(sum(lnew));
op4 = 0;
for k = 1:3
     corruptp(k) = (1 - fi(k)).^lnew; % total corruption cost
     op4 = op4 + r(k)*(1 - corruptp(k));
end
end
%% Distribution center service time
function DCservicetimePop = myf5(consumer,t,Tmin,Tmax,dcp) 
i = 1;
pos = 1;
omega = zeros(3,length(consumer));
    while pos <= length(consumer)
        while consumer(pos) == 0
            if(pos < length(consumer))
                i = i + 1;
                pos = pos + 1;
            else
                break
            end        
        end
        if(pos == length(consumer) && consumer(pos) == 0)
            i = i + 1;
            pos = pos;
            omega(i,pos) = 0;
            break
        end 
        if(t(dcp(i),consumer(pos)) <= Tmin(consumer(pos)))
            omega(i,pos) = 1;
        elseif (t(dcp(i),consumer(pos)) > Tmin(consumer(pos)) && t(dcp(i),consumer(pos)) < Tmax(consumer(pos)))
            omega(i,pos) = ((Tmax(consumer(pos)) - t(dcp(i),consumer(pos)))/(Tmax(consumer(pos)) - Tmin(consumer(pos))));
        else 
            omega(i,pos) = 0;
        end
        pos = pos + 1;       
    end
    DCservicetimePop = sum(omega);
end




%% random row
function [r1,r2] = randm()
    r1 = round(random("Uniform",1,10)-0.5);
    r2 = round(random("Uniform",1,10)-0.5);
end

%% choosing
function check = choosing(r1,r2)
    if(r1 == r2)
        check = 0;
    else 
        check = 1;
    end   
end
function [parent1,parent2] = choosing2(population,r1,r2)
    parent1 = population(r1,(1:23));
    parent2 = population(r2,(1:23));
end
%% split Dc from array parrent
function [x,y] = splitintodcarray(parent1,parent2)
    for i = 1:9
        x(i) = parent1(i);
        y(i) = parent2(i);
    end
end
%% crossover for DC
function [childrendc1,childrendc2] = crs(x,y)
    for i = 1:5
        childrendc1(i) = x(i); 
        childrendc1(9 + 1-i) = y(9 + 1-i);
        childrendc2(i) = y(i);
        childrendc2(9 + 1-i) = x(9 + 1-i);
    end
end
%Muitation for DC
%randomize
function ran = randomize(ub)
    ran = round(random('Uniform',0,ub)-0.5);
end
%% split DC
function [dc1,dc2,dc3] = splitdc(childrendc1,childrendc2)
ran3 = round(random('Uniform',1,3)-0.5);
if (ran3 == 1)
    for i = 1:3
      dc1(i) = childrendc1(i);
      dc2(i) = childrendc1(i + 3);
      dc3(i) = childrendc1(i + 6);
    end
end
if(ran3 == 2)
    for i = 1:3
      dc1(i) = childrendc2(i);
      dc2(i) = childrendc2(i + 3);
      dc3(i) = childrendc2(i + 6);
    end
end
end
%% changing dc to num dc
function [dcn1, dcn2, dcn3] = dctonum(dc1,dc2,dc3)
dcn1 = 0;
dcn2 = 0;
dcn3 = 0;
    for i = 0:2
     dcn1 = dcn1 + (2^i)*dc1(i+1);
     dcn2 = dcn2 + (2^i)*dc2(i+1);
     dcn3 = dcn3 + (2^i)*dc3(i+1);
    end
end

function out1 = myfunc3(dcn1,dcn2,dcn3)
             if(dcn1 == dcn2)
                out1 =0;
            elseif (dcn2 == dcn3)
                out1 = 0;
            elseif (dcn1 == dcn3)
                out1 = 0;
            else 
                out1 = 1;
             end 
end

function feasible = myfunc4(dcn1,dcn2,dcn3)
        if(dcn1 == 0)
            feasible = 0;
        elseif (dcn2 == 0)
            feasible = 0;
        elseif (dcn3 == 0)
            feasible = 0;
        else 
            feasible = 1;
        end
end

%% split customer array
function [p1,p2] = splitcus(parent1,parent2)
    for i = 1:14
      p1(i) = parent1(i + 9);
      p2(i) = parent2(i + 9);
    end
end
%% randomize position
function ran2 = randomized(ub)
         ran2 = round(random('Uniform',1,ub)-0.5);
end
%% crossover
function v1 = cross(p1,p2,ran2)
         r4 = round(random('Uniform',1,3)-0.5);
            if (r4 == 1)
                v1 = cat(2,p1(ran2 :end),p1(1:ran2 -1)); 
            elseif (r4 == 2)
                v1 = cat(2,p2(ran2 :end),p2(1:ran2 -1));
            end
end
%% muitation
function v5 = transfering(v1,r3) 
    for i = 1:14
        v5(i) = v1(r3(i));
    end       
end
%% fixed cost and operation cost
function o1 = cost1(A,B,dcn)
k = 0 ;
k1 = 0;
     for i = 1:3 
         k = k + A(1,dcn(i));
         k1 = k1 + B(1,dcn(i));
         o1 = k + k1; %total fixed and operation cost
     end
end

%% transportation cost
function o2 = cost2(C1,C2,s,L1,L2,dcn,demand,v5)
trans1 = 0;
for h = 1:3
    for i = 1:3
       trans1 = trans1 + C1(s(h),dcn(i))*L1(s(h),dcn(i));
    end
end
trans2 = zeros(3,length(v5));
pos = 1;
i = 1;
while pos <= length(v5) 
    while v5(pos) == 0
        if(pos < length(v5))
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(v5) && v5(pos) == 0)
        i = i + 1;
        pos = pos;
        trans2(i,pos) = 0;
        break
    end
    trans2(i,pos) = C2(dcn(i),v5(pos))*L2(dcn(i),v5(pos))*demand(v5(pos));
    pos = pos + 1;
end
trans2 = sum(sum(trans2));
o2 = trans1 + trans2; %total of transportation cost
end
%% Delayed penalty cost
function o3 = cost3(f,dcn,t,T,v5)
i = 1;
pos = 1;
p1 = zeros(3, length(v5));
while pos <= length(v5)
    while v5(pos) == 0
        if(pos < length(v5))
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(v5) && v5(pos) == 0)
        i = i + 1;
        pos = pos;
        p1(i,pos) = 0;
        break
    end  
    if t(dcn(i),v5(pos)) > T(v5(pos))
        p1(i,pos) = f(dcn(i),v5(pos)) * (t(dcn(i),v5(pos)) - T(v5(pos)));
    else
        p1(i,pos) = 0;
    end
     pos = pos + 1;
end
o3 = sum(sum(p1));
end

%% Freshgood corruption cost
function [o4,corrupt] = cost4(t,dcn,alpha1,alpha2,L1,L2,v,fi,v5,s,r)
lnew = zeros(3,length(v5));
pos = 1;
i = 1;
h = 1;
while pos <= length(v5)
    while v5(pos) == 0
        if(pos < length(v5))
            h = h + 1;
            i = i + 1;
            pos = pos + 1;
        else
            break
        end        
    end
    if(pos == length(v5) && v5(pos) == 0)
        i = i + 1;
        h = h + 1;
        pos = pos;
        lnew(i,pos) = 0;
        break
    end 
        lnew(i,pos) = (alpha1*((L1(s(h),dcn(i))/v) + t(dcn(i),v5(pos))) + alpha2*(L1(s(h),dcn(i)) + L2(dcn(i),v5(pos))));
        pos = pos + 1; 
end
lnew = sum(sum(lnew));
o4 = 0;
for k = 1:3
     corrupt(k) = (1 - fi(k)).^lnew; % total corruption cost
     o4 = o4 + r(k)*(1 - corrupt(k));
end
end
%% distribution center service time
function DCservicetime = sv(v5,t,Tmin,Tmax,dcn) 
i = 1;
pos = 1;
omega = zeros(3,length(v5));
    while pos <= length(v5)
        while v5(pos) == 0
            if(pos < length(v5))
                i = i + 1;
                pos = pos + 1;
            else
                break
            end        
        end
        if(pos == length(v5) && v5(pos) == 0)
            i = i + 1;
            pos = pos;
            omega(i,pos) = 0;
            break
        end 
        if(t(dcn(i),v5(pos)) <= Tmin(v5(pos)))
            omega(i,pos) = 1;
        elseif (t(dcn(i),v5(pos)) > Tmin(v5(pos)) && t(dcn(i),v5(pos)) < Tmax(v5(pos)))
            omega(i,pos) = ((Tmax(v5(pos)) - t(dcn(i),v5(pos)))/(Tmax(v5(pos)) - Tmin(v5(pos))));
        else 
            omega(i,pos) = 0;
        end
        pos = pos + 1;       
    end
    DCservicetime = sum(omega);
end






