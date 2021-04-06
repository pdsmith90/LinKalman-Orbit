function [] = OPKshell()
%%
intxt_physprop='OPK PhysProps.txt';
intxt_scenario='OPK Inputs.txt';
% Open phys props .txt and define global variables
inid=fopen(intxt_physprop);
in_physprop=textscan(inid,'%f%f%f%f%f%f%f','Commentstyle','%','Delimiter','\n','MultipleDelimsAsOne',1,'CollectOutput',0);
fclose(inid);
global mu J2 Re sig
mu=in_physprop{1};
J2=in_physprop{2};
Re=in_physprop{3};

sig.r=in_physprop{4}.*1e-3;% GPS position accuracy error, m to km
sig.v=in_physprop{5}.*1e-3;% GPS velocity accuracy error, m/s to km/s
% %unused gps data
% gerr.pp=in_physprop{6};% GPS measurement precision error
% gerr.hz=in_physprop{7};% GPS acquision rate

% Open input .txt and get cell array of input conditions
inid=fopen(intxt_scenario);
in_scenario=textscan(inid,'%f%s%s%s%s%s%f%f','Commentstyle','%','Delimiter','\t','MultipleDelimsAsOne',1,'CollectOutput',0);
fclose(inid);
% Preallocate
out=cell(1,length(in_scenario{1}));
timesv=cell(1,length(in_scenario{1}));
% Loop for each input condition
if length(in_scenario{1})==1 
    disp(strcat(['No Inputs Specified in ',intxt,', stopping']));return;
else
    for i=1:length(in_scenario{1})
    run=in_scenario{1}(i);
    % convert strings to arrays
    for j=[3,5,6]
    in_scenario{j}{i}=eval(in_scenario{j}{i});
    end
    % convert start time, dt, run time to a linspace of serial dates
    timesv{i}=datenum(in_scenario{6}{i})+[0,linspace(in_scenario{7}(i),in_scenario{8}(i),(in_scenario{8}(i)./in_scenario{7}(i)))./86400]';
       
    % Run propagation, save results to file
    out{i} = ADCSrun(in_scenario{2}{i},in_scenario{3}{i},in_scenario{4}{i},in_scenario{5}{i},timesv{i});
    mkdir('Results\');
    dlmwrite(strcat('Results\',intxt_scenario(1:(end-4)),'_',num2str(run),'.csv'),...
        'serial date,ri,rj,rk,rmi,rmj,rmk,rgi,rgj,rgk,vi,vj,vk,vmi,vmj,vmk,vgi,vgj,vgk,e,a,i,om,w,nu','delimiter','')
    dlmwrite(strcat('Results\',intxt_scenario(1:(end-4)),'_',num2str(run),'.csv'), [timesv{i},out{i}],'-append','delimiter',',');
    display(strcat('Data saved to ','Results\',intxt_scenario(1:(end-4)),'_',num2str(run),'.csv'))
    % Plot
    ADCSplot(strcat('Results\',intxt_scenario(1:(end-4)),'_',num2str(run),'.csv'),i);
    end
end
end


function out = ADCSrun(oitype,oi,kitype,ki,timesv) % include input for initial time
global mu J2 Re sig
%preallocate
emp=zeros(length(timesv),3);
r=emp;v=emp;
rm=emp;vm=emp;
rk=emp;vk=emp;
oe=[emp,emp];

% make sure input orbital, kalman is in ECI
if not(strcmp(oitype,'ECI'))
oi = CoordTransform( oi,oitype,'ECI' );
end
if not(strcmp(kitype,'ECI'))
ki = CoordTransform( ki,kitype,'ECI' );
end

% set first set of data points as inputs
r(1,:)=oi(1,:);v(1,:)=oi(2,:);
rm(1,:)=oi(1,:);vm(1,:)=oi(2,:);
rk(1,:)=ki(1,:);vk(1,:)=ki(2,:);

% loop orbit, kalman
dt=(timesv(2)-timesv(1))*86400;
for s=1:(length(timesv)-1)
%%    [ r((s+1),:),v((s+1),:),rm((s+1),:),vm((s+1),:) ] = OrbitPropagate( r(s,:),v(s,:),dt );
    rvinput = [r(s,:)';v(s,:)'];
    orbitFunc = @(rvinput,J2,mu,Re) [
    rvinput(4)
    rvinput(5)
    rvinput(6)
    (-rvinput(1)*mu/norm(rvinput(1:3))^3+(-.5*J2*mu*Re^2)*(-15*rvinput(1)*rvinput(3)^2/norm(rvinput(1:3))^7+3*rvinput(1)/norm(rvinput(1:3))^5))
    (-rvinput(2)*mu/norm(rvinput(1:3))^3+(-.5*J2*mu*Re^2)*(-15*rvinput(2)*rvinput(3)^2/norm(rvinput(1:3))^7+3*rvinput(2)/norm(rvinput(1:3))^5))
    (-rvinput(3)*mu/norm(rvinput(1:3))^3+(J2*mu*Re^2/2)*(15*rvinput(3)^3/norm(rvinput(1:3))^7 - 9*rvinput(3)/norm(rvinput(1:3))^5))
    ];

    %% [rvout] = orbitRK4(orbitFunc,rvinput,dt);
    k1 = dt.*orbitFunc(rvinput,J2,mu,Re);
    k2 = dt.*orbitFunc(rvinput+k1/2,J2,mu,Re);
    k3 = dt.*orbitFunc(rvinput+k2/2,J2,mu,Re);
    k4 = dt.*orbitFunc(rvinput+k3,J2,mu,Re);

    rvout = rvinput + (1/6).*(k1+2.*k2+2.*k3+k4);
    %%
    
    r((s+1),:) = rvout(1:3);
    v((s+1),:) = rvout(4:6);
    
    % add in gps measurement noise
    rm((s+1),:)=r((s+1),:).*(1+sig.r.*(rand(1)-0.5));
    vm((s+1),:)=v((s+1),:).*(1+sig.v.*(rand(1)-0.5));
%%
    
    [ rk((s+1),:),vk((s+1),:) ] = LinKalman( rm(s,:),vm(s,:),dt );
end

% get set of orbital elements
for s=1:length(timesv)
    oe(s,:)=CoordTransform([r(s,:),v(s,:)],'ECI','OE');
end
out=[r,rm,rk,v,vm,vk,oe];
end
%%
function [timesv,r,v,oe] = ADCSplot(incsv,run)
global Re
%read the results .csv
out=csvread(incsv,1,0);
timesv=out(:,1);
r=out(:,2:4);v=out(:,11:13);
rm=out(:,5:7);v=out(:,14:16);
rg=out(:,8:10);v=out(:,17:19);
oe=out(:,20:25);

%plot 1, orbit r
figure(run*5-4)
hold on
plot3(r(:,1),r(:,2),r(:,3))
[x,y,z]=sphere;
xe=Re.*x;
ye=Re.*y;
ze=Re.*z;
surf(xe,ye,ze)
alpha(0.4)
xlabel('i_E_C_I')
ylabel('j_E_C_I')
zlabel('k_E_C_I')
title(['Orbit, Run # ',num2str(run)]);
axis equal
view(0,25)
hold off
saveas((run*5-4),[incsv(1:(end-4)),'_3d orbit'],'png')

%plot 2, orbital elements
figure(run*5-3)
subplot(6,1,1)
plot(timesv,oe(:,1))
title(['Orbital Elements, Run # ',num2str(run)]);
ylabel('e')
datetick('x',0)
subplot(6,1,2)
plot(timesv,oe(:,2))
ylabel('a')
datetick('x',0)
subplot(6,1,3)
plot(timesv,oe(:,3))
ylabel('i')
datetick('x',0)
subplot(6,1,4)
plot(timesv,oe(:,4))
ylabel('OMEGA')
datetick('x',0)
subplot(6,1,5)
plot(timesv,oe(:,5))
ylabel('w')
datetick('x',0)
subplot(6,1,6)
plot(timesv,oe(:,6))
ylabel('nu')
datetick('x',0)
xlabel('Time')
saveas((run*5-3),[incsv(1:(end-4)),'_oe'],'png')

end