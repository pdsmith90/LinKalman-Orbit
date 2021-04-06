function [ routclean,voutclean, routmeas,voutmeas ] = OrbitPropagate( r0,v0,dt )
global mu J2 Re sig

% a2b=-mu/norm(ro).*r
% aJ2=-(3.*J2.*mu.*Re.^2)/(2*norm(r)^5).*((1-(5.*r(3)^2)/norm(r)^2)+[0,0,2*r(3)])
% ORBIT = @(a2b,aJ2) a2b+aJ2;

% Function Handles

%the ORBIT function handle in terms of d(r,v)/dt = [dr/dt, dv/dt] = [v,a]
%with inputs of [ro,vo]

% acceleration due to J2
% 11/01/11 line 14/15 were corrected. mu term was added.
% 11/09/11 Runge Kutta code specific to Orbit Propagation was added.
%           Equations were cleaned up for simplicity
rvinput = [r0';v0'];
orbitFunc = @(rvinput,J2,mu,Re) [
    rvinput(4)
    rvinput(5)
    rvinput(6)
    (-rvinput(1)*mu/norm(rvinput(1:3))^3+(-.5*J2*mu*Re^2)*(-15*rvinput(1)*rvinput(3)^2/norm(rvinput(1:3))^7+3*rvinput(1)/norm(rvinput(1:3))^5))
    (-rvinput(2)*mu/norm(rvinput(1:3))^3+(-.5*J2*mu*Re^2)*(-15*rvinput(2)*rvinput(3)^2/norm(rvinput(1:3))^7+3*rvinput(2)/norm(rvinput(1:3))^5))
    (-rvinput(3)*mu/norm(rvinput(1:3))^3+(J2*mu*Re^2/2)*(15*rvinput(3)^3/norm(rvinput(1:3))^7 - 9*rvinput(3)/norm(rvinput(1:3))^5))
    ];

[rvout] = orbitRK4(orbitFunc,rvinput,dt);
routclean = rvout(1:3);
voutclean = rvout(4:6);

% add in gps measurement noise
routmeas=routclean.*(1+sig.r.*(rand(1)-0.5));
voutmeas=voutclean.*(1+sig.v.*(rand(1)-0.5));

end

function [rvout] = orbitRK4(orbitFunc,rvinput,dt)
global mu J2 Re

k1 = dt.*orbitFunc(rvinput,J2,mu,Re);
k2 = dt.*orbitFunc(rvinput+k1/2,J2,mu,Re);
k3 = dt.*orbitFunc(rvinput+k2/2,J2,mu,Re);
k4 = dt.*orbitFunc(rvinput+k3,J2,mu,Re);

rvout = rvinput + (1/6).*(k1+2.*k2+2.*k3+k4);
end
