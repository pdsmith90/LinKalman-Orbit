function [ out ] = CoordTransform( in,intype,outtype )
% intype supported:  string: 'ECI','OE'
% outtype supported: string: 'ECI','OE','LVLH'
% ECI:  [r,v]==[rx,ry,rz,vx,vy,vz]
% OE:   [e,a,i,om,w,nu], with angles given in degrees
% LVLH: [i,j,k] in ECI

global mu

if strcmp(intype,'OE') && strcmp(outtype,'OE')
    out=in;
    return
elseif strcmp(intype,'OE') % OE to ECI
    e=in(1);
    a=in(2);
    i=in(3).*pi./180;
    om=in(4).*pi./180;
    w=in(5).*pi./180;
    nu=in(6).*pi./180;    
    p=a*(1-e^2);
    rpqw=[(p*cos(nu)/(1+e*cos(nu))), (p*sin(nu)/(1+e*cos(nu))), 0];
    vpqw=[(-sqrt(mu/p)*sin(nu)),(sqrt(mu/p)*(e+cos(nu))),0];
    rot=[[cos(om)*cos(w)-sin(om)*sin(w)*cos(i), -cos(om)*sin(w)-sin(om)*cos(w)*cos(i), sin(om)*sin(i)];
         [sin(om)*cos(w)+cos(om)*sin(w)*cos(i), -sin(om)*sin(w)+cos(om)*cos(w)*cos(i), -cos(om)*sin(i)];
         [sin(w)*sin(i), cos(w)*sin(i), cos(i)]];
    for i=1:size(rpqw,1)
    r=rpqw*rot;
    v=vpqw*rot;
    in=[r;v];
    end
end
% At this point, in variable is in the ECI frame
r=in(1:3);
v=in(4:6);

if strcmp(outtype,'ECI')
    out=in;
    return
elseif strcmp(outtype,'OE') % ECI to OE
    h=cross(r,v);
    n=cross([0,0,1],h);
    e=cross(in(4:6),h)/mu-r/norm(r);
    xi=.5*norm(v)^2-mu/norm(r);
    a=-mu/(2*xi);
    i=acos(h(3)/norm(h));
    om=acos(n(1)/norm(n));
    if n(2)<0
        om=2*pi-om;
    end
    w=acos(dot(n,e)/(norm(n)*norm(e)));
    if e(3)<0
        w=2*pi-w;
    end
    nu=acos(dot(e,r)/(norm(e)*norm(r)));
    if dot(r,v)<0
        nu=2*pi-nu;
    end
    angles=[i,om,w,nu];
    angles=angles.*180./pi;
    out=[norm(e),a,angles];
    return
elseif strcmp(outtype,'LVLH') % ECI to LVLH
    i=r/norm(r);
    h=cross(r,v);
    k=h/norm(h);
    j=cross(k,i);
    out=[i,j,k];
    return
end
end