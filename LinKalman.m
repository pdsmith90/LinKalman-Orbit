function [ r1,v1 ] = LinKalman( r0,v0,dt )
global mu J2 Re sig
persistent Q R F P Xp
%% Initialization
if isempty(Q)==1 % Calculate process noise matrix.
    Q=[sig.r.^2 0 0 0 0 0;
         0 sig.r.^2 0 0 0 0;
         0 0 sig.r.^2 0 0 0;
         0 0 0 sig.v.^2 0 0;
         0 0 0 0 sig.v.^2 0;
         0 0 0 0 0 sig.v.^2];
end
if isempty(R)==1 % Calculate measurement noise matrix.
    R=[(sig.r.*(rand(1)-0.5)).^2 0 0 0 0 0;
        0 (sig.r.*(rand(1)-0.5)).^2 0 0 0 0;
        0 0 (sig.r.*(rand(1)-0.5)).^2 0 0 0;
        0 0 0 (sig.v.*(rand(1)-0.5)).^2 0 0;
        0 0 0 0 (sig.v.*(rand(1)-0.5)).^2 0;
        0 0 0 0 0 (sig.v.*(rand(1)-0.5)).^2];
end
if isempty(F)==1 % Calculate Jacobian matrix.
    F= eye(6,6);
end
if isempty(P)==1 % Calculate error state transition matrix.
    P= zeros(6,6);
end

if isempty(Xp)==1 %Define an initial state vector
    Xp = [r0,v0]';
end
%% Predict
%Step 1: Project the state ahead
    %Define the state transition matrix A
    a = 350 + Re; %Average spacecraft position in orbit in km
    A = [1 0 0 dt 0 0;
         0 1 0 0 dt 0;
         0 0 1 0 0 dt;
         -mu./a.^3.*dt 0 0 1 0 0;
         0 -mu./a.^3.*dt 0 0 1 0;
         0 0 -mu./a.^3.*dt 0 0 1];
     
     %Define the control input model
     B = zeros(6,6);

     %Define the control inputs
     up = zeros(6,1);

     %Project the state ahead
     Xp = A*Xp+B*up;
     
 %Step 2: Project the error covariance ahead
     P = A*P*A' + Q;

%% Correct
%Step 1: Compute the Kalman Gain
    K = P*F'*(F*P*F'+R)^-1;
    
%Step 2:Update the estimate with measurement zk
    Xp = Xp + K*([r0,v0]'-F*Xp);
    
%Step 3: Update the error covariance
    P = (eye(6,6) - K*F)*P;
    
%Output results
r1 = Xp(1:3,1);
v1 = Xp(4:6,1);

end