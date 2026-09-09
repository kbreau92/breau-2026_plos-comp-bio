clear all



k3 = .8; %maximum rate constant for CELSR-FZS VANGL binding 
k4 = .6; %maximum rate constant for CELSR-FZS-VANGL dissociation

h3 = 0.2; %concentration of complexed FZD at which mutual inhibition is half max
h4 = 0.2;  %concentration of CELSR-FZS-VANGL at which auto-activation is half max

l = 2; %Hill coefficient for regulation of CELSR-FZS VANGL binding
ll = 2; %Hill coefficient for regulation of CELSR-FZS-VANGL dissociation

%Total concentrations of CELSR, VANGL and FZS
CT = 1;

tfinal = 50;
dt = .01;
nits = ceil(tfinal/dt);

%output vectors
tout = [0:dt:tfinal];
cf2out = zeros(1,nits);
cv2out = zeros(1,nits);

%Initial conditions
cv1 = 0.05;
cv2 = 0.0;

cv1out(1,1) = cv1;
cv2out(1,1) = cv2;

for i = 1:nits %iterate equations

    %cocentrations of free CELSR, VANGL and FZS
    c = CT  - cv1 -cv2;
  
    cv1_n = cv1 + dt*(k3/((cv2/h3)^l+1)*c - k4*cv1/(1+(cv1/h4)^ll));
    cv2_n = cv2 + dt*(k3/((cv1/h3)^l+1)*c - k4*cv2/(1+(cv2/h4)^ll));
   
    cv1 = cv1_n;
    cv2 = cv2_n;

    cv1out(1,i+1) = cv1;
    cv2out(1,i+1) = cv2;
end



figure(2)
clf 
plot(tout, cv1out, tout, cv2out)
axis([0 tfinal 0 CT])
xlabel('Time')
ylabel('CFV tetramer')
legend('CFV1','CFV2')


