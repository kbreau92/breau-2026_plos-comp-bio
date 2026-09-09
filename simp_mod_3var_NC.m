clear all



k3 = .8; %maximum rate for CELSR-FZS VANGL binding 
k4 = .6; %maximum rate for CELSR-FZS-VANGL dissociation

h3 =  0.2; %concentration of CELSR-FZD-VANGL_i at which FZD-VANGL_j binding to CELSR is half max
h4 = 0.2;  %concentration of CELSR-FZD-VANGL_i at which rate constant for CELSR-FZS-VANGL_i dissociation is half max

l = 2; %Hill coefficient for regulation of CELSR FZD_VANGL binding
ll = 2; %Hill coefficient for regulation of CELSR-FZD-VANGL dissociation

%Total concentrations of CELSR
CT = 1;


tfinal = 50;
dt = .01;
nits = ceil(tfinal/dt);

cv2vals = [0:.005:CT];
[m,n] = size(cv2vals);
cv1out = zeros(1,n);

cv1 = 0.1;

for nn = 1:n  %interate through values of cv2

cv2 = cv2vals(nn);

    for i = 1:nits %integrate cv1 equation to find the steady state 

        c = CT  - cv1 -cv2;
        cv1 = cv1 + dt*(k3/((cv2/h3)^l+1)*c - k4*cv1/(1+(cv1/h4)^ll));
  
    end

    cv1out(1,nn) = cv1;
end


%Plot nullclines
figure(1)
clf 
plot(cv2vals, cv1out)
hold on
plot(cv1out,cv2vals)
hold on 
xlabel('CFV_1')
ylabel('CFV_2')



