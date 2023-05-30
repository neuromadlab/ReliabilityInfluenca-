%% Modeling Influenca Data 2021 - Model Fitting (Basic Model)
% contact: nils.kroemer@uni-tuebingen.de @cornu_copiae
% Parameter Estimation for Basic Model (fmincon call) 

function L = fit_model_basic(x,D)

%Model Parameters, original units 
alpha = x(1); 
beta = x(2);

%um einen Parameter erweitern wegen 2 Learning Rates 
if length(x) > 2
    gamma = x(3);
else
    gamma = 1;
end

%Model Parameters, transformed for unconstrained search
%alpha = 1/(1+exp(-x(1)));
%beta  = exp(x(2)); 
%gamma = exp(x(3));

%Parse Data
rew = D(:,2);
action = D(:,1);
outcome = D(:,3);
rew_grid = D(:,4:5);
%win = D(:,6); %line 60 T.win  

%Initialization of values for the model
L = 0;                          %initial value sum of squared error
RPE_alpha = zeros(length(action),1);  %initial RPEs
Q = zeros(length(action),2);    %initial Q values
est_prob = 0.5 * ones(length(action),1); % initial estimates of probability 

%loop through trials
for t = 1:length(action)   
       
    % Compute outcome prediction error
    % rew. pred. error is difference between was I right (1) or wrong(0) and est. prob
    % for that reward
    RPE_alpha(t,1) = outcome(t,1) - est_prob(t,1);
    %disp(RPE_alpha(t,1));
    % draw_blue immer 0 oder 1 und est_prob für den ersten Durchgang immer
    % 0,5; damit RPE_alpha im ersten Durchgang immer größer 0
 
    %updates Q values
    if t < length(action)          
        %Basic Model 
        est_prob(t+1,1) = est_prob(t,1) + alpha * RPE_alpha(t,1);
        Q(t+1,1) = max(min(gamma*(est_prob(t+1,1)-0.5)+0.5,1),0)*(rew_grid(t+1,1)/50);
        Q(t+1,2) = max(min(gamma*((1-est_prob(t+1,1))-0.5)+0.5,1),0)*(rew_grid(t+1,2)/50);        
        %disp(Q(t+1,2));         
    end    
    %calculates loglikelihood based on discrepancies between observed and
    %predicted choices
    L = L + log(exp(Q(t,action(t,1)) .* beta) / (exp(Q(t,1) .* beta) + exp(Q(t,2) .* beta)));
    %disp(Q(t,action(t,1)));
end

L=-L; %negative due to function minimization
 
end


