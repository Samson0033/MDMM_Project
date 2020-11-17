
% This code replicates the empirical results of SGU (2018) using 
% data on UEMOA, CEMAC and CMA countries 

% By Samson M'boueke

% Some of the functions used are borrowed from the authors
%%%% CODE NOT FINAL


clear all; close all; clc; dbstop if error; 
rng('default');


[~,sheet_name]=xlsfinfo('data_uemoa_cemac_cma.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('data_uemoa_cemac_cma.xlsx',sheet_name{k});
end

nc = size(data,2); %number of countries

% Taking the logs of the variables, except TB

for j=1:nc
data_new{1,j} = [log(data{1,j}(:,2)) data{1,j}(:,3) log(data{1,j}(:,4)) log(data{1,j}(:,5))...
             log(data{1,j}(:,6)) log(data{1,j}(:,7))];
end

% Quadratic detrending of the variables, as specified in SGU

for j=1:nc
data_new{1,j} = [detrend(data_new{1,j}(:,1),2)... 
                 detrend((data{1,j}(:,3))./(data{1,j}(:,4)-detrend(data{1,j}(:,4),2)),2)... 
                 detrend(data_new{1,j}(:,3),2) detrend(data_new{1,j}(:,4),2) detrend(data_new{1,j}(:,5),2)... 
                 detrend(data_new{1,j}(:,6),2)];
end

nLags=1;
useConstant=1;
useTrend=0;
nv = size(data_new{1,1},2);


% Estimate the reduced form via OLS for
% countries with no missing observations
cc = [2:3 6:8 9 14 17:18];

% for j=cc
% reduced_form(j) = estimatevarols(data_new{1,j},nLags,useConstant,useTrend);

% Recover structural parameters
% end


% Estimate the structural form via OLS for
% countries with no missing observations

R = ones(nv);
R(1,2:end) = 0; %restriction on A (tot univariate)
cholesky = 1;
T = 11; % length of impulse responses


for k=cc
%construct current and lagged data for country k
d = data_new{k};
y0(k,:) = d(1,:);
 [A(:,:,k),PI(:,:,k),R2(:,k),~,~,stdu(:,k),U{k}] = svar_estim(d,R,cholesky);

% variance decomposition
[~,vdx]=variance_decomposition(eye(nv),A(:,:,k),PI(:,:,k));
v_share(k,:) = vdx(1,:)*100;

%Population variance
[~,sigx(:,:,k)]=mom(eye(nv),A(:,:,k),PI(:,:,k)*PI(:,:,k)');
PI1 = PI(:,:,k);
PI1(:,2:end) = 0;
[~,ssigx(:,:,k)]=mom(eye(nv),A(:,:,k),PI1*PI1');
var_svar_tot(:,k) = diag(ssigx(:,:,k)); %variances conditional on ttot shocks implied by empirical SVAR model
var_svar(:,k) = diag(sigx(:,:,k)); %total variances  implied by empirical SVAR model

%Impulse Responses
x0 = PI(:,1,k)/PI(1,1,k);
x =ir(eye(nv),A(:,:,k),x0,T);
IR(:,:,k) = x(:,1:nv);

end

