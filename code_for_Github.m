%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This m file documents the codes used to generate table and predicted
% solvent content for the manuscript. These codes require the presence of
% PLS toolbox v.9.3.1 (Eigenvector Research, Inc. Manson, WA) in order to 
% be fully functional. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snv + CLS w 25°C PC
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_25C)]);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_25C)]);
CLS_pred_API1=snv(APIdata_1.data)*pinv([snv(K_25C(:,135:end))]);
CLS_pred_API2=snv(APIdata_2.data)*pinv([snv(K_25C)]);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% snv + CLS w 60°C PC
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_60C)]);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_60C)]);
CLS_pred_API1=snv(APIdata_1.data)*pinv([snv(K_60C(:,135:end))]);
CLS_pred_API2=snv(APIdata_2.data)*pinv([snv(K_60C)]);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% snv + PACLS(>25°C)
RES=snv(RESIDUAL_CAL)-Y_RESIDUAL*snv(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_API1=snv(APIdata_1.data)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
CLS_pred_API2=snv(APIdata_2.data)*pinv([snv(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% snv + PACLS(>60°C)
RES=snv(RESIDUAL_CAL1)-Y_RESIDUAL1*snv(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_API1=snv(APIdata_1.data)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
CLS_pred_API2=snv(APIdata_2.data)*pinv([snv(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% GLS (γ = 0.7) + snv + PACLS(>25°C)
RES=snv(RESIDUAL_CAL)-Y_RESIDUAL*snv(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);

ind=[0 1e-4 1e-3 1e-2 0.1 0.3 0.5 0.7 0.9 1 1.1 1.3 1.5 10 100];
for i=1:15
    X_OUT{1,i}=GLSW(APIdata_1.data,PCA_diff.loads{2,1},ind(i));
    CLS_pred_API1=snv(X_OUT{1,i})*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
    RMSEP_API1(i)=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
end

X_OUT=GLSW(APIdata_1.data,PCA_diff.loads{2,1},0.7);
CLS_pred_API1=snv(X_OUT)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);

for i=1:15
    X_OUT{1,i}=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},ind(i));
    CLS_pred_API2=snv(X_OUT{1,i})*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
    RMSEP_API2(i)=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
end

X_OUT=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},0.7);
CLS_pred_API2=snv(X_OUT)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% GLS (γ = 0.7) + snv + PACLS(>60°C)
RES=snv(RESIDUAL_CAL1)-Y_RESIDUAL1*snv(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=snv(data_70C.data)*pinv([snv(K_25C);V(:,1:5)']);
CLS_pred_80C=snv(data_80C.data)*pinv([snv(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);

ind=[0 1e-4 1e-3 1e-2 0.1 0.3 0.5 0.7 0.9 1 1.1 1.3 1.5 10 100];
for i=1:15
    X_OUT{1,i}=GLSW(APIdata_1.data,PCA_diff.loads{2,1},ind(i));
    CLS_pred_API1=snv(X_OUT{1,i})*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
    RMSEP_API1(i)=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);
end

X_OUT=GLSW(APIdata_1.data,PCA_diff.loads{2,1},0.7);
CLS_pred_API1=snv(X_OUT)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);

for i=1:15
    X_OUT{1,i}=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},ind(i));
    CLS_pred_API2=snv(X_OUT{1,i})*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
    RMSEP_API2(i)=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
end

X_OUT=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},0.7);
CLS_pred_API2=snv(X_OUT)*pinv([snv(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% GLS (γ = 0.7) + PACLS (>25°C) w/o snv
RES=(RESIDUAL_CAL)-Y_RESIDUAL*(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=(data_70C.data)*pinv([(K_25C);V(:,1:5)']);
CLS_pred_80C=(data_80C.data)*pinv([(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);

X_OUT=GLSW(APIdata_1.data,PCA_diff.loads{2,1},0.7);
CLS_pred_API1=(X_OUT)*pinv([(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);

X_OUT=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},0.7);
CLS_pred_API2=(X_OUT)*pinv([(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% GLS (γ = 0.7) + PACLS (>60°C) w/o snv
RES=(RESIDUAL_CAL1)-Y_RESIDUAL1*(K_25C);
[U,S,V]=svd(RES);
CLS_pred_70C=(data_70C.data)*pinv([(K_25C);V(:,1:5)']);
CLS_pred_80C=(data_80C.data)*pinv([(K_25C);V(:,1:5)']);
RMSEP_70C=sqrt(sum((CLS_pred_70C(ind_70C,1)-Y_CAL_GC_70C).^2)/24);
RMSEP_80C=sqrt(sum((CLS_pred_80C(ind_80C,1)-Y_CAL_GC_80C).^2)/9);

X_OUT=GLSW(APIdata_1.data,PCA_diff.loads{2,1},0.7);
CLS_pred_API1=(X_OUT)*pinv([(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API1=sqrt(sum((CLS_pred_API1(ind_API1,1)-Y_GC_API1).^2)/4);

X_OUT=GLSW(APIdata_2.data(:,135:end),PCA_diff.loads{2,1},0.7);
CLS_pred_API2=(X_OUT)*pinv([(K_25C(:,135:end));V(135:end,1:5)']);
RMSEP_API2=sqrt(sum((CLS_pred_API2(ind_API2,1)-Y_GC_API2).^2)/8);

function X_OUT=GLSW(X_IN,K,Gamma)
% Purpose: remove K by GLS shrinking per Jour of Chemometrics. 2010,24:288-299
% Made by: H.Martens 2003
% Modified by: Pete Shi 2024 
% Input:
%   X_IN(nObj x nXVar) input spectra
%   K (nXVar x nK) spectra of bad components
%   Gamma (scalar) squared degree of removal
%
% Output:
%   X_OUT(nObj x nVar) = X after partial removal of P
%

[nXVar,nK]=size(K);
[nObj,nXVar]=size(X_IN);

[P,S,V]=svd(K,0);

Sigma=P*P'*Gamma + eye(nXVar);

X_OUT=X_IN*(Sigma^(-0.5));
