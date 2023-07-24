% Based on Chengcheng Huang code

function [Corr, Cov, rate2, rate1, var2, var1]=spkCntCorr(s1,s2,N1,N2,dim,Nc) 
% s1, s2: spike times of neuron Pop1 & Pop2, respectively, 2xN (1st row contains spike times, 2nd row contains indices of neurons that spike)
% N1, N2: total # of neurons in Pop1 & Pop2, respectively
% Nc: 1x2 # of sampled neurons, e.g. Nc=[500, 500];
% only neuron in the center square [.25, .75]x[.25, .75] are sampled 
% Corr= spike count correlations
% Cov: same as Corr for covariance 
% rate1, rate2: mean rate (Hz) of the Nc sampled neurons 
% var1, var2: variance of the Nc sampled neurons 

T=floor(min([max(s1(1,:)) max(s2(1,:))]));

s1=s1(:,s1(1,:)<=T);
s2=s2(:,s2(1,:)<=T);

switch dim
    case '1D'
        I1=transpose(unique(s1(2,:))); % sorted indices, I0 from 1 to Np
        I2=transpose(unique(s2(2,:)));
        
        I1=I1(I1<=N1);
        I2=I2(I2<=N2);
        
        Ic1=randsample(I1,Nc(1));
        Ic2=randsample(I2,Nc(2));

    case '2D'
        N11=round(sqrt(N1));
        N21=round(sqrt(N2));
        if size(s1,1)==3
        s1(2,:)=(s1(2,:)-1)*N11+s1(3,:); % x=ceil((I)/Nx1), y=mod(I-1,Nx1)+1
        end
        if size(s2,1)==3
        s2(2,:)=(s2(2,:)-1)*N21+s2(3,:);
        end
        I1=transpose(unique(s1(2,:))); %%% OG: indices of all the neurons which spiked (excitatory and inhibitory)
        I2=transpose(unique(s2(2,:)));
        
        I1=I1(I1<=N1); %%% OG: only select excitatory neurons
        I2=I2(I2<=N2&I2>0);
        
        Ix10=(ceil(I1/N11))/N11;
        Iy10=(mod((I1-1),N11)+1)/N11;
        
        Ix20=(ceil(I2/N21))/N21;
        Iy20=(mod((I2-1),N21)+1)/N21;
        I1=I1(Ix10<0.75 & Ix10>0.25 & Iy10<0.75 & Iy10>0.25);
        I2=I2(Ix20<0.75 & Ix20>0.25 & Iy20<0.75 & Iy20>0.25);
        
        Ic1=randsample(I1,Nc(1)); %%% OG: Get the indices of Nc(1) randomly selected excitatory neurons from L1
        Ic2=randsample(I2,Nc(2)); %%% OG: Get the indices of Nc(2) randomly selected excitatory neurons from L2
end

% compute spike counts using sliding window 
Tw = 100; % sliding window size NB (OG): it is 200ms in corr_d.m
Tburn = 1000; 
time = 0:1:T; %%% OG: 1ms timewindows

%%% OG: compute instantaneous firing rate per 1ms time window
re1=zeros(Nc(1),length(time)); %%% ... of Nc(1) randomly selected excitatory neurons in L1 (re1)
re2=zeros(Nc(2),length(time)); %%% ... of Nc(2) randomly selected excitatory neurons in L2 (re2)

%%% OG: get the number of spikes per 1ms timewindow (loop over each
%%% selected neuron), and then multiply by 1e3 to get a firing rate in [Hz]
for mm=1:Nc(1)
    re1(mm,:)=hist(s1(1,Ic1(mm)-1/4<s1(2,:) & s1(2,:)<=Ic1(mm)+1/4),time)*1e3;
end
for mm=1:Nc(2)
    re2(mm,:)=hist(s2(1,Ic2(mm)-1/4<s2(2,:) & s2(2,:)<=Ic2(mm)+1/4),time)*1e3;
end
re2_s=imfilter(re2(:,Tburn+1:end),ones(1,Tw)/Tw);re2_s=re2_s(:,Tw/2-1:end-Tw/2);
re1_s=imfilter(re1(:,Tburn+1:end),ones(1,Tw)/Tw);re1_s=re1_s(:,Tw/2-1:end-Tw/2);

ind1=mean(re1_s,2)>2;
ind2=mean(re2_s,2)>2;
re1_s=re1_s(ind1,:);
re2_s=re2_s(ind2,:);

switch dim 
    case '1D'
        Ix1=Ic1(ind1)/N1; Nc(1)=length(Ix1);
        Ix2=Ic2(ind2)/N2; Nc(2)=length(Ix2);
        D = pdist2([Ix2;Ix1],[Ix2;Ix1],'euclidean');
        D(D>0.5)=1-D(D>0.5); % periodic boundary condition
    case '2D'
        Ix1=(ceil((Ic1(ind1))/N11))/N11;
        Iy1=(mod((Ic1(ind1)-1),N11)+1)/N11;
        Nc(1)=length(Ix1);
        Ix2=(ceil((Ic2(ind2))/N21))/N21;
        Iy2=(mod((Ic2(ind2)-1),N21)+1)/N21;
        Nc(2)=length(Ix2);
        D = pdist2([[Ix2;Ix1],[Iy2;Iy1]],[[Ix2;Ix1],[Iy2;Iy1]],'euclidean');
end

Cov = cov([re2_s; re1_s]'); % L2 first, then L1
Var = diag(Cov);
var2 = Var(1:Nc(2));
var1 = Var(Nc(2)+1:end);
rate1 = mean(re1_s,2);
rate2 = mean(re2_s,2);

Corr = Cov./sqrt(Var*Var'); % L2 first, then L1

end





