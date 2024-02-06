clear all
clc
close all


%% LOAD FILE
%% Directory Selection

selectedDirectory = uigetdir();
prompt={'Type file name without extension','Type file name to process','Velocity trial','Patient name','Patology'}; 
titles='Files to process';
answer=inputdlg(prompt,titles);
name = answer{1}; %file name without extension: TRIAL
A = str2num(answer{2});%file name to process: number of trial
Velocity=answer{3};%Velocity trial: M or L for NORMATIVO and M for ATASSIA or PARAPARESI
ID=answer{4};%ID: patient's name
Patology=answer{5};%Patology: NORMATIVO, ATASSIA, PARAPARESI

%Variable initialization

index_near_right_2=zeros(length(A),1);
index_near_right_1=zeros(length(A),1);
index_near_right_single_stance=zeros(length(A),1);
index_near_right_toe_off=zeros(length(A),1);
index_near_left_toe_off=zeros(length(A),1);

time_index=cell(length(A),1);
time_vector=cell(length(A),1); 

time_foot_strike_1= zeros(length(A),1);
time_foot_strike_2=zeros(length(A),1);
time_single_stance=zeros(length(A),1);
time_right_toe_off=zeros(length(A),1);
time_left_toe_off=zeros(length(A),1);

R_RF=cell(length(A),1);
normalized_R_RF=cell(length(A),1);

R_VL=cell(length(A),1);
normalized_R_VL=cell(length(A),1);

R_VM=cell(length(A),1);
normalized_R_VM=cell(length(A),1);

R_G=cell(length(A),1);
normalized_R_G=cell(length(A),1);

R_BFCL=cell(length(A),1);
normalized_R_BFCL=cell(length(A),1);

R_TA=cell(length(A),1);
normalized_R_TA=cell(length(A),1);

R_PL=cell(length(A),1);
normalized_R_PL=cell(length(A),1);

R_SO=cell(length(A),1);
normalized_R_SO=cell(length(A),1);

R_GL=cell(length(A),1);
normalized_R_GL=cell(length(A),1);

R_GMA=cell(length(A),1);
normalized_R_GMA=cell(length(A),1);

num_samples=zeros(length(A),1); 
time_window_length=zeros(length(A),1);
new_time_vector=cell(length(A),1);

for i = 1:length(A)

    for k=A(i)
        iesimo=k;
    
    if strcmp(Velocity,"M") || isempty(Velocity) 
    nametdf=strcat(name,num2str(iesimo),'M.tdf');
    nametxt=strcat(name,num2str(iesimo),'M_txt');
    elseif strcmp(Velocity,"L" )
        nametdf=strcat(name,num2str(iesimo),'L.tdf');
        nametxt=strcat(name,num2str(iesimo),'L_txt');
    end
    end

%% import file TDF
[startTime,frequency,emgMap,labels,emgData] = tdfReadDataEmg (nametdf);


%% import file TXT
filename =nametxt;

    data=importdata(filename);

        if isstruct(data)
            
            dati = data.data; 
        else
            
            dati = data;
        end

%% SIGNALS PROCESSING: parameters 

 dt = 1/frequency; %Sampling time
 time_acq=(0:dt:(length(emgData(1,:))-1)*dt); %Acquisition time

%Filter parameters
order=4;%fourth order bandpass ButterWorth filter
order_low=2;%second order low pass ButterWorth filter
low_pass_cutoff=10;
high_pass_cutoff=400;
nyq_f=frequency/2;
bpWn(1)=low_pass_cutoff/nyq_f;
bpWn(2)=high_pass_cutoff/nyq_f;

%% Definition of time istants of gait cycle

%Gait Cycle events
right_foot_strike_1=data.data(5,1); 
left_toe_off=data.data(6,1);
right_single_stance=data.data(7,1);
right_toe_off=data.data(8,1);
right_foot_strike_2=data.data(9,1); 

%Calculation of the difference between instants of time
diff_right_1 = abs(time_acq - right_foot_strike_1);
diff_right_2 = abs(time_acq - right_foot_strike_2);
diff_right_single_stance = abs(time_acq - right_single_stance);
diff_left_toe_off = abs(time_acq - left_toe_off);
diff_right_toe_off= abs(time_acq - right_toe_off);

%Instant index with the smallest difference
index_near_right_1(i,1) = find(diff_right_1 == min(diff_right_1));
index_near_right_2(i,1) = find(diff_right_2 == min(diff_right_2));
index_near_right_single_stance(i,1) = find(diff_right_single_stance == min(diff_right_single_stance));
index_near_left_toe_off(i,1) = find(diff_left_toe_off == min(diff_left_toe_off));
index_near_right_toe_off(i,1) = find(diff_right_toe_off == min(diff_right_toe_off));

%Nearest instant to index found
time_foot_strike_1 (i,1)= time_acq(index_near_right_1(i,1));
time_foot_strike_2 (i,1)= time_acq(index_near_right_2(i,1));
time_single_stance (i,1)= time_acq(index_near_right_single_stance(i,1));
time_left_toe_off (i,1)= time_acq(index_near_left_toe_off(i,1));
time_right_toe_off (i,1)= time_acq(index_near_right_toe_off(i,1));

time_index{i}=index_near_right_1(i,1):index_near_right_2(i,1);
time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1),index_near_right_2(i,1)-index_near_right_1(i,1) +1);

%% SIGNALS PROCESSING: 

%% RIGHT RECTUS FEMORIS

%Raw signal
R_RF_raw=emgData(1,:);
%Centered signal
offset_R_RF=mean(emgData(1,:));
R_RF_cen=emgData(1,:)-(offset_R_RF);
%Band-pass filter
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_RF_filt=filtfilt(bpB,bpA,R_RF_cen);
%Rectified signal
R_RF_abs=abs(R_RF_filt);
%Envelope
cutoff_env=5/nyq_f; %5Hz cutoff frequency; cutoff normalized with respect to Nyquist frequency
[enB,enA]=butter(order_low,cutoff_env);
R_RF_env=filtfilt(enB,enA,R_RF_abs);

R_RF{i}=R_RF_env(time_index{i});


%% VASTUS LATERALIS

R_VL_raw=emgData(2,:);
offset_R_VL=mean(emgData(2,:));
R_VL_cen=emgData(2,:)-(offset_R_VL);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_VL_filt=filtfilt(bpB,bpA,R_VL_cen);
R_VL_abs=abs(R_VL_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_VL_env=filtfilt(enB,enA,R_VL_abs);

R_VL{i}=R_VL_env(time_index{i});


%% VASTUS MEDIALIS

R_VM_raw=emgData(3,:);
offset_R_VM=mean(emgData(3,:));
R_VM_cen=emgData(3,:)-(offset_R_VM);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_VM_filt=filtfilt(bpB,bpA,R_VM_cen);
R_VM_abs=abs(R_VM_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_VM_env=filtfilt(enB,enA,R_VM_abs);

R_VM{i}=R_VM_env(time_index{i});

%% GLUTEUS MEDIUS

R_G_raw=emgData(4,:);
offset_R_G=mean(emgData(4,:));
R_G_cen=emgData(4,:)-(offset_R_G);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_G_filt=filtfilt(bpB,bpA,R_G_cen);
R_G_abs=abs(R_G_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_G_env=filtfilt(enB,enA,R_G_abs);

R_G{i}=R_G_env(time_index{i});

%% BICEPS FEMORIS CAPUT LONGUS

R_BFCL_raw=emgData(7,:);
offset_R_BFCL=mean(emgData(7,:));
R_BFCL_cen=emgData(7,:)-(offset_R_BFCL);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_BFCL_filt=filtfilt(bpB,bpA,R_BFCL_cen);
R_BFCL_abs=abs(R_BFCL_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_BFCL_env=filtfilt(enB,enA,R_BFCL_abs);

R_BFCL{i}=R_BFCL_env(time_index{i});


%% TIBIALIS ANTERIOR

R_TA_raw=emgData(8,:);
offset_R_TA=mean(emgData(8,:));
R_TA_cen=emgData(8,:)-(offset_R_TA);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_TA_filt=filtfilt(bpB,bpA,R_TA_cen);
R_TA_abs=abs(R_TA_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_TA_env=filtfilt(enB,enA,R_TA_abs);

R_TA{i}=R_TA_env(time_index{i});


%% PERONEUS LONGUS

R_PL_raw=emgData(9,:);
offset_R_PL=mean(emgData(9,:));
R_PL_cen=emgData(9,:)-(offset_R_PL);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_PL_filt=filtfilt(bpB,bpA,R_PL_cen);
R_PL_abs=abs(R_PL_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_PL_env=filtfilt(enB,enA,R_PL_abs);

R_PL{i}=R_PL_env(time_index{i});


%% SOLEUS

R_SO_raw=emgData(10,:);
offset_R_SO=mean(emgData(10,:));
R_SO_cen=emgData(10,:)-(offset_R_SO);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_SO_filt=filtfilt(bpB,bpA,R_SO_cen);
R_SO_abs=abs(R_SO_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_SO_env=filtfilt(enB,enA,R_SO_abs);

R_SO{i}=R_SO_env(time_index{i});

%% GASTROCNEMIUS LATERALIS

R_GL_raw=emgData(12,:);
offset_R_GL=mean(emgData(12,:));
R_GL_cen=emgData(12,:)-(offset_R_GL);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_GL_filt=filtfilt(bpB,bpA,R_GL_cen);
R_GL_abs=abs(R_GL_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_GL_env=filtfilt(enB,enA,R_GL_abs);

R_GL{i}=R_GL_env(time_index{i});


%% GLUTEUS MAXIMUS

R_GMA_raw=emgData(13,:);
offset_R_GMA=mean(emgData(13,:));
R_GMA_cen=emgData(13,:)-(offset_R_GMA);
[bpB,bpA]=butter(order,bpWn,'bandpass');
R_GMA_filt=filtfilt(bpB,bpA,R_GMA_cen);
R_GMA_abs=abs(R_GMA_filt);
cutoff_env=5/nyq_f; 
[enB,enA]=butter(order_low,cutoff_env);
R_GMA_env=filtfilt(enB,enA,R_GMA_abs);

R_GMA{i}=R_GMA_env(time_index{i});


%% SAMPLE NUMBER CALCULATION

time_window_length(i,1)=(time_foot_strike_2(i,1)-time_foot_strike_1(i,1));
num_samples(i,1)=time_window_length(i,1)*frequency;

end


%% 

max_percentage=0.2;

%% RIGHT RECTUS FEMORIS

% Normalization and threshold
max_global_RF = -Inf;  

for i = 1:length(A)
    max_local_RF = max(R_RF{i});
    if max_local_RF > max_global_RF
        max_global_RF = max_local_RF;
    end
end

for i=1:length(A)
normalized_R_RF{i}=(R_RF{i})./(max_global_RF);
end

threshold_R_RF=max_percentage*max_global_RF;
threshold_R_RF_norm = (threshold_R_RF) / (max_global_RF);



% Signal above threshold
for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_RF{i}(j)>threshold_R_RF_norm
            R_RF_values{i,j}=normalized_R_RF{i}(j);
        else
            R_RF_values{i,j}=0;
    
        end
    end
end

R_RF_above_threshold = cell(length(A), 1);
non_empty_cells_RF=~cellfun('isempty',R_RF_values);

for i = 1:length(A)
    non_empty_values_RF=R_RF_values(i,non_empty_cells_RF(i,:));    
    non_empty_values_RF_double = zeros(1, numel(non_empty_values_RF)); 
    
    for k = 1:length(non_empty_values_RF)
        if isnumeric(non_empty_values_RF{k})
            non_empty_values_RF_double(k) = double(non_empty_values_RF{k});
        end
    end

    R_RF_above_threshold{i} = non_empty_values_RF_double;
end


%% VASTUS LATERALIS
% Normalization and threshold


max_global_VL = -Inf;  

for i = 1:length(A)
    max_local_VL = max(R_VL{i});
    if max_local_VL > max_global_VL
        max_global_VL = max_local_VL;
    end
end

for i=1:length(A)
normalized_R_VL{i}=(R_VL{i})./(max_global_VL);
end

threshold_R_VL=max_percentage*max_global_VL;
threshold_R_VL_norm = (threshold_R_VL ) / (max_global_VL);

%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_VL{i}(j)>threshold_R_VL_norm
            R_VL_values{i,j}=normalized_R_VL{i}(j);
        else
            R_VL_values{i,j}=0;
        end
    end
end

R_VL_above_threshold=cell(length(A),1);
non_empty_cells_VL=~cellfun('isempty',R_VL_values);

for i = 1:length(A)
    non_empty_values_VL=R_VL_values(i,non_empty_cells_VL(i,:));    
    non_empty_values_VL_double = zeros(1, numel(non_empty_values_VL)); 
    
    for k = 1:length(non_empty_values_VL)
        if isnumeric(non_empty_values_VL{k})
            non_empty_values_VL_double(k) = double(non_empty_values_VL{k});
        end
    end

    R_VL_above_threshold{i} = non_empty_values_VL_double;
end


%% VASTUS MEDIALIS
% Normalization and threshold


max_global_VM = -Inf;  

for i = 1:length(A)
    max_local_VM = max(R_VM{i});
    if max_local_VM > max_global_VM
        max_global_VM = max_local_VM;
    end
end

for i=1:length(A)
normalized_R_VM{i}=(R_VM{i})./(max_global_VM);
end

threshold_R_VM=max_percentage*max_global_VM;
threshold_R_VM_norm = (threshold_R_VM) / (max_global_VM);


%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_VM{i}(j)>threshold_R_VM_norm
            R_VM_values{i,j}=normalized_R_VM{i}(j);
        else
            R_VM_values{i,j}=0;
        end
    end
end

R_VM_above_threshold=cell(length(A),1);
non_empty_cells_VM=~cellfun('isempty',R_VM_values);

for i = 1:length(A)
    non_empty_values_VM=R_VM_values(i,non_empty_cells_VM(i,:));    
    non_empty_values_VM_double = zeros(1, numel(non_empty_values_VM)); 
    
    for k = 1:length(non_empty_values_VM)
        if isnumeric(non_empty_values_VM{k})
            non_empty_values_VM_double(k) = double(non_empty_values_VM{k});
        end
    end

    R_VM_above_threshold{i} = non_empty_values_VM_double;
end



%% GLUTEUS MEDIUS

% Normalization and threshold

max_global_G = -Inf;  

for i = 1:length(A)
    max_local_G = max(R_G{i});
    if max_local_G > max_global_G
        max_global_G = max_local_G;
    end
end

for i=1:length(A)
normalized_R_G{i}=(R_G{i})./(max_global_G);
end

threshold_R_G=max_percentage*max_global_G;
threshold_R_G_norm = (threshold_R_G) / (max_global_G);

%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_G{i}(j)>threshold_R_G_norm
            R_G_values{i,j}=normalized_R_G{i}(j);
        else
            R_G_values{i,j}=0;
        end
    end
end

R_G_above_threshold=cell(length(A),1);
non_empty_cells_G=~cellfun('isempty',R_G_values);

for i = 1:length(A)
    non_empty_values_G=R_G_values(i,non_empty_cells_G(i,:));    
    non_empty_values_G_double = zeros(1, numel(non_empty_values_G)); 
    
    for k = 1:length(non_empty_values_G)
        if isnumeric(non_empty_values_G{k})
            non_empty_values_G_double(k) = double(non_empty_values_G{k});
        end
    end

    R_G_above_threshold{i} = non_empty_values_G_double;
end


%% BICEPS FEMORIS CAPUT LONGUS

%Normalization and threshold

max_global_BFCL = -Inf;  

for i = 1:length(A)
    max_local_BFCL = max(R_BFCL{i});
    if max_local_BFCL > max_global_BFCL
        max_global_BFCL = max_local_BFCL;
    end
end

for i=1:length(A)
normalized_R_BFCL{i}=(R_BFCL{i})./(max_global_BFCL);
end

threshold_R_BFCL=max_percentage*max_global_BFCL;
threshold_R_BFCL_norm = (threshold_R_BFCL) / (max_global_BFCL);


%signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_BFCL{i}(j)>threshold_R_BFCL_norm
            R_BFCL_values{i,j}=normalized_R_BFCL{i}(j);
        else
            R_BFCL_values{i,j}=0;
        end
    end
end

R_BFCL_above_threshold=cell(length(A),1);
non_empty_cells_BFCL=~cellfun('isempty',R_BFCL_values);

for i = 1:length(A)
    non_empty_values_BFCL=R_BFCL_values(i,non_empty_cells_BFCL(i,:));    
    non_empty_values_BFCL_double = zeros(1, numel(non_empty_values_BFCL)); 
    
    for k = 1:length(non_empty_values_BFCL)
        if isnumeric(non_empty_values_BFCL{k})
            non_empty_values_BFCL_double(k) = double(non_empty_values_BFCL{k});
        end
    end

    R_BFCL_above_threshold{i} = non_empty_values_BFCL_double;
end


%% TIBIALIS ANTERIOR

% Normalization and threshold

max_global_TA = -Inf;  

for i = 1:length(A)
    max_local_TA = max(R_TA{i});
    if max_local_TA > max_global_TA
        max_global_TA = max_local_TA;
    end
end

for i=1:length(A)
normalized_R_TA{i}=(R_TA{i})./(max_global_TA);
end

threshold_R_TA=max_percentage*max_global_TA;
threshold_R_TA_norm = (threshold_R_TA) / (max_global_TA);


%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_TA{i}(j)>threshold_R_TA_norm
            R_TA_values{i,j}=normalized_R_TA{i}(j);
        else
            R_TA_values{i,j}=0;
        end
    end
end

R_TA_above_threshold=cell(length(A),1);
non_empty_cells_TA=~cellfun('isempty',R_TA_values);

for i = 1:length(A)
    non_empty_values_TA=R_TA_values(i,non_empty_cells_TA(i,:));    
    non_empty_values_TA_double = zeros(1, numel(non_empty_values_TA)); 
    
    for k = 1:length(non_empty_values_TA)
        if isnumeric(non_empty_values_TA{k})
            non_empty_values_TA_double(k) = double(non_empty_values_TA{k});
        end
    end

    R_TA_above_threshold{i} = non_empty_values_TA_double;
end


%% PERONEUS LONGUS

% Normalization and threshold

max_global_PL = -Inf;  

for i = 1:length(A)
    max_local_PL = max(R_PL{i});
    if max_local_PL > max_global_PL
        max_global_PL = max_local_PL;
    end
end

for i=1:length(A)
normalized_R_PL{i}=(R_PL{i})./(max_global_PL);
end

threshold_R_PL=max_percentage*max_global_PL;
threshold_R_PL_norm = (threshold_R_PL) / (max_global_PL);

%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_PL{i}(j)>threshold_R_PL_norm
            R_PL_values{i,j}=normalized_R_PL{i}(j);
        else
            R_PL_values{i,j}=0;
        end
    end
end

R_PL_above_threshold=cell(length(A),1);
non_empty_cells_PL=~cellfun('isempty',R_PL_values);

for i = 1:length(A)
    non_empty_values_PL=R_PL_values(i,non_empty_cells_PL(i,:));    
    non_empty_values_PL_double = zeros(1, numel(non_empty_values_PL)); 
    
    for k = 1:length(non_empty_values_PL)
        if isnumeric(non_empty_values_PL{k})
            non_empty_values_PL_double(k) = double(non_empty_values_PL{k});
        end
    end

    R_PL_above_threshold{i} = non_empty_values_PL_double;
end



%% SOLEUS

% Normalization and threshold

max_global_SO= -Inf;  

for i = 1:length(A)
    max_local_SO = max(R_SO{i});
    if max_local_SO > max_global_SO
        max_global_SO = max_local_SO;
    end
end

for i=1:length(A)
normalized_R_SO{i}=(R_SO{i})./(max_global_SO);
end

threshold_R_SO=max_percentage*max_global_SO;
threshold_R_SO_norm = (threshold_R_SO) / (max_global_SO);



%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_SO{i}(j)>threshold_R_SO_norm
            R_SO_values{i,j}=normalized_R_SO{i}(j);
        else
            R_SO_values{i,j}=0;
        end
    end
end

R_SO_above_threshold=cell(length(A),1);
non_empty_cells_SO=~cellfun('isempty',R_SO_values);

for i = 1:length(A)
    non_empty_values_SO=R_SO_values(i,non_empty_cells_SO(i,:));    
    non_empty_values_SO_double = zeros(1, numel(non_empty_values_SO)); 
    
    for k = 1:length(non_empty_values_SO)
        if isnumeric(non_empty_values_SO{k})
            non_empty_values_SO_double(k) = double(non_empty_values_SO{k});
        end
    end

    R_SO_above_threshold{i} = non_empty_values_SO_double;
end


%% GASTROCNEMIUS LATERALIS

% Normalization and threshold

max_global_GL = -Inf;  

for i = 1:length(A)
    max_local_GL = max(R_GL{i});
    if max_local_GL > max_global_GL
        max_global_GL = max_local_GL;
    end
end

for i=1:length(A)
normalized_R_GL{i}=(R_GL{i})./(max_global_GL);
end

threshold_R_GL=max_percentage*max_global_GL;
threshold_R_GL_norm = (threshold_R_GL) / (max_global_GL);



%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_GL{i}(j)>threshold_R_GL_norm
            R_GL_values{i,j}=normalized_R_GL{i}(j);
        else
            R_GL_values{i,j}=0;
        end
    end
end

R_GL_above_threshold=cell(length(A),1);
non_empty_cells_GL=~cellfun('isempty',R_GL_values);

for i = 1:length(A)
    non_empty_values_GL=R_GL_values(i,non_empty_cells_GL(i,:));    
    non_empty_values_GL_double = zeros(1, numel(non_empty_values_GL)); 
    
    for k = 1:length(non_empty_values_GL)
        if isnumeric(non_empty_values_GL{k})
            non_empty_values_GL_double(k) = double(non_empty_values_GL{k});
        end
    end

    R_GL_above_threshold{i} = non_empty_values_GL_double;
end

%% GLUTEUS MAXIMUS

%Normalization and threshold

max_global_GMA = -Inf;  

for i = 1:length(A)
    max_local_GMA = max(R_GMA{i});
    if max_local_GMA > max_global_GMA
        max_global_GMA = max_local_GMA;
    end
end

for i=1:length(A)
normalized_R_GMA{i}=(R_GMA{i})./(max_global_GMA);
end

threshold_R_GMA=max_percentage*max_global_GMA;
threshold_R_GMA_norm = (threshold_R_GMA ) / (max_global_GMA);

%Signal above threshold

for i=1:length(A)
    for j=1:length(time_index{i})
        if normalized_R_GMA{i}(j)>threshold_R_GMA_norm
            R_GMA_values{i,j}=normalized_R_GMA{i}(j);
        else
            R_GMA_values{i,j}=0;
        end
    end
end

R_GMA_above_threshold=cell(length(A),1);
non_empty_cells_GMA=~cellfun('isempty',R_GMA_values);

for i = 1:length(A)
    non_empty_values_GMA=R_GMA_values(i,non_empty_cells_GMA(i,:));    
    non_empty_values_GMA_double = zeros(1, numel(non_empty_values_GMA)); 
    
    for k = 1:length(non_empty_values_GMA)
        if isnumeric(non_empty_values_GMA{k})
            non_empty_values_GMA_double(k) = double(non_empty_values_GMA{k});
        end
    end

    R_GMA_above_threshold{i} = non_empty_values_GMA_double;
end

%% INTERPOLATION

%Mean of samples of all trials to be used as uniform value of samples for interpolation
mean_n_samples=mean(num_samples(:,1)); 

%Variable initialitation

R_RF_unif=cell(length(A),1); 
R_RF_unif_uniformed=cell(length(A),1); 
R_VL_unif=cell(length(A),1); 
R_VL_unif_uniformed=cell(length(A),1); 
R_VM_unif=cell(length(A),1); 
R_VM_unif_uniformed=cell(length(A),1); 
R_G_unif=cell(length(A),1); 
R_G_unif_uniformed=cell(length(A),1); 
R_BFCL_unif=cell(length(A),1); 
R_BFCL_unif_uniformed=cell(length(A),1); 
R_TA_unif=cell(length(A),1); 
R_TA_unif_uniformed=cell(length(A),1); 
R_PL_unif=cell(length(A),1); 
R_PL_unif_uniformed=cell(length(A),1); 
R_SO_unif=cell(length(A),1); 
R_SO_unif_uniformed=cell(length(A),1); 
R_GL_unif=cell(length(A),1); 
R_GL_unif_uniformed=cell(length(A),1); 
R_GMA_unif=cell(length(A),1); 
R_GMA_unif_uniformed=cell(length(A),1); 


coef=cell(length(A),1);
xval=cell(length(A),1);
yval=cell(length(A),1);


%Interpolation:RF
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_RF_above_threshold{i},5); 
    R_RF_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_RF_unif_uniformed{i}=yval{i};
    R_RF_unif_uniformed{i} = max(0, R_RF_unif_uniformed{i});
    R_RF_unif_uniformed{i} = min(1, R_RF_unif_uniformed{i});
    

end


%Vector column with RF results 
results_RF=[];
for j=1:numel(R_RF_unif_uniformed)
    vector_RF=R_RF_unif_uniformed{j}';
    results_RF=[results_RF;vector_RF];
end

%Interpolation:VL
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_VL_above_threshold{i},5); 
    R_VL_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_VL_unif_uniformed{i}=yval{i};
    R_VL_unif_uniformed{i} = max(0, R_VL_unif_uniformed{i});
    R_VL_unif_uniformed{i} = min(1, R_VL_unif_uniformed{i});
   

end

%Vector column with VL results
results_VL=[];
for j=1:numel(R_VL_unif_uniformed)
    vector_VL=R_VL_unif_uniformed{j}';
    results_VL=[results_VL;vector_VL];
end

%Interpolation:VM
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_VM_above_threshold{i},5);  
    R_VM_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_VM_unif_uniformed{i}=yval{i};
    R_VM_unif_uniformed{i} = max(0, R_VM_unif_uniformed{i});
    R_VM_unif_uniformed{i} = min(1, R_VM_unif_uniformed{i});
  
end

%Vector column with VM results
results_VM=[];
for j=1:numel(R_VM_unif_uniformed)
    vector_VM=R_VM_unif_uniformed{j}';
    results_VM=[results_VM;vector_VM];
end

%Interpolation:G
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_G_above_threshold{i},5);  
    R_G_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_G_unif_uniformed{i}=yval{i};
    R_G_unif_uniformed{i} = max(0, R_G_unif_uniformed{i});
    R_G_unif_uniformed{i} = min(1, R_G_unif_uniformed{i});
   

end

%Vector column with G results
results_G=[];
for j=1:numel(R_G_unif_uniformed)
    vector_G=R_G_unif_uniformed{j}';
    results_G=[results_G;vector_G];
end

%Interpolation:BFCL
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_BFCL_above_threshold{i},5);  
    R_BFCL_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_BFCL_unif_uniformed{i}=yval{i};
    R_BFCL_unif_uniformed{i} = max(0, R_BFCL_unif_uniformed{i});
    R_BFCL_unif_uniformed{i} = min(1, R_BFCL_unif_uniformed{i});
   

end

%Vector column with BFCL results
results_BFCL=[];
for j=1:numel(R_BFCL_unif_uniformed)
    vector_BFCL=R_BFCL_unif_uniformed{j}';
    results_BFCL=[results_BFCL;vector_BFCL];
end

%Interpolation:TA
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_TA_above_threshold{i},5);  
    R_TA_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_TA_unif_uniformed{i}=yval{i};
    R_TA_unif_uniformed{i} = max(0, R_TA_unif_uniformed{i});
    R_TA_unif_uniformed{i} = min(1, R_TA_unif_uniformed{i});
   

end

%Vector column with TA results
results_TA=[];
for j=1:numel(R_TA_unif_uniformed)
    vector_TA=R_TA_unif_uniformed{j}';
    results_TA=[results_TA;vector_TA];
end

%Interpolation:PL
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_PL_above_threshold{i},5);  
    R_PL_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_PL_unif_uniformed{i}=yval{i};
    R_PL_unif_uniformed{i} = max(0, R_PL_unif_uniformed{i});
    R_PL_unif_uniformed{i} = min(1, R_PL_unif_uniformed{i});
    

end

%Vector column with PL results
results_PL=[];
for j=1:numel(R_PL_unif_uniformed)
    vector_PL=R_PL_unif_uniformed{j}';
    results_PL=[results_PL;vector_PL];
end

%Interpolation SO
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_SO_above_threshold{i},5);  
    R_SO_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_SO_unif_uniformed{i}=yval{i};
    R_SO_unif_uniformed{i} = max(0, R_SO_unif_uniformed{i});
    R_SO_unif_uniformed{i} = min(1, R_SO_unif_uniformed{i});
   
end

%Vector column with SO results
results_SO=[];
for j=1:numel(R_SO_unif_uniformed)
    vector_SO=R_SO_unif_uniformed{j}';
    results_SO=[results_SO;vector_SO];
end

%Interpolation:GL
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_GL_above_threshold{i},5);  
    R_GL_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_GL_unif_uniformed{i}=yval{i};
    R_GL_unif_uniformed{i} = max(0, R_GL_unif_uniformed{i});
    R_GL_unif_uniformed{i} = min(1, R_GL_unif_uniformed{i});
   
end

%Vector column with GL results
results_GL=[];
for j=1:numel(R_GL_unif_uniformed)
    vector_GL=R_GL_unif_uniformed{j}';
    results_GL=[results_GL;vector_GL];
end

%Interpolation:GMA
for i=1:length(A)
    new_time_vector{i}=linspace(time_foot_strike_1(i,1),time_foot_strike_2(i,1) , mean_n_samples); 
    coef{i}=polyfit(time_vector{i},R_GMA_above_threshold{i},5);  
    R_GMA_unif{i}=poly2sym(coef{i});
    xval{i}=new_time_vector{i};
    yval{i}=polyval(coef{i},xval{i});
    R_GMA_unif_uniformed{i}=yval{i};
    R_GMA_unif_uniformed{i} = max(0, R_GMA_unif_uniformed{i});
    R_GMA_unif_uniformed{i} = min(1, R_GMA_unif_uniformed{i});
   

end

%Vector column with GMA results
results_GMA=[];
for j=1:numel(R_GMA_unif_uniformed)
    vector_GMA=R_GMA_unif_uniformed{j}';
    results_GMA=[results_GMA;vector_GMA];
end


%% COLUMN ON-OFF

% Column ON-OFF RF
for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_RF_unif_uniformed{i}(j) > threshold_R_RF_norm
            R_RF_tot{i,j}='ON';
        else
            R_RF_tot{i,j}='OFF';
        end
    end
end

R_RF_transpose=transpose(R_RF_tot);
non_empty_cells=~cellfun('isempty',R_RF_transpose);
R_RF_TOT=R_RF_transpose(non_empty_cells); %Final column for the database

%Column ON-OFF VL

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_VL_unif_uniformed{i}(j) > threshold_R_VL_norm
            R_VL_tot{i,j}='ON';
        else
            R_VL_tot{i,j}='OFF';
        end
    end
end



%Column ON-OFF VM

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_VM_unif_uniformed{i}(j)>threshold_R_VM_norm
            R_VM_tot{i,j}='ON';
        else
            R_VM_tot{i,j}='OFF';
        end
    end
end



%Column ON-OFF G

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_G_unif_uniformed{i}(j)>threshold_R_G_norm
            R_G_tot{i,j}='ON';
        else
            R_G_tot{i,j}='OFF';
        end
    end
end


%Column ON-OFF

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_BFCL_unif_uniformed{i}(j)>threshold_R_BFCL_norm
            R_BFCL_tot{i,j}='ON';
        else
            R_BFCL_tot{i,j}='OFF';
        end
    end
end


%Column ON-OFF TA

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_TA_unif_uniformed{i}(j)>threshold_R_TA_norm
            R_TA_tot{i,j}='ON';
        else
            R_TA_tot{i,j}='OFF';
        end
    end
end

%Column ON-OFF PL

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_PL_unif_uniformed{i}(j)>threshold_R_PL_norm
            R_PL_tot{i,j}='ON';
        else
            R_PL_tot{i,j}='OFF';
        end
    end
end


%Column ON-OFF SO

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_SO_unif_uniformed{i}(j)>threshold_R_SO_norm
            R_SO_tot{i,j}='ON';
        else
            R_SO_tot{i,j}='OFF';
        end
    end
end


%Column ON-OFF GL

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_GL_unif_uniformed{i}(j)>threshold_R_GL_norm
            R_GL_tot{i,j}='ON';
        else
            R_GL_tot{i,j}='OFF';
        end
    end
end



%Column ON-OFF GMA

for i=1:length(A)
    for j=1:floor(mean_n_samples)
        if R_GMA_unif_uniformed{i}(j)>threshold_R_GMA_norm
            R_GMA_tot{i,j}='ON';
        else
            R_GMA_tot{i,j}='OFF';
        end
    end
end



%ID column: second column of the final matrix

results_ID=strings(12,1);
for j=1:12
    results_ID{j}=ID;
end

%Walk column:third column of the final matrix
results_Velocity=strings(12,1);
for j=1:12
    results_Velocity{j}=Velocity;
end

if strcmp(Patology,"NORMATIVO")

 values_WALK=[1,1,1,1,2,2,2,2,3,3,3,3];
 results_WALK=strcat(string(values_WALK(:)),results_Velocity);

elseif ~strcmp(Patology, "NORMATIVO" ) %when patology is ATASSIA or PARAPARESI, Velocity is M
 values_WALK=[1,1,1,1,2,2,2,2,3,3,3,3];
 results_WALK=strcat(string(values_WALK(:)),results_Velocity);

end


%Patology column
results_Patology=strings(12,1);
for j=1:12
    results_Patology{j}=Patology; 
end

%Percentage muscle's activation column

phase_names={'HEEL CONTACT','SINGLE STANCE','TERMINAL STANCE','SWING'};
phases=cell(length(A),length(phase_names));

for i=1:length(A)

for j=1:length(new_time_vector{i})

    if new_time_vector{i}(j)< time_left_toe_off(i)
        phases{i,1}{end+1}=new_time_vector{i}(j);
    elseif new_time_vector{i}(j) >= time_left_toe_off(i)
        if new_time_vector{i}(j)< time_single_stance(i)
             phases{i,2}{end+1}=new_time_vector{i}(j);
             elseif new_time_vector{i}(j) >= time_single_stance(i)
                 if new_time_vector{i}(j)<= time_right_toe_off(i)
                      phases{i,3}{end+1}=new_time_vector{i}(j);
                      elseif new_time_vector{i}(j)> time_right_toe_off(i)
                         phases{i,4}{end+1}=new_time_vector{i}(j);

                end
        end
    end
end
end


numPhases = size(phases, 2);
for i=1:length(A)

    for j=1:numel(phases{i,1})
    results_heel_contact{i,j}=phase_names{1};
    end

    for j=1:numel(phases{i,2})
    results_single_stance{i,j}=phase_names{2};
    end

    for j=1:numel(phases{i,3})
    results_terminal_stance{i,j}=phase_names{3};
    end

    for j=1:numel(phases{i,4})
    results_swing{i,j}=phase_names{4};
    end
end

values_phases_1=horzcat(results_heel_contact(1,:),results_single_stance(1,:),results_terminal_stance(1,:),results_swing(1,:));
values_phases_transpose_1=transpose(values_phases_1);
non_empty_cells_1 = ~cellfun('isempty', values_phases_transpose_1);%find the non-empty cells
results_phases_1 = values_phases_transpose_1(non_empty_cells_1);

values_phases_2=horzcat(results_heel_contact(2,:),results_single_stance(2,:),results_terminal_stance(2,:),results_swing(2,:));
values_phases_transpose_2=transpose(values_phases_2);
non_empty_cells_2 = ~cellfun('isempty', values_phases_transpose_2);%find the non-empty cells
results_phases_2 = values_phases_transpose_2(non_empty_cells_2);

values_phases_3=horzcat(results_heel_contact(3,:),results_single_stance(3,:),results_terminal_stance(3,:),results_swing(3,:));
values_phases_transpose_3=transpose(values_phases_3);
non_empty_cells_3 = ~cellfun('isempty', values_phases_transpose_3);%find the non-empty cells
results_phases_3 = values_phases_transpose_3(non_empty_cells_3);


%% R_RF
R_RF_new_1=[R_RF_tot(1,:)',results_phases_1(:,1)];

R_RF_1={};

for riga = 1:size(R_RF_new_1, 1)
    if strcmp(strtrim(R_RF_new_1{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_1{riga,2}), 'HEEL CONTACT')
        R_RF_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_1{riga,2}), 'HEEL CONTACT')
        R_RF_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_1{riga,2}), 'SINGLE STANCE')
        R_RF_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_1{riga,2}), 'SINGLE STANCE')
        R_RF_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_1{riga,2}), 'TERMINAL STANCE')
        R_RF_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_1{riga,2}), 'TERMINAL STANCE')
        R_RF_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_1{riga,2}), 'SWING')
        R_RF_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_RF_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_1{riga,2}), 'SWING')
        R_RF_1{end+1} = 'NO SW' ;
        end
end

num_HC_RF_1 = sum(strcmp(R_RF_1, 'HC'));
num_NO_HC_RF_1 = sum(strcmp(R_RF_1, 'NO HC'));
percent_HC_RF_1 = [num2str(round(num_HC_RF_1 / (num_HC_RF_1 + num_NO_HC_RF_1) * 100)),'%'];

num_ST_RF_1 = sum(strcmp(R_RF_1, 'ST'));
num_NO_ST_RF_1 = sum(strcmp(R_RF_1, 'NO ST'));
percent_ST_RF_1 = [num2str(round(num_ST_RF_1 / (num_ST_RF_1 + num_NO_ST_RF_1) * 100)),'%'];

num_TER_RF_1 = sum(strcmp(R_RF_1, 'TER'));
num_NO_TER_RF_1 = sum(strcmp(R_RF_1, 'NO TER'));
percent_TER_RF_1 = [num2str(round(num_TER_RF_1/ (num_TER_RF_1 + num_NO_TER_RF_1) * 100)),'%'];

num_SW_RF_1 = sum(strcmp(R_RF_1, 'SW'));
num_NO_SW_RF_1 = sum(strcmp(R_RF_1, 'NO SW'));
percent_SW_RF_1= [num2str(round(num_SW_RF_1 / (num_SW_RF_1 + num_NO_SW_RF_1) * 100)),'%'];
R_RF_PERC_1={percent_HC_RF_1;percent_ST_RF_1;percent_TER_RF_1;percent_SW_RF_1};


R_RF_new_2=[R_RF_tot(2,:)',results_phases_2(:,1)];
R_RF_2={};

for riga = 1:size(R_RF_new_2, 1)
    if strcmp(strtrim(R_RF_new_2{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_2{riga,2}), 'HEEL CONTACT')
        R_RF_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_2{riga,2}), 'HEEL CONTACT')
        R_RF_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_2{riga,2}), 'SINGLE STANCE')
        R_RF_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_2{riga,2}), 'SINGLE STANCE')
        R_RF_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_2{riga,2}), 'TERMINAL STANCE')
        R_RF_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_2{riga,2}), 'TERMINAL STANCE')
        R_RF_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_2{riga,2}), 'SWING')
        R_RF_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_RF_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_2{riga,2}), 'SWING')
        R_RF_2{end+1} = 'NO SW' ;
        end
end

num_HC_RF_2 = sum(strcmp(R_RF_2, 'HC'));
num_NO_HC_RF_2 = sum(strcmp(R_RF_2, 'NO HC'));
percent_HC_RF_2 = [num2str(round(num_HC_RF_2 / (num_HC_RF_2 + num_NO_HC_RF_2) * 100)),'%'];

num_ST_RF_2= sum(strcmp(R_RF_2, 'ST'));
num_NO_ST_RF_2 = sum(strcmp(R_RF_2, 'NO ST'));
percent_ST_RF_2 =[num2str(round(num_ST_RF_2 / (num_ST_RF_2 + num_NO_ST_RF_2) * 100)),'%'];

num_TER_RF_2 = sum(strcmp(R_RF_2, 'TER'));
num_NO_TER_RF_2 = sum(strcmp(R_RF_2, 'NO TER'));
percent_TER_RF_2 = [num2str(round(num_TER_RF_2/ (num_TER_RF_2 + num_NO_TER_RF_2) * 100)),'%'];

num_SW_RF_2 = sum(strcmp(R_RF_2, 'SW'));
num_NO_SW_RF_2 = sum(strcmp(R_RF_2, 'NO SW'));
percent_SW_RF_2= [num2str(round(num_SW_RF_2 / (num_SW_RF_2 + num_NO_SW_RF_2) * 100)),'%'];
R_RF_PERC_2={percent_HC_RF_2;percent_ST_RF_2;percent_TER_RF_2;percent_SW_RF_2};

R_RF_new_3=[R_RF_tot(3,:)',results_phases_3(:,1)];
R_RF_3={};

for riga = 1:size(R_RF_new_3, 1)
    if strcmp(strtrim(R_RF_new_3{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_3{riga,2}), 'HEEL CONTACT')
        R_RF_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_3{riga,2}), 'HEEL CONTACT')
        R_RF_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_3{riga,2}), 'SINGLE STANCE')
        R_RF_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_3{riga,2}), 'SINGLE STANCE')
        R_RF_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_3{riga,2}), 'TERMINAL STANCE')
        R_RF_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_3{riga,2}), 'TERMINAL STANCE')
        R_RF_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'ON') && strcmp(strtrim(R_RF_new_3{riga,2}), 'SWING')
        R_RF_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_RF_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_RF_new_3{riga,2}), 'SWING')
        R_RF_3{end+1} = 'NO SW' ;
        end
end

num_HC_RF_3 = sum(strcmp(R_RF_3, 'HC'));
num_NO_HC_RF_3 = sum(strcmp(R_RF_3, 'NO HC'));
percent_HC_RF_3 = [num2str(round(num_HC_RF_3 / (num_HC_RF_3 + num_NO_HC_RF_3) * 100)),'%'];


num_ST_RF_3= sum(strcmp(R_RF_3, 'ST'));
num_NO_ST_RF_3 = sum(strcmp(R_RF_3, 'NO ST'));
percent_ST_RF_3 =[num2str(round(num_ST_RF_3 / (num_ST_RF_3 + num_NO_ST_RF_3) * 100)),'%'];

num_TER_RF_3 = sum(strcmp(R_RF_3, 'TER'));
num_NO_TER_RF_3 = sum(strcmp(R_RF_3, 'NO TER'));
percent_TER_RF_3 = [num2str(round(num_TER_RF_3/ (num_TER_RF_3 + num_NO_TER_RF_3) * 100)),'%'];


num_SW_RF_3 = sum(strcmp(R_RF_3, 'SW'));
num_NO_SW_RF_3 = sum(strcmp(R_RF_3, 'NO SW'));
percent_SW_RF_3= [num2str(round(num_SW_RF_3 / (num_SW_RF_3 + num_NO_SW_RF_3) * 100)),'%'];
R_RF_PERC_3={percent_HC_RF_3;percent_ST_RF_3;percent_TER_RF_3;percent_SW_RF_3};

R_RF_PERC=[R_RF_PERC_1;R_RF_PERC_2;R_RF_PERC_3];

%% R_VL
R_VL_new_1=[R_VL_tot(1,:)',results_phases_1(:,1)];

R_VL_1={};

for riga = 1:size(R_VL_new_1, 1)
    if strcmp(strtrim(R_VL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_1{riga,2}), 'HEEL CONTACT')
        R_VL_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_1{riga,2}), 'HEEL CONTACT')
        R_VL_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_1{riga,2}), 'SINGLE STANCE')
        R_VL_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_1{riga,2}), 'SINGLE STANCE')
        R_VL_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_1{riga,2}), 'TERMINAL STANCE')
        R_VL_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_1{riga,2}), 'TERMINAL STANCE')
        R_VL_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_1{riga,2}), 'SWING')
        R_VL_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_VL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_1{riga,2}), 'SWING')
        R_VL_1{end+1} = 'NO SW' ;
        end
end

num_HC_VL_1 = sum(strcmp(R_VL_1, 'HC'));
num_NO_HC_VL_1 = sum(strcmp(R_VL_1, 'NO HC'));
percent_HC_VL_1 = [num2str(round(num_HC_VL_1 / (num_HC_VL_1 + num_NO_HC_VL_1) * 100)),'%'];

num_ST_VL_1 = sum(strcmp(R_VL_1, 'ST'));
num_NO_VL_1 = sum(strcmp(R_VL_1, 'NO ST'));
percent_ST_VL_1 =[num2str(round(num_ST_VL_1 / (num_ST_VL_1 + num_NO_VL_1) * 100)),'%'];

num_TER_VL_1 = sum(strcmp(R_VL_1, 'TER'));
num_NO_TER_VL_1 = sum(strcmp(R_VL_1, 'NO TER'));
percent_TER_VL_1 = [num2str(round(num_TER_VL_1/ (num_TER_VL_1 + num_NO_TER_VL_1) * 100)),'%'];

num_SW_VL_1 = sum(strcmp(R_VL_1, 'SW'));
num_NO_SW_VL_1 = sum(strcmp(R_VL_1, 'NO SW'));
percent_SW_VL_1= [num2str(round(num_SW_VL_1 / (num_SW_VL_1 + num_NO_SW_VL_1) * 100)),'%'];
R_VL_PERC_1={percent_HC_VL_1;percent_ST_VL_1;percent_TER_VL_1;percent_SW_VL_1};


R_VL_new_2=[R_VL_tot(2,:)',results_phases_2(:,1)];
R_VL_2={};

for riga = 1:size(R_VL_new_2, 1)
    if strcmp(strtrim(R_VL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_2{riga,2}), 'HEEL CONTACT')
        R_VL_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_2{riga,2}), 'HEEL CONTACT')
        R_VL_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_2{riga,2}), 'SINGLE STANCE')
        R_VL_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_2{riga,2}), 'SINGLE STANCE')
        R_VL_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_2{riga,2}), 'TERMINAL STANCE')
        R_VL_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_2{riga,2}), 'TERMINAL STANCE')
        R_VL_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_2{riga,2}), 'SWING')
        R_VL_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_VL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_2{riga,2}), 'SWING')
        R_VL_2{end+1} = 'NO SW' ;
        end
end

num_HC_VL_2 = sum(strcmp(R_VL_2, 'HC'));
num_NO_HC_VL_2 = sum(strcmp(R_VL_2, 'NO HC'));
percent_HC_VL_2 = [num2str(round(num_HC_VL_2 / (num_HC_VL_2 + num_NO_HC_VL_2) * 100)),'%'];

num_ST_VL_2= sum(strcmp(R_VL_2, 'ST'));
num_NO_ST_VL_2 = sum(strcmp(R_VL_2, 'NO ST'));
percent_ST_VL_2 =[num2str(round(num_ST_VL_2 / (num_ST_VL_2 + num_NO_ST_VL_2) * 100)),'%'];

num_TER_VL_2 = sum(strcmp(R_VL_2, 'TER'));
num_NO_TER_VL_2 = sum(strcmp(R_VL_2, 'NO TER'));
percent_TER_VL_2 = [num2str(round(num_TER_VL_2/ (num_TER_VL_2 + num_NO_TER_VL_2) * 100)),'%'];

num_SW_VL_2 = sum(strcmp(R_VL_2, 'SW'));
num_NO_SW_VL_2 = sum(strcmp(R_VL_2, 'NO SW'));
percent_SW_VL_2= [num2str(round(num_SW_VL_2 / (num_SW_VL_2 + num_NO_SW_VL_2) * 100)),'%'];
R_VL_PERC_2={percent_HC_VL_2;percent_ST_VL_2;percent_TER_VL_2;percent_SW_VL_2};

R_VL_new_3=[R_VL_tot(3,:)',results_phases_3(:,1)];
R_VL_3={};

for riga = 1:size(R_VL_new_3, 1)
    if strcmp(strtrim(R_VL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_3{riga,2}), 'HEEL CONTACT')
        R_VL_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_3{riga,2}), 'HEEL CONTACT')
        R_VL_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_3{riga,2}), 'SINGLE STANCE')
        R_VL_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_3{riga,2}), 'SINGLE STANCE')
        R_VL_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_3{riga,2}), 'TERMINAL STANCE')
        R_VL_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_3{riga,2}), 'TERMINAL STANCE')
        R_VL_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VL_new_3{riga,2}), 'SWING')
        R_VL_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_VL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VL_new_3{riga,2}), 'SWING')
        R_VL_3{end+1} = 'NO SW' ;
        end
end

num_HC_VL_3 = sum(strcmp(R_VL_3, 'HC'));
num_NO_HC_VL_3 = sum(strcmp(R_VL_3, 'NO HC'));
percent_HC_VL_3 = [num2str(round(num_HC_VL_3 / (num_HC_VL_3 + num_NO_HC_VL_3) * 100)),'%'];


num_ST_VL_3= sum(strcmp(R_VL_3, 'ST'));
num_NO_ST_VL_3 = sum(strcmp(R_VL_3, 'NO ST'));
percent_ST_VL_3 =[num2str(round(num_ST_VL_3 / (num_ST_VL_3 + num_NO_ST_VL_3) * 100)),'%'];

num_TER_VL_3 = sum(strcmp(R_VL_3, 'TER'));
num_NO_TER_VL_3 = sum(strcmp(R_VL_3, 'NO TER'));
percent_TER_VL_3 = [num2str(round(num_TER_VL_3/ (num_TER_VL_3 + num_NO_TER_VL_3) * 100)),'%'];


num_SW_VL_3 = sum(strcmp(R_VL_3, 'SW'));
num_NO_SW_VL_3 = sum(strcmp(R_VL_3, 'NO SW'));
percent_SW_VL_3= [num2str(round(num_SW_VL_3 / (num_SW_VL_3 + num_NO_SW_VL_3) * 100)),'%'];
R_VL_PERC_3={percent_HC_VL_3;percent_ST_VL_3;percent_TER_VL_3;percent_SW_VL_3};

R_VL_PERC=[R_VL_PERC_1;R_VL_PERC_2;R_VL_PERC_3];

%% R_VM
R_VM_new_1=[R_VM_tot(1,:)',results_phases_1(:,1)];

R_VM_1={};

for riga = 1:size(R_VM_new_1, 1)
    if strcmp(strtrim(R_VM_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_1{riga,2}), 'HEEL CONTACT')
        R_VM_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_1{riga,2}), 'HEEL CONTACT')
        R_VM_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_1{riga,2}), 'SINGLE STANCE')
        R_VM_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_1{riga,2}), 'SINGLE STANCE')
        R_VM_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_1{riga,2}), 'TERMINAL STANCE')
        R_VM_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_1{riga,2}), 'TERMINAL STANCE')
        R_VM_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_1{riga,2}), 'SWING')
        R_VM_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_VM_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_1{riga,2}), 'SWING')
        R_VM_1{end+1} = 'NO SW' ;
        end
end

num_HC_VM_1 = sum(strcmp(R_VM_1, 'HC'));
num_NO_HC_VM_1 = sum(strcmp(R_VM_1, 'NO HC'));
percent_HC_VM_1 = [num2str(round(num_HC_VM_1 / (num_HC_VM_1 + num_NO_HC_VM_1) * 100)),'%'];

num_ST_VM_1 = sum(strcmp(R_VM_1, 'ST'));
num_NO_VM_1 = sum(strcmp(R_VM_1, 'NO ST'));
percent_ST_VM_1 =[num2str(round(num_ST_VM_1 / (num_ST_VM_1 + num_NO_VM_1) * 100)),'%'];

num_TER_VM_1 = sum(strcmp(R_VM_1, 'TER'));
num_NO_TER_VM_1 = sum(strcmp(R_VM_1, 'NO TER'));
percent_TER_VM_1 = [num2str(round(num_TER_VM_1/ (num_TER_VM_1 + num_NO_TER_VM_1) * 100)),'%'];

num_SW_VM_1 = sum(strcmp(R_VM_1, 'SW'));
num_NO_SW_VM_1 = sum(strcmp(R_VM_1, 'NO SW'));
percent_SW_VM_1= [num2str(round(num_SW_VM_1 / (num_SW_VM_1 + num_NO_SW_VM_1) * 100)),'%'];
R_VM_PERC_1={percent_HC_VM_1;percent_ST_VM_1;percent_TER_VM_1;percent_SW_VM_1};


R_VM_new_2=[R_VM_tot(2,:)',results_phases_2(:,1)];
R_VM_2={};

for riga = 1:size(R_VM_new_2, 1)
    if strcmp(strtrim(R_VM_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_2{riga,2}), 'HEEL CONTACT')
        R_VM_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_2{riga,2}), 'HEEL CONTACT')
        R_VM_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_2{riga,2}), 'SINGLE STANCE')
        R_VM_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_2{riga,2}), 'SINGLE STANCE')
        R_VM_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_2{riga,2}), 'TERMINAL STANCE')
        R_VM_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_2{riga,2}), 'TERMINAL STANCE')
        R_VM_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_2{riga,2}), 'SWING')
        R_VM_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_VM_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_2{riga,2}), 'SWING')
        R_VM_2{end+1} = 'NO SW' ;
        end
end

num_HC_VM_2 = sum(strcmp(R_VM_2, 'HC'));
num_NO_HC_VM_2 = sum(strcmp(R_VM_2, 'NO HC'));
percent_HC_VM_2 = [num2str(round(num_HC_VM_2 / (num_HC_VM_2 + num_NO_HC_VM_2) * 100)),'%'];

num_ST_VM_2= sum(strcmp(R_VM_2, 'ST'));
num_NO_ST_VM_2 = sum(strcmp(R_VM_2, 'NO ST'));
percent_ST_VM_2 =[num2str(round(num_ST_VM_2 / (num_ST_VM_2 + num_NO_ST_VM_2) * 100)),'%'];

num_TER_VM_2 = sum(strcmp(R_VM_2, 'TER'));
num_NO_TER_VM_2 = sum(strcmp(R_VM_2, 'NO TER'));
percent_TER_VM_2 = [num2str(round(num_TER_VM_2/ (num_TER_VM_2 + num_NO_TER_VM_2) * 100)),'%'];

num_SW_VM_2 = sum(strcmp(R_VM_2, 'SW'));
num_NO_SW_VM_2 = sum(strcmp(R_VM_2, 'NO SW'));
percent_SW_VM_2= [num2str(round(num_SW_VM_2 / (num_SW_VM_2 + num_NO_SW_VM_2) * 100)),'%'];
R_VM_PERC_2={percent_HC_VM_2;percent_ST_VM_2;percent_TER_VM_2;percent_SW_VM_2};

R_VM_new_3=[R_VM_tot(3,:)',results_phases_3(:,1)];
R_VM_3={};

for riga = 1:size(R_VM_new_3, 1)
    if strcmp(strtrim(R_VM_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_3{riga,2}), 'HEEL CONTACT')
        R_VM_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_3{riga,2}), 'HEEL CONTACT')
        R_VM_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_3{riga,2}), 'SINGLE STANCE')
        R_VM_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_3{riga,2}), 'SINGLE STANCE')
        R_VM_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_3{riga,2}), 'TERMINAL STANCE')
        R_VM_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_3{riga,2}), 'TERMINAL STANCE')
        R_VM_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'ON') && strcmp(strtrim(R_VM_new_3{riga,2}), 'SWING')
        R_VM_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_VM_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_VM_new_3{riga,2}), 'SWING')
        R_VM_3{end+1} = 'NO SW' ;
        end
end

num_HC_VM_3 = sum(strcmp(R_VM_3, 'HC'));
num_NO_HC_VM_3 = sum(strcmp(R_VM_3, 'NO HC'));
percent_HC_VM_3 = [num2str(round(num_HC_VM_3 / (num_HC_VM_3 + num_NO_HC_VM_3) * 100)),'%'];


num_ST_VM_3= sum(strcmp(R_VM_3, 'ST'));
num_NO_ST_VM_3 = sum(strcmp(R_VM_3, 'NO ST'));
percent_ST_VM_3 =[num2str(round(num_ST_VM_3 / (num_ST_VM_3 + num_NO_ST_VM_3) * 100)),'%'];

num_TER_VM_3 = sum(strcmp(R_VM_3, 'TER'));
num_NO_TER_VM_3 = sum(strcmp(R_VM_3, 'NO TER'));
percent_TER_VM_3 = [num2str(round(num_TER_VM_3/ (num_TER_VM_3 + num_NO_TER_VM_3) * 100)),'%'];


num_SW_VM_3 = sum(strcmp(R_VM_3, 'SW'));
num_NO_SW_VM_3 = sum(strcmp(R_VM_3, 'NO SW'));
percent_SW_VM_3= [num2str(round(num_SW_VM_3 / (num_SW_VM_3 + num_NO_SW_VM_3) * 100)),'%'];
R_VM_PERC_3={percent_HC_VM_3;percent_ST_VM_3;percent_TER_VM_3;percent_SW_VM_3};

R_VM_PERC=[R_VM_PERC_1;R_VM_PERC_2;R_VM_PERC_3];


%% R_G
R_G_new_1=[R_G_tot(1,:)',results_phases_1(:,1)];

R_G_1={};

for riga = 1:size(R_G_new_1, 1)
    if strcmp(strtrim(R_G_new_1{riga,1}), 'ON') && strcmp(strtrim(R_G_new_1{riga,2}), 'HEEL CONTACT')
        R_G_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_1{riga,2}), 'HEEL CONTACT')
        R_G_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'ON') && strcmp(strtrim(R_G_new_1{riga,2}), 'SINGLE STANCE')
        R_G_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_1{riga,2}), 'SINGLE STANCE')
        R_G_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'ON') && strcmp(strtrim(R_G_new_1{riga,2}), 'TERMINAL STANCE')
        R_G_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_1{riga,2}), 'TERMINAL STANCE')
        R_G_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_G_new_1{riga,1}), 'ON') && strcmp(strtrim(R_G_new_1{riga,2}), 'SWING')
        R_G_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_G_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_1{riga,2}), 'SWING')
        R_G_1{end+1} = 'NO SW' ;
        end
end

num_HC_G_1 = sum(strcmp(R_G_1, 'HC'));
num_NO_HC_G_1 = sum(strcmp(R_G_1, 'NO HC'));
percent_HC_G_1 = [num2str(round(num_HC_G_1 / (num_HC_G_1 + num_NO_HC_G_1) * 100)),'%'];

num_ST_G_1 = sum(strcmp(R_G_1, 'ST'));
num_NO_G_1 = sum(strcmp(R_G_1, 'NO ST'));
percent_ST_G_1 =[num2str(round(num_ST_G_1 / (num_ST_G_1 + num_NO_G_1) * 100)),'%'];

num_TER_G_1 = sum(strcmp(R_G_1, 'TER'));
num_NO_TER_G_1 = sum(strcmp(R_G_1, 'NO TER'));
percent_TER_G_1 = [num2str(round(num_TER_G_1/ (num_TER_G_1 + num_NO_TER_G_1) * 100)),'%'];

num_SW_G_1 = sum(strcmp(R_G_1, 'SW'));
num_NO_SW_G_1 = sum(strcmp(R_G_1, 'NO SW'));
percent_SW_G_1= [num2str(round(num_SW_G_1 / (num_SW_G_1 + num_NO_SW_G_1) * 100)),'%'];
R_G_PERC_1={percent_HC_G_1;percent_ST_G_1;percent_TER_G_1;percent_SW_G_1};


R_G_new_2=[R_G_tot(2,:)',results_phases_2(:,1)];
R_G_2={};

for riga = 1:size(R_G_new_2, 1)
    if strcmp(strtrim(R_G_new_2{riga,1}), 'ON') && strcmp(strtrim(R_G_new_2{riga,2}), 'HEEL CONTACT')
        R_G_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_2{riga,2}), 'HEEL CONTACT')
        R_G_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'ON') && strcmp(strtrim(R_G_new_2{riga,2}), 'SINGLE STANCE')
        R_G_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_2{riga,2}), 'SINGLE STANCE')
        R_G_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'ON') && strcmp(strtrim(R_G_new_2{riga,2}), 'TERMINAL STANCE')
        R_G_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_2{riga,2}), 'TERMINAL STANCE')
        R_G_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_G_new_2{riga,1}), 'ON') && strcmp(strtrim(R_G_new_2{riga,2}), 'SWING')
        R_G_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_G_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_2{riga,2}), 'SWING')
        R_G_2{end+1} = 'NO SW' ;
        end
end

num_HC_G_2 = sum(strcmp(R_G_2, 'HC'));
num_NO_HC_G_2 = sum(strcmp(R_G_2, 'NO HC'));
percent_HC_G_2 = [num2str(round(num_HC_G_2 / (num_HC_G_2 + num_NO_HC_G_2) * 100)),'%'];

num_ST_G_2= sum(strcmp(R_G_2, 'ST'));
num_NO_ST_G_2 = sum(strcmp(R_G_2, 'NO ST'));
percent_ST_G_2 =[num2str(round(num_ST_G_2 / (num_ST_G_2 + num_NO_ST_G_2) * 100)),'%'];

num_TER_G_2 = sum(strcmp(R_G_2, 'TER'));
num_NO_TER_G_2 = sum(strcmp(R_G_2, 'NO TER'));
percent_TER_G_2 = [num2str(round(num_TER_G_2/ (num_TER_G_2 + num_NO_TER_G_2) * 100)),'%'];

num_SW_G_2 = sum(strcmp(R_G_2, 'SW'));
num_NO_SW_G_2 = sum(strcmp(R_G_2, 'NO SW'));
percent_SW_G_2= [num2str(round(num_SW_G_2 / (num_SW_G_2 + num_NO_SW_G_2) * 100)),'%'];
R_G_PERC_2={percent_HC_G_2;percent_ST_G_2;percent_TER_G_2;percent_SW_G_2};

R_G_new_3=[R_G_tot(3,:)',results_phases_3(:,1)];
R_G_3={};

for riga = 1:size(R_G_new_3, 1)
    if strcmp(strtrim(R_G_new_3{riga,1}), 'ON') && strcmp(strtrim(R_G_new_3{riga,2}), 'HEEL CONTACT')
        R_G_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_3{riga,2}), 'HEEL CONTACT')
        R_G_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'ON') && strcmp(strtrim(R_G_new_3{riga,2}), 'SINGLE STANCE')
        R_G_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_3{riga,2}), 'SINGLE STANCE')
        R_G_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'ON') && strcmp(strtrim(R_G_new_3{riga,2}), 'TERMINAL STANCE')
        R_G_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_3{riga,2}), 'TERMINAL STANCE')
        R_G_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_G_new_3{riga,1}), 'ON') && strcmp(strtrim(R_G_new_3{riga,2}), 'SWING')
        R_G_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_G_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_G_new_3{riga,2}), 'SWING')
        R_G_3{end+1} = 'NO SW' ;
        end
end

num_HC_G_3 = sum(strcmp(R_G_3, 'HC'));
num_NO_HC_G_3 = sum(strcmp(R_G_3, 'NO HC'));
percent_HC_G_3 = [num2str(round(num_HC_G_3 / (num_HC_G_3 + num_NO_HC_G_3) * 100)),'%'];


num_ST_G_3= sum(strcmp(R_G_3, 'ST'));
num_NO_ST_G_3 = sum(strcmp(R_G_3, 'NO ST'));
percent_ST_G_3 =[num2str(round(num_ST_G_3 / (num_ST_G_3 + num_NO_ST_G_3) * 100)),'%'];

num_TER_G_3 = sum(strcmp(R_G_3, 'TER'));
num_NO_TER_G_3 = sum(strcmp(R_G_3, 'NO TER'));
percent_TER_G_3 = [num2str(round(num_TER_G_3/ (num_TER_G_3 + num_NO_TER_G_3) * 100)),'%'];


num_SW_G_3 = sum(strcmp(R_G_3, 'SW'));
num_NO_SW_G_3 = sum(strcmp(R_G_3, 'NO SW'));
percent_SW_G_3= [num2str(round(num_SW_G_3 / (num_SW_G_3 + num_NO_SW_G_3) * 100)),'%'];
R_G_PERC_3={percent_HC_G_3;percent_ST_G_3;percent_TER_G_3;percent_SW_G_3};

R_G_PERC=[R_G_PERC_1;R_G_PERC_2;R_G_PERC_3];

%% R_BFCL
R_BFCL_new_1=[R_BFCL_tot(1,:)',results_phases_1(:,1)];

R_BFCL_1={};

for riga = 1:size(R_BFCL_new_1, 1)
    if strcmp(strtrim(R_BFCL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'HEEL CONTACT')
        R_BFCL_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'HEEL CONTACT')
        R_BFCL_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'SINGLE STANCE')
        R_BFCL_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'SINGLE STANCE')
        R_BFCL_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'TERMINAL STANCE')
        R_BFCL_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'TERMINAL STANCE')
        R_BFCL_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'SWING')
        R_BFCL_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_BFCL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_1{riga,2}), 'SWING')
        R_BFCL_1{end+1} = 'NO SW' ;
        end
end

num_HC_BFCL_1 = sum(strcmp(R_BFCL_1, 'HC'));
num_NO_HC_BFCL_1 = sum(strcmp(R_BFCL_1, 'NO HC'));
percent_HC_BFCL_1 = [num2str(round(num_HC_BFCL_1 / (num_HC_BFCL_1 + num_NO_HC_BFCL_1) * 100)),'%'];

num_ST_BFCL_1 = sum(strcmp(R_BFCL_1, 'ST'));
num_NO_BFCL_1 = sum(strcmp(R_BFCL_1, 'NO ST'));
percent_ST_BFCL_1 =[num2str(round(num_ST_BFCL_1 / (num_ST_BFCL_1 + num_NO_BFCL_1) * 100)),'%'];

num_TER_BFCL_1 = sum(strcmp(R_BFCL_1, 'TER'));
num_NO_TER_BFCL_1 = sum(strcmp(R_BFCL_1, 'NO TER'));
percent_TER_BFCL_1 = [num2str(round(num_TER_BFCL_1/ (num_TER_BFCL_1 + num_NO_TER_BFCL_1) * 100)),'%'];

num_SW_BFCL_1 = sum(strcmp(R_BFCL_1, 'SW'));
num_NO_SW_BFCL_1 = sum(strcmp(R_BFCL_1, 'NO SW'));
percent_SW_BFCL_1= [num2str(round(num_SW_BFCL_1 / (num_SW_BFCL_1 + num_NO_SW_BFCL_1) * 100)),'%'];
R_BFCL_PERC_1={percent_HC_BFCL_1;percent_ST_BFCL_1;percent_TER_BFCL_1;percent_SW_BFCL_1};


R_BFCL_new_2=[R_BFCL_tot(2,:)',results_phases_2(:,1)];
R_BFCL_2={};

for riga = 1:size(R_BFCL_new_2, 1)
    if strcmp(strtrim(R_BFCL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'HEEL CONTACT')
        R_BFCL_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'HEEL CONTACT')
        R_BFCL_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'SINGLE STANCE')
        R_BFCL_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'SINGLE STANCE')
        R_BFCL_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'TERMINAL STANCE')
        R_BFCL_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'TERMINAL STANCE')
        R_BFCL_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'SWING')
        R_BFCL_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_BFCL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_2{riga,2}), 'SWING')
        R_BFCL_2{end+1} = 'NO SW' ;
        end
end

num_HC_BFCL_2 = sum(strcmp(R_BFCL_2, 'HC'));
num_NO_HC_BFCL_2 = sum(strcmp(R_BFCL_2, 'NO HC'));
percent_HC_BFCL_2 = [num2str(round(num_HC_BFCL_2 / (num_HC_BFCL_2 + num_NO_HC_BFCL_2) * 100)),'%'];

num_ST_BFCL_2= sum(strcmp(R_BFCL_2, 'ST'));
num_NO_ST_BFCL_2 = sum(strcmp(R_BFCL_2, 'NO ST'));
percent_ST_BFCL_2 =[num2str(round(num_ST_BFCL_2 / (num_ST_BFCL_2 + num_NO_ST_BFCL_2) * 100)),'%'];

num_TER_BFCL_2 = sum(strcmp(R_BFCL_2, 'TER'));
num_NO_TER_BFCL_2 = sum(strcmp(R_BFCL_2, 'NO TER'));
percent_TER_BFCL_2 = [num2str(round(num_TER_BFCL_2/ (num_TER_BFCL_2 + num_NO_TER_BFCL_2) * 100)),'%'];

num_SW_BFCL_2 = sum(strcmp(R_BFCL_2, 'SW'));
num_NO_SW_BFCL_2 = sum(strcmp(R_BFCL_2, 'NO SW'));
percent_SW_BFCL_2= [num2str(round(num_SW_BFCL_2 / (num_SW_BFCL_2 + num_NO_SW_BFCL_2) * 100)),'%'];
R_BFCL_PERC_2={percent_HC_BFCL_2;percent_ST_BFCL_2;percent_TER_BFCL_2;percent_SW_BFCL_2};

R_BFCL_new_3=[R_BFCL_tot(3,:)',results_phases_3(:,1)];
R_BFCL_3={};

for riga = 1:size(R_BFCL_new_3, 1)
    if strcmp(strtrim(R_BFCL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'HEEL CONTACT')
        R_BFCL_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'HEEL CONTACT')
        R_BFCL_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'SINGLE STANCE')
        R_BFCL_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'SINGLE STANCE')
        R_BFCL_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'TERMINAL STANCE')
        R_BFCL_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'TERMINAL STANCE')
        R_BFCL_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'SWING')
        R_BFCL_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_BFCL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_BFCL_new_3{riga,2}), 'SWING')
        R_BFCL_3{end+1} = 'NO SW' ;
        end
end

num_HC_BFCL_3 = sum(strcmp(R_BFCL_3, 'HC'));
num_NO_HC_BFCL_3 = sum(strcmp(R_BFCL_3, 'NO HC'));
percent_HC_BFCL_3 = [num2str(round(num_HC_BFCL_3 / (num_HC_BFCL_3 + num_NO_HC_BFCL_3) * 100)),'%'];


num_ST_BFCL_3= sum(strcmp(R_BFCL_3, 'ST'));
num_NO_ST_BFCL_3 = sum(strcmp(R_BFCL_3, 'NO ST'));
percent_ST_BFCL_3 =[num2str(round(num_ST_BFCL_3 / (num_ST_BFCL_3 + num_NO_ST_BFCL_3) * 100)),'%'];

num_TER_BFCL_3 = sum(strcmp(R_BFCL_3, 'TER'));
num_NO_TER_BFCL_3 = sum(strcmp(R_BFCL_3, 'NO TER'));
percent_TER_BFCL_3 = [num2str(round(num_TER_BFCL_3/ (num_TER_BFCL_3 + num_NO_TER_BFCL_3) * 100)),'%'];


num_SW_BFCL_3 = sum(strcmp(R_BFCL_3, 'SW'));
num_NO_SW_BFCL_3 = sum(strcmp(R_BFCL_3, 'NO SW'));
percent_SW_BFCL_3= [num2str(round(num_SW_BFCL_3 / (num_SW_BFCL_3 + num_NO_SW_BFCL_3) * 100)),'%'];
R_BFCL_PERC_3={percent_HC_BFCL_3;percent_ST_BFCL_3;percent_TER_BFCL_3;percent_SW_BFCL_3};

R_BFCL_PERC=[R_BFCL_PERC_1;R_BFCL_PERC_2;R_BFCL_PERC_3];

%% R_TA
R_TA_new_1=[R_TA_tot(1,:)',results_phases_1(:,1)];

R_TA_1={};

for riga = 1:size(R_TA_new_1, 1)
    if strcmp(strtrim(R_TA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_1{riga,2}), 'HEEL CONTACT')
        R_TA_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_1{riga,2}), 'HEEL CONTACT')
        R_TA_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_1{riga,2}), 'SINGLE STANCE')
        R_TA_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_1{riga,2}), 'SINGLE STANCE')
        R_TA_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_1{riga,2}), 'TERMINAL STANCE')
        R_TA_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_1{riga,2}), 'TERMINAL STANCE')
        R_TA_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_1{riga,2}), 'SWING')
        R_TA_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_TA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_1{riga,2}), 'SWING')
        R_TA_1{end+1} = 'NO SW' ;
        end
end

num_HC_TA_1 = sum(strcmp(R_TA_1, 'HC'));
num_NO_HC_TA_1 = sum(strcmp(R_TA_1, 'NO HC'));
percent_HC_TA_1 = [num2str(round(num_HC_TA_1 / (num_HC_TA_1 + num_NO_HC_TA_1) * 100)),'%'];

num_ST_TA_1 = sum(strcmp(R_TA_1, 'ST'));
num_NO_TA_1 = sum(strcmp(R_TA_1, 'NO ST'));
percent_ST_TA_1 =[num2str(round(num_ST_TA_1 / (num_ST_TA_1 + num_NO_TA_1) * 100)),'%'];

num_TER_TA_1 = sum(strcmp(R_TA_1, 'TER'));
num_NO_TER_TA_1 = sum(strcmp(R_TA_1, 'NO TER'));
percent_TER_TA_1 = [num2str(round(num_TER_TA_1/ (num_TER_TA_1 + num_NO_TER_TA_1) * 100)),'%'];

num_SW_TA_1 = sum(strcmp(R_TA_1, 'SW'));
num_NO_SW_TA_1 = sum(strcmp(R_TA_1, 'NO SW'));
percent_SW_TA_1= [num2str(round(num_SW_TA_1 / (num_SW_TA_1 + num_NO_SW_TA_1) * 100)),'%'];
R_TA_PERC_1={percent_HC_TA_1;percent_ST_TA_1;percent_TER_TA_1;percent_SW_TA_1};


R_TA_new_2=[R_TA_tot(2,:)',results_phases_2(:,1)];
R_TA_2={};

for riga = 1:size(R_TA_new_2, 1)
    if strcmp(strtrim(R_TA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_2{riga,2}), 'HEEL CONTACT')
        R_TA_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_2{riga,2}), 'HEEL CONTACT')
        R_TA_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_2{riga,2}), 'SINGLE STANCE')
        R_TA_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_2{riga,2}), 'SINGLE STANCE')
        R_TA_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_2{riga,2}), 'TERMINAL STANCE')
        R_TA_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_2{riga,2}), 'TERMINAL STANCE')
        R_TA_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_2{riga,2}), 'SWING')
        R_TA_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_TA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_2{riga,2}), 'SWING')
        R_TA_2{end+1} = 'NO SW' ;
        end
end

num_HC_TA_2 = sum(strcmp(R_TA_2, 'HC'));
num_NO_HC_TA_2 = sum(strcmp(R_TA_2, 'NO HC'));
percent_HC_TA_2 = [num2str(round(num_HC_TA_2 / (num_HC_TA_2 + num_NO_HC_TA_2) * 100)),'%'];

num_ST_TA_2= sum(strcmp(R_TA_2, 'ST'));
num_NO_ST_TA_2 = sum(strcmp(R_TA_2, 'NO ST'));
percent_ST_TA_2 =[num2str(round(num_ST_TA_2 / (num_ST_TA_2 + num_NO_ST_TA_2) * 100)),'%'];

num_TER_TA_2 = sum(strcmp(R_TA_2, 'TER'));
num_NO_TER_TA_2 = sum(strcmp(R_TA_2, 'NO TER'));
percent_TER_TA_2 = [num2str(round(num_TER_TA_2/ (num_TER_TA_2 + num_NO_TER_TA_2) * 100)),'%'];

num_SW_TA_2 = sum(strcmp(R_TA_2, 'SW'));
num_NO_SW_TA_2 = sum(strcmp(R_TA_2, 'NO SW'));
percent_SW_TA_2= [num2str(round(num_SW_TA_2 / (num_SW_TA_2 + num_NO_SW_TA_2) * 100)),'%'];
R_TA_PERC_2={percent_HC_TA_2;percent_ST_TA_2;percent_TER_TA_2;percent_SW_TA_2};

R_TA_new_3=[R_TA_tot(3,:)',results_phases_3(:,1)];
R_TA_3={};

for riga = 1:size(R_TA_new_3, 1)
    if strcmp(strtrim(R_TA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_3{riga,2}), 'HEEL CONTACT')
        R_TA_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_3{riga,2}), 'HEEL CONTACT')
        R_TA_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_3{riga,2}), 'SINGLE STANCE')
        R_TA_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_3{riga,2}), 'SINGLE STANCE')
        R_TA_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_3{riga,2}), 'TERMINAL STANCE')
        R_TA_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_3{riga,2}), 'TERMINAL STANCE')
        R_TA_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_TA_new_3{riga,2}), 'SWING')
        R_TA_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_TA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_TA_new_3{riga,2}), 'SWING')
        R_TA_3{end+1} = 'NO SW' ;
        end
end

num_HC_TA_3 = sum(strcmp(R_TA_3, 'HC'));
num_NO_HC_TA_3 = sum(strcmp(R_TA_3, 'NO HC'));
percent_HC_TA_3 = [num2str(round(num_HC_TA_3 / (num_HC_TA_3 + num_NO_HC_TA_3) * 100)),'%'];


num_ST_TA_3= sum(strcmp(R_TA_3, 'ST'));
num_NO_ST_TA_3 = sum(strcmp(R_TA_3, 'NO ST'));
percent_ST_TA_3 =[num2str(round(num_ST_TA_3 / (num_ST_TA_3 + num_NO_ST_TA_3) * 100)),'%'];

num_TER_TA_3 = sum(strcmp(R_TA_3, 'TER'));
num_NO_TER_TA_3 = sum(strcmp(R_TA_3, 'NO TER'));
percent_TER_TA_3 = [num2str(round(num_TER_TA_3/ (num_TER_TA_3 + num_NO_TER_TA_3) * 100)),'%'];


num_SW_TA_3 = sum(strcmp(R_TA_3, 'SW'));
num_NO_SW_TA_3 = sum(strcmp(R_TA_3, 'NO SW'));
percent_SW_TA_3= [num2str(round(num_SW_TA_3 / (num_SW_TA_3 + num_NO_SW_TA_3) * 100)),'%'];
R_TA_PERC_3={percent_HC_TA_3;percent_ST_TA_3;percent_TER_TA_3;percent_SW_TA_3};

R_TA_PERC=[R_TA_PERC_1;R_TA_PERC_2;R_TA_PERC_3];

%% R_PL
R_PL_new_1=[R_PL_tot(1,:)',results_phases_1(:,1)];

R_PL_1={};

for riga = 1:size(R_PL_new_1, 1)
    if strcmp(strtrim(R_PL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_1{riga,2}), 'HEEL CONTACT')
        R_PL_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_1{riga,2}), 'HEEL CONTACT')
        R_PL_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_1{riga,2}), 'SINGLE STANCE')
        R_PL_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_1{riga,2}), 'SINGLE STANCE')
        R_PL_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_1{riga,2}), 'TERMINAL STANCE')
        R_PL_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_1{riga,2}), 'TERMINAL STANCE')
        R_PL_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_1{riga,2}), 'SWING')
        R_PL_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_PL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_1{riga,2}), 'SWING')
        R_PL_1{end+1} = 'NO SW' ;
        end
end

num_HC_PL_1 = sum(strcmp(R_PL_1, 'HC'));
num_NO_HC_PL_1 = sum(strcmp(R_PL_1, 'NO HC'));
percent_HC_PL_1 = [num2str(round(num_HC_PL_1 / (num_HC_PL_1 + num_NO_HC_PL_1) * 100)),'%'];

num_ST_PL_1 = sum(strcmp(R_PL_1, 'ST'));
num_NO_PL_1 = sum(strcmp(R_PL_1, 'NO ST'));
percent_ST_PL_1 =[num2str(round(num_ST_PL_1 / (num_ST_PL_1 + num_NO_PL_1) * 100)),'%'];

num_TER_PL_1 = sum(strcmp(R_PL_1, 'TER'));
num_NO_TER_PL_1 = sum(strcmp(R_PL_1, 'NO TER'));
percent_TER_PL_1 = [num2str(round(num_TER_PL_1/ (num_TER_PL_1 + num_NO_TER_PL_1) * 100)),'%'];

num_SW_PL_1 = sum(strcmp(R_PL_1, 'SW'));
num_NO_SW_PL_1 = sum(strcmp(R_PL_1, 'NO SW'));
percent_SW_PL_1= [num2str(round(num_SW_PL_1 / (num_SW_PL_1 + num_NO_SW_PL_1) * 100)),'%'];
R_PL_PERC_1={percent_HC_PL_1;percent_ST_PL_1;percent_TER_PL_1;percent_SW_PL_1};


R_PL_new_2=[R_PL_tot(2,:)',results_phases_2(:,1)];
R_PL_2={};

for riga = 1:size(R_PL_new_2, 1)
    if strcmp(strtrim(R_PL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_2{riga,2}), 'HEEL CONTACT')
        R_PL_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_2{riga,2}), 'HEEL CONTACT')
        R_PL_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_2{riga,2}), 'SINGLE STANCE')
        R_PL_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_2{riga,2}), 'SINGLE STANCE')
        R_PL_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_2{riga,2}), 'TERMINAL STANCE')
        R_PL_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_2{riga,2}), 'TERMINAL STANCE')
        R_PL_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_2{riga,2}), 'SWING')
        R_PL_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_PL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_2{riga,2}), 'SWING')
        R_PL_2{end+1} = 'NO SW' ;
        end
end

num_HC_PL_2 = sum(strcmp(R_PL_2, 'HC'));
num_NO_HC_PL_2 = sum(strcmp(R_PL_2, 'NO HC'));
percent_HC_PL_2 = [num2str(round(num_HC_PL_2 / (num_HC_PL_2 + num_NO_HC_PL_2) * 100)),'%'];

num_ST_PL_2= sum(strcmp(R_PL_2, 'ST'));
num_NO_ST_PL_2 = sum(strcmp(R_PL_2, 'NO ST'));
percent_ST_PL_2 =[num2str(round(num_ST_PL_2 / (num_ST_PL_2 + num_NO_ST_PL_2) * 100)),'%'];

num_TER_PL_2 = sum(strcmp(R_PL_2, 'TER'));
num_NO_TER_PL_2 = sum(strcmp(R_PL_2, 'NO TER'));
percent_TER_PL_2 = [num2str(round(num_TER_PL_2/ (num_TER_PL_2 + num_NO_TER_PL_2) * 100)),'%'];

num_SW_PL_2 = sum(strcmp(R_PL_2, 'SW'));
num_NO_SW_PL_2 = sum(strcmp(R_PL_2, 'NO SW'));
percent_SW_PL_2= [num2str(round(num_SW_PL_2 / (num_SW_PL_2 + num_NO_SW_PL_2) * 100)),'%'];
R_PL_PERC_2={percent_HC_PL_2;percent_ST_PL_2;percent_TER_PL_2;percent_SW_PL_2};

R_PL_new_3=[R_PL_tot(3,:)',results_phases_3(:,1)];
R_PL_3={};

for riga = 1:size(R_PL_new_3, 1)
    if strcmp(strtrim(R_PL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_3{riga,2}), 'HEEL CONTACT')
        R_PL_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_3{riga,2}), 'HEEL CONTACT')
        R_PL_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_3{riga,2}), 'SINGLE STANCE')
        R_PL_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_3{riga,2}), 'SINGLE STANCE')
        R_PL_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_3{riga,2}), 'TERMINAL STANCE')
        R_PL_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_3{riga,2}), 'TERMINAL STANCE')
        R_PL_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_PL_new_3{riga,2}), 'SWING')
        R_PL_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_PL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_PL_new_3{riga,2}), 'SWING')
        R_PL_3{end+1} = 'NO SW' ;
        end
end

num_HC_PL_3 = sum(strcmp(R_PL_3, 'HC'));
num_NO_HC_PL_3 = sum(strcmp(R_PL_3, 'NO HC'));
percent_HC_PL_3 = [num2str(round(num_HC_PL_3 / (num_HC_PL_3 + num_NO_HC_PL_3) * 100)),'%'];


num_ST_PL_3= sum(strcmp(R_PL_3, 'ST'));
num_NO_ST_PL_3 = sum(strcmp(R_PL_3, 'NO ST'));
percent_ST_PL_3 =[num2str(round(num_ST_PL_3 / (num_ST_PL_3 + num_NO_ST_PL_3) * 100)),'%'];

num_TER_PL_3 = sum(strcmp(R_PL_3, 'TER'));
num_NO_TER_PL_3 = sum(strcmp(R_PL_3, 'NO TER'));
percent_TER_PL_3 = [num2str(round(num_TER_PL_3/ (num_TER_PL_3 + num_NO_TER_PL_3) * 100)),'%'];


num_SW_PL_3 = sum(strcmp(R_PL_3, 'SW'));
num_NO_SW_PL_3 = sum(strcmp(R_PL_3, 'NO SW'));
percent_SW_PL_3= [num2str(round(num_SW_PL_3 / (num_SW_PL_3 + num_NO_SW_PL_3) * 100)),'%'];
R_PL_PERC_3={percent_HC_PL_3;percent_ST_PL_3;percent_TER_PL_3;percent_SW_PL_3};

R_PL_PERC=[R_PL_PERC_1;R_PL_PERC_2;R_PL_PERC_3];

%% R_SO
R_SO_new_1=[R_SO_tot(1,:)',results_phases_1(:,1)];

R_SO_1={};

for riga = 1:size(R_SO_new_1, 1)
    if strcmp(strtrim(R_SO_new_1{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_1{riga,2}), 'HEEL CONTACT')
        R_SO_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_1{riga,2}), 'HEEL CONTACT')
        R_SO_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_1{riga,2}), 'SINGLE STANCE')
        R_SO_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_1{riga,2}), 'SINGLE STANCE')
        R_SO_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_1{riga,2}), 'TERMINAL STANCE')
        R_SO_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_1{riga,2}), 'TERMINAL STANCE')
        R_SO_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_1{riga,2}), 'SWING')
        R_SO_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_SO_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_1{riga,2}), 'SWING')
        R_SO_1{end+1} = 'NO SW' ;
        end
end

num_HC_SO_1 = sum(strcmp(R_SO_1, 'HC'));
num_NO_HC_SO_1 = sum(strcmp(R_SO_1, 'NO HC'));
percent_HC_SO_1 = [num2str(round(num_HC_SO_1 / (num_HC_SO_1 + num_NO_HC_SO_1) * 100)),'%'];

num_ST_SO_1 = sum(strcmp(R_SO_1, 'ST'));
num_NO_SO_1 = sum(strcmp(R_SO_1, 'NO ST'));
percent_ST_SO_1 =[num2str(round(num_ST_SO_1 / (num_ST_SO_1 + num_NO_SO_1) * 100)),'%'];

num_TER_SO_1 = sum(strcmp(R_SO_1, 'TER'));
num_NO_TER_SO_1 = sum(strcmp(R_SO_1, 'NO TER'));
percent_TER_SO_1 = [num2str(round(num_TER_SO_1/ (num_TER_SO_1 + num_NO_TER_SO_1) * 100)),'%'];

num_SW_SO_1 = sum(strcmp(R_SO_1, 'SW'));
num_NO_SW_SO_1 = sum(strcmp(R_SO_1, 'NO SW'));
percent_SW_SO_1= [num2str(round(num_SW_SO_1 / (num_SW_SO_1 + num_NO_SW_SO_1) * 100)),'%'];
R_SO_PERC_1={percent_HC_SO_1;percent_ST_SO_1;percent_TER_SO_1;percent_SW_SO_1};


R_SO_new_2=[R_SO_tot(2,:)',results_phases_2(:,1)];
R_SO_2={};

for riga = 1:size(R_SO_new_2, 1)
    if strcmp(strtrim(R_SO_new_2{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_2{riga,2}), 'HEEL CONTACT')
        R_SO_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_2{riga,2}), 'HEEL CONTACT')
        R_SO_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_2{riga,2}), 'SINGLE STANCE')
        R_SO_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_2{riga,2}), 'SINGLE STANCE')
        R_SO_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_2{riga,2}), 'TERMINAL STANCE')
        R_SO_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_2{riga,2}), 'TERMINAL STANCE')
        R_SO_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_2{riga,2}), 'SWING')
        R_SO_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_SO_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_2{riga,2}), 'SWING')
        R_SO_2{end+1} = 'NO SW' ;
        end
end

num_HC_SO_2 = sum(strcmp(R_SO_2, 'HC'));
num_NO_HC_SO_2 = sum(strcmp(R_SO_2, 'NO HC'));
percent_HC_SO_2 = [num2str(round(num_HC_SO_2 / (num_HC_SO_2 + num_NO_HC_SO_2) * 100)),'%'];

num_ST_SO_2= sum(strcmp(R_SO_2, 'ST'));
num_NO_ST_SO_2 = sum(strcmp(R_SO_2, 'NO ST'));
percent_ST_SO_2 =[num2str(round(num_ST_SO_2 / (num_ST_SO_2 + num_NO_ST_SO_2) * 100)),'%'];

num_TER_SO_2 = sum(strcmp(R_SO_2, 'TER'));
num_NO_TER_SO_2 = sum(strcmp(R_SO_2, 'NO TER'));
percent_TER_SO_2 = [num2str(round(num_TER_SO_2/ (num_TER_SO_2 + num_NO_TER_SO_2) * 100)),'%'];

num_SW_SO_2 = sum(strcmp(R_SO_2, 'SW'));
num_NO_SW_SO_2 = sum(strcmp(R_SO_2, 'NO SW'));
percent_SW_SO_2= [num2str(round(num_SW_SO_2 / (num_SW_SO_2 + num_NO_SW_SO_2) * 100)),'%'];
R_SO_PERC_2={percent_HC_SO_2;percent_ST_SO_2;percent_TER_SO_2;percent_SW_SO_2};

R_SO_new_3=[R_SO_tot(3,:)',results_phases_3(:,1)];
R_SO_3={};

for riga = 1:size(R_SO_new_3, 1)
    if strcmp(strtrim(R_SO_new_3{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_3{riga,2}), 'HEEL CONTACT')
        R_SO_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_3{riga,2}), 'HEEL CONTACT')
        R_SO_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_3{riga,2}), 'SINGLE STANCE')
        R_SO_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_3{riga,2}), 'SINGLE STANCE')
        R_SO_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_3{riga,2}), 'TERMINAL STANCE')
        R_SO_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_3{riga,2}), 'TERMINAL STANCE')
        R_SO_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'ON') && strcmp(strtrim(R_SO_new_3{riga,2}), 'SWING')
        R_SO_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_SO_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_SO_new_3{riga,2}), 'SWING')
        R_SO_3{end+1} = 'NO SW' ;
        end
end

num_HC_SO_3 = sum(strcmp(R_SO_3, 'HC'));
num_NO_HC_SO_3 = sum(strcmp(R_SO_3, 'NO HC'));
percent_HC_SO_3 = [num2str(round(num_HC_SO_3 / (num_HC_SO_3 + num_NO_HC_SO_3) * 100)),'%'];


num_ST_SO_3= sum(strcmp(R_SO_3, 'ST'));
num_NO_ST_SO_3 = sum(strcmp(R_SO_3, 'NO ST'));
percent_ST_SO_3 =[num2str(round(num_ST_SO_3 / (num_ST_SO_3 + num_NO_ST_SO_3) * 100)),'%'];

num_TER_SO_3 = sum(strcmp(R_SO_3, 'TER'));
num_NO_TER_SO_3 = sum(strcmp(R_SO_3, 'NO TER'));
percent_TER_SO_3 = [num2str(round(num_TER_SO_3/ (num_TER_SO_3 + num_NO_TER_SO_3) * 100)),'%'];


num_SW_SO_3 = sum(strcmp(R_SO_3, 'SW'));
num_NO_SW_SO_3 = sum(strcmp(R_SO_3, 'NO SW'));
percent_SW_SO_3= [num2str(round(num_SW_SO_3 / (num_SW_SO_3 + num_NO_SW_SO_3) * 100)),'%'];
R_SO_PERC_3={percent_HC_SO_3;percent_ST_SO_3;percent_TER_SO_3;percent_SW_SO_3};

R_SO_PERC=[R_SO_PERC_1;R_SO_PERC_2;R_SO_PERC_3];

%% R_GL
R_GL_new_1=[R_GL_tot(1,:)',results_phases_1(:,1)];

R_GL_1={};

for riga = 1:size(R_GL_new_1, 1)
    if strcmp(strtrim(R_GL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_1{riga,2}), 'HEEL CONTACT')
        R_GL_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_1{riga,2}), 'HEEL CONTACT')
        R_GL_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_1{riga,2}), 'SINGLE STANCE')
        R_GL_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_1{riga,2}), 'SINGLE STANCE')
        R_GL_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_1{riga,2}), 'TERMINAL STANCE')
        R_GL_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_1{riga,2}), 'TERMINAL STANCE')
        R_GL_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_1{riga,2}), 'SWING')
        R_GL_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_GL_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_1{riga,2}), 'SWING')
        R_GL_1{end+1} = 'NO SW' ;
        end
end

num_HC_GL_1 = sum(strcmp(R_GL_1, 'HC'));
num_NO_HC_GL_1 = sum(strcmp(R_GL_1, 'NO HC'));
percent_HC_GL_1 = [num2str(round(num_HC_GL_1 / (num_HC_GL_1 + num_NO_HC_GL_1) * 100)),'%'];

num_ST_GL_1 = sum(strcmp(R_GL_1, 'ST'));
num_NO_GL_1 = sum(strcmp(R_GL_1, 'NO ST'));
percent_ST_GL_1 =[num2str(round(num_ST_GL_1 / (num_ST_GL_1 + num_NO_GL_1) * 100)),'%'];

num_TER_GL_1 = sum(strcmp(R_GL_1, 'TER'));
num_NO_TER_GL_1 = sum(strcmp(R_GL_1, 'NO TER'));
percent_TER_GL_1 = [num2str(round(num_TER_GL_1/ (num_TER_GL_1 + num_NO_TER_GL_1) * 100)),'%'];

num_SW_GL_1 = sum(strcmp(R_GL_1, 'SW'));
num_NO_SW_GL_1 = sum(strcmp(R_GL_1, 'NO SW'));
percent_SW_GL_1= [num2str(round(num_SW_GL_1 / (num_SW_GL_1 + num_NO_SW_GL_1) * 100)),'%'];
R_GL_PERC_1={percent_HC_GL_1;percent_ST_GL_1;percent_TER_GL_1;percent_SW_GL_1};


R_GL_new_2=[R_GL_tot(2,:)',results_phases_2(:,1)];
R_GL_2={};

for riga = 1:size(R_GL_new_2, 1)
    if strcmp(strtrim(R_GL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_2{riga,2}), 'HEEL CONTACT')
        R_GL_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_2{riga,2}), 'HEEL CONTACT')
        R_GL_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_2{riga,2}), 'SINGLE STANCE')
        R_GL_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_2{riga,2}), 'SINGLE STANCE')
        R_GL_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_2{riga,2}), 'TERMINAL STANCE')
        R_GL_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_2{riga,2}), 'TERMINAL STANCE')
        R_GL_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_2{riga,2}), 'SWING')
        R_GL_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_GL_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_2{riga,2}), 'SWING')
        R_GL_2{end+1} = 'NO SW' ;
        end
end

num_HC_GL_2 = sum(strcmp(R_GL_2, 'HC'));
num_NO_HC_GL_2 = sum(strcmp(R_GL_2, 'NO HC'));
percent_HC_GL_2 = [num2str(round(num_HC_GL_2 / (num_HC_GL_2 + num_NO_HC_GL_2) * 100)),'%'];

num_ST_GL_2= sum(strcmp(R_GL_2, 'ST'));
num_NO_ST_GL_2 = sum(strcmp(R_GL_2, 'NO ST'));
percent_ST_GL_2 =[num2str(round(num_ST_GL_2 / (num_ST_GL_2 + num_NO_ST_GL_2) * 100)),'%'];

num_TER_GL_2 = sum(strcmp(R_GL_2, 'TER'));
num_NO_TER_GL_2 = sum(strcmp(R_GL_2, 'NO TER'));
percent_TER_GL_2 = [num2str(round(num_TER_GL_2/ (num_TER_GL_2 + num_NO_TER_GL_2) * 100)),'%'];

num_SW_GL_2 = sum(strcmp(R_GL_2, 'SW'));
num_NO_SW_GL_2 = sum(strcmp(R_GL_2, 'NO SW'));
percent_SW_GL_2= [num2str(round(num_SW_GL_2 / (num_SW_GL_2 + num_NO_SW_GL_2) * 100)),'%'];
R_GL_PERC_2={percent_HC_GL_2;percent_ST_GL_2;percent_TER_GL_2;percent_SW_GL_2};

R_GL_new_3=[R_GL_tot(3,:)',results_phases_3(:,1)];
R_GL_3={};

for riga = 1:size(R_GL_new_3, 1)
    if strcmp(strtrim(R_GL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_3{riga,2}), 'HEEL CONTACT')
        R_GL_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_3{riga,2}), 'HEEL CONTACT')
        R_GL_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_3{riga,2}), 'SINGLE STANCE')
        R_GL_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_3{riga,2}), 'SINGLE STANCE')
        R_GL_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_3{riga,2}), 'TERMINAL STANCE')
        R_GL_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_3{riga,2}), 'TERMINAL STANCE')
        R_GL_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GL_new_3{riga,2}), 'SWING')
        R_GL_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_GL_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GL_new_3{riga,2}), 'SWING')
        R_GL_3{end+1} = 'NO SW' ;
        end
end

num_HC_GL_3 = sum(strcmp(R_GL_3, 'HC'));
num_NO_HC_GL_3 = sum(strcmp(R_GL_3, 'NO HC'));
percent_HC_GL_3 = [num2str(round(num_HC_GL_3 / (num_HC_GL_3 + num_NO_HC_GL_3) * 100)),'%'];


num_ST_GL_3= sum(strcmp(R_GL_3, 'ST'));
num_NO_ST_GL_3 = sum(strcmp(R_GL_3, 'NO ST'));
percent_ST_GL_3 =[num2str(round(num_ST_GL_3 / (num_ST_GL_3 + num_NO_ST_GL_3) * 100)),'%'];

num_TER_GL_3 = sum(strcmp(R_GL_3, 'TER'));
num_NO_TER_GL_3 = sum(strcmp(R_GL_3, 'NO TER'));
percent_TER_GL_3 = [num2str(round(num_TER_GL_3/ (num_TER_GL_3 + num_NO_TER_GL_3) * 100)),'%'];


num_SW_GL_3 = sum(strcmp(R_GL_3, 'SW'));
num_NO_SW_GL_3 = sum(strcmp(R_GL_3, 'NO SW'));
percent_SW_GL_3= [num2str(round(num_SW_GL_3 / (num_SW_GL_3 + num_NO_SW_GL_3) * 100)),'%'];
R_GL_PERC_3={percent_HC_GL_3;percent_ST_GL_3;percent_TER_GL_3;percent_SW_GL_3};

R_GL_PERC=[R_GL_PERC_1;R_GL_PERC_2;R_GL_PERC_3];

%% R_GMA
R_GMA_new_1=[R_GMA_tot(1,:)',results_phases_1(:,1)];

R_GMA_1={};

for riga = 1:size(R_GMA_new_1, 1)
    if strcmp(strtrim(R_GMA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'HEEL CONTACT')
        R_GMA_1{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'HEEL CONTACT')
        R_GMA_1{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'SINGLE STANCE')
        R_GMA_1{end+1} = 'ST';
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'SINGLE STANCE')
        R_GMA_1{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'TERMINAL STANCE')
        R_GMA_1{end+1} = 'TER';
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'TERMINAL STANCE')
        R_GMA_1{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'SWING')
        R_GMA_1{end+1} = 'SW';
    elseif strcmp(strtrim(R_GMA_new_1{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_1{riga,2}), 'SWING')
        R_GMA_1{end+1} = 'NO SW' ;
        end
end

num_HC_GMA_1 = sum(strcmp(R_GMA_1, 'HC'));
num_NO_HC_GMA_1 = sum(strcmp(R_GMA_1, 'NO HC'));
percent_HC_GMA_1 = [num2str(round(num_HC_GMA_1 / (num_HC_GMA_1 + num_NO_HC_GMA_1) * 100)),'%'];

num_ST_GMA_1 = sum(strcmp(R_GMA_1, 'ST'));
num_NO_GMA_1 = sum(strcmp(R_GMA_1, 'NO ST'));
percent_ST_GMA_1 =[num2str(round(num_ST_GMA_1 / (num_ST_GMA_1 + num_NO_GMA_1) * 100)),'%'];

num_TER_GMA_1 = sum(strcmp(R_GMA_1, 'TER'));
num_NO_TER_GMA_1 = sum(strcmp(R_GMA_1, 'NO TER'));
percent_TER_GMA_1 = [num2str(round(num_TER_GMA_1/ (num_TER_GMA_1 + num_NO_TER_GMA_1) * 100)),'%'];

num_SW_GMA_1 = sum(strcmp(R_GMA_1, 'SW'));
num_NO_SW_GMA_1 = sum(strcmp(R_GMA_1, 'NO SW'));
percent_SW_GMA_1= [num2str(round(num_SW_GMA_1 / (num_SW_GMA_1 + num_NO_SW_GMA_1) * 100)),'%'];
R_GMA_PERC_1={percent_HC_GMA_1;percent_ST_GMA_1;percent_TER_GMA_1;percent_SW_GMA_1};


R_GMA_new_2=[R_GMA_tot(2,:)',results_phases_2(:,1)];
R_GMA_2={};

for riga = 1:size(R_GMA_new_2, 1)
    if strcmp(strtrim(R_GMA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'HEEL CONTACT')
        R_GMA_2{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'HEEL CONTACT')
        R_GMA_2{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'SINGLE STANCE')
        R_GMA_2{end+1} = 'ST';
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'SINGLE STANCE')
        R_GMA_2{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'TERMINAL STANCE')
        R_GMA_2{end+1} = 'TER';
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'TERMINAL STANCE')
        R_GMA_2{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'SWING')
        R_GMA_2{end+1} = 'SW';
    elseif strcmp(strtrim(R_GMA_new_2{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_2{riga,2}), 'SWING')
        R_GMA_2{end+1} = 'NO SW' ;
        end
end

num_HC_GMA_2 = sum(strcmp(R_GMA_2, 'HC'));
num_NO_HC_GMA_2 = sum(strcmp(R_GMA_2, 'NO HC'));
percent_HC_GMA_2 = [num2str(round(num_HC_GMA_2 / (num_HC_GMA_2 + num_NO_HC_GMA_2) * 100)),'%'];

num_ST_GMA_2= sum(strcmp(R_GMA_2, 'ST'));
num_NO_ST_GMA_2 = sum(strcmp(R_GMA_2, 'NO ST'));
percent_ST_GMA_2 =[num2str(round(num_ST_GMA_2 / (num_ST_GMA_2 + num_NO_ST_GMA_2) * 100)),'%'];

num_TER_GMA_2 = sum(strcmp(R_GMA_2, 'TER'));
num_NO_TER_GMA_2 = sum(strcmp(R_GMA_2, 'NO TER'));
percent_TER_GMA_2 = [num2str(round(num_TER_GMA_2/ (num_TER_GMA_2 + num_NO_TER_GMA_2) * 100)),'%'];

num_SW_GMA_2 = sum(strcmp(R_GMA_2, 'SW'));
num_NO_SW_GMA_2 = sum(strcmp(R_GMA_2, 'NO SW'));
percent_SW_GMA_2= [num2str(round(num_SW_GMA_2 / (num_SW_GMA_2 + num_NO_SW_GMA_2) * 100)),'%'];
R_GMA_PERC_2={percent_HC_GMA_2;percent_ST_GMA_2;percent_TER_GMA_2;percent_SW_GMA_2};

R_GMA_new_3=[R_GMA_tot(3,:)',results_phases_3(:,1)];
R_GMA_3={};

for riga = 1:size(R_GMA_new_3, 1)
    if strcmp(strtrim(R_GMA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'HEEL CONTACT')
        R_GMA_3{end+1} = 'HC';       
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'HEEL CONTACT')
        R_GMA_3{end+1} = 'NO HC' ;
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'SINGLE STANCE')
        R_GMA_3{end+1} = 'ST';
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'SINGLE STANCE')
        R_GMA_3{end+1} = 'NO ST' ;
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'TERMINAL STANCE')
        R_GMA_3{end+1} = 'TER';
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'TERMINAL STANCE')
        R_GMA_3{end+1} = 'NO TER' ;
        elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'ON') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'SWING')
        R_GMA_3{end+1} = 'SW';
    elseif strcmp(strtrim(R_GMA_new_3{riga,1}), 'OFF') && strcmp(strtrim(R_GMA_new_3{riga,2}), 'SWING')
        R_GMA_3{end+1} = 'NO SW' ;
        end
end

num_HC_GMA_3 = sum(strcmp(R_GMA_3, 'HC'));
num_NO_HC_GMA_3 = sum(strcmp(R_GMA_3, 'NO HC'));
percent_HC_GMA_3 = [num2str(round(num_HC_GMA_3 / (num_HC_GMA_3 + num_NO_HC_GMA_3) * 100)),'%'];


num_ST_GMA_3= sum(strcmp(R_GMA_3, 'ST'));
num_NO_ST_GMA_3 = sum(strcmp(R_GMA_3, 'NO ST'));
percent_ST_GMA_3 =[num2str(round(num_ST_GMA_3 / (num_ST_GMA_3 + num_NO_ST_GMA_3) * 100)),'%'];

num_TER_GMA_3 = sum(strcmp(R_GMA_3, 'TER'));
num_NO_TER_GMA_3 = sum(strcmp(R_GMA_3, 'NO TER'));
percent_TER_GMA_3 = [num2str(round(num_TER_GMA_3/ (num_TER_GMA_3 + num_NO_TER_GMA_3) * 100)),'%'];


num_SW_GMA_3 = sum(strcmp(R_GMA_3, 'SW'));
num_NO_SW_GMA_3 = sum(strcmp(R_GMA_3, 'NO SW'));
percent_SW_GMA_3= [num2str(round(num_SW_GMA_3 / (num_SW_GMA_3 + num_NO_SW_GMA_3) * 100)),'%'];
R_GMA_PERC_3={percent_HC_GMA_3;percent_ST_GMA_3;percent_TER_GMA_3;percent_SW_GMA_3};

R_GMA_PERC=[R_GMA_PERC_1;R_GMA_PERC_2;R_GMA_PERC_3];



% Phases column
results_phases={'HEEL CONTACT';'SINGLE STANCE';'TERMINAL STANCE';'SWING';'HEEL CONTACT';'SINGLE STANCE';'TERMINAL STANCE';'SWING';'HEEL CONTACT';'SINGLE STANCE';'TERMINAL STANCE';'SWING'};

%% FINAL MATRIX 

results_matrix=[results_ID, results_WALK,results_Patology,results_phases,R_RF_PERC,R_VL_PERC,R_VM_PERC,R_G_PERC,R_BFCL_PERC,R_TA_PERC,R_PL_PERC,R_SO_PERC,R_GL_PERC,R_GMA_PERC];


output_folder=uigetdir('','Select destination folder CSV'); 
if output_folder==0
    disp('Operation canceled.')%if the folder is empty
else
    input_user=inputdlg('Enter a name for the CSV file to save:','Save as');

    if isempty(input_user)
    disp('Operation canceled.')%if the user clicks on 'Annulla'
else
    file_name=input_user{1};
    file_path=fullfile(output_folder,[file_name,'.csv']);

writematrix(results_matrix, file_path);

    end
end














