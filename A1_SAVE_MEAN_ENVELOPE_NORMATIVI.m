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
Velocity=answer{3};%Velocity trial: M or L for NORMATIVO and blank space for ATASSIA or PARAPARESI
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

%%plot
% figure(1); 
% plot(time_acq,R_RF_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right rectus femoris','FontWeight','bold')
% title('EMG Right Rectus Femoris: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_RF_cen , 'Color','b');
% hold on
% plot(time_acq, R_RF_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_RF_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_RF_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified', 'EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure
% plot(time_acq,R_RF_raw,'g')
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% xlabel('Time','FontWeight','bold')
% ylabel('Right Rectus Femoris','FontWeight','bold')
% title('EMG Right Rectus Femoris: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_RF_cen , 'Color','b');
% hold on
% plot(time_acq, R_RF_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_RF_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_RF_env , 'Color' , 'k');
% hold on
% 
% legend('EMG raw','EMG centered','EMG filtered', 'EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(2); 
% plot(time_acq,R_VL_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right vastus lateralis','FontWeight','bold')
% title('EMG Right Vastus Lateralis: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_VL_cen , 'Color','b');
% hold on
% plot(time_acq, R_VL_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_VL_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_VL_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2, time_foot_strike_2], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified', 'EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure
% plot(time_acq,R_VL_raw,'g')
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% xlabel('Time','FontWeight','bold')
% ylabel('Right Vastus Lateralis','FontWeight','bold')
% title('EMG Right Vastus Lateralis: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_VL_cen , 'Color','b');
% hold on
% plot(time_acq, R_VL_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_VL_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_VL_env , 'Color' , 'k');
% hold on
% 
% legend('EMG raw','EMG centered','EMG filtered', 'EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(3); 
% plot(time_acq,R_VM_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right vastus medialis','FontWeight','bold')
% title('EMG Right Vastus Medialis: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_VM_cen , 'Color','b');
% hold on
% plot(time_acq, R_VM_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_VM_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_VM_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified', 'EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window gait cycle
% figure
% plot(time_acq,R_VM_raw,'g')
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% xlabel('Time','FontWeight','bold')
% ylabel('Right vastus medialis','FontWeight','bold')
% title('EMG Right Vastus Medialis: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_VM_cen , 'Color','b');
% hold on
% plot(time_acq, R_VM_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_VM_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_VM_env , 'Color' , 'k');
% hold on
% 
% legend('EMG raw','EMG centered','EMG filtered', 'EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(4); 
% plot(time_acq,R_G_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right gluteus medius','FontWeight','bold')
% title('EMG Right Gluteus Medius: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_G_cen , 'Color','b');
% hold on
% plot(time_acq, R_G_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_G_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_G_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified', 'EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure
% plot(time_acq,R_G_raw,'g')
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% xlabel('Time','FontWeight','bold')
% ylabel('Right Gluteus Medius','FontWeight','bold')
% title('EMG Right Gluteus Medius: gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_G_cen , 'Color','b');
% hold on
% plot(time_acq, R_G_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_G_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_G_env , 'Color' , 'k');
% hold on
% 
% legend('EMG raw','EMG centered','EMG filtered', 'EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(7); 
% hold on; 
% plot(time_acq,R_BFCL_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right biceps femoris caput longus','FontWeight','bold')
% title('EMG Right Biceps Femoris Caput Longus:gait cycle','FontSize',15)
% plot(time_acq, R_BFCL_cen , 'Color','b');
% hold on
% plot(time_acq, R_BFCL_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_BFCL_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_BFCL_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered', 'EMG filtered','EMG rectified','EMG envelope', 'Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% plot(time_acq,R_BFCL_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% 
% xlabel('Time','FontWeight','bold')
% ylabel('Right biceps femoris caput longus','FontWeight','bold')
% title('EMG Right Biceps Femoris Caput Longus:gait cycle','FontSize',15)
% hold on
% plot(time_acq, R_BFCL_cen , 'Color','b');
% hold on
% plot(time_acq, R_BFCL_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_BFCL_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_BFCL_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered', 'EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(8); 
% hold on; 
% plot(time_acq,R_TA_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right tibialis anterior','FontWeight','bold')
% title('EMG Right Tibialis Anterior:gait cycle','FontSize',15)
% plot(time_acq, R_TA_cen , 'Color','b');
% hold on
% plot(time_acq, R_TA_filt , 'Color', 'r');
% hold on
% plot(time_acq,R_TA_abs , 'Color' , 'y');
% hold on
% plot(time_acq,R_TA_env , 'Color' , 'k');
% hold on
% 
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% hold on; 
% plot(time_acq,R_TA_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% 
% xlabel('Time','FontWeight','bold')
% ylabel('Right tibialis anterior','FontWeight','bold')
% title('EMG Right Tibialis Anterior:gait cycle','FontSize',15)
% plot(time_acq, R_TA_cen , 'Color','b');
% hold on
% plot(time_acq, R_TA_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_TA_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_TA_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Location','best')
% hold off

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

%%plot
% figure(9); 
% hold on; 
% plot(time_acq,R_PL_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right peroneus longus','FontWeight','bold')
% title('EMG Right Peroneus Longus:gait cycle','FontSize',15)
% plot(time_acq, R_PL_cen , 'Color','b');
% hold on
% plot(time_acq, R_PL_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_PL_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_PL_env , 'Color' , 'k');
% hold on
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% hold on; 
% plot(time_acq,R_PL_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% 
% xlabel('Time','FontWeight','bold')
% ylabel('Right peroneus longus','FontWeight','bold')
% title('EMG Right Peroneus Longus:gait cycle','FontSize',15)
% plot(time_acq, R_PL_cen , 'Color','b');
% hold on
% plot(time_acq, R_PL_filt , 'Color','r');
% hold on
% plot(time_acq, R_PL_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_PL_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Location','best')
% hold off

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

%%plot
% figure(10); 
% hold on; 
% plot(time_acq,R_SO_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right soleus','FontWeight','bold')
% title('EMG Right Soleus:gait cycle','FontSize',15)
% plot(time_acq, R_SO_cen , 'Color','b');
% hold on
% plot(time_acq, R_SO_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_SO_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_SO_env , 'Color' , 'k');
% hold on
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% hold on; 
% plot(time_acq,R_SO_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% 
% xlabel('Time','FontWeight','bold')
% ylabel('Right soleus','FontWeight','bold')
% title('EMG Right Soleus:gait cycle','FontSize',15)
% plot(time_acq, R_SO_cen , 'Color','b');
% hold on
% plot(time_acq, R_SO_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_SO_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_SO_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope','Location','best')
% hold off

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

%%plot
% figure(12); 
% hold on; 
% plot(time_acq,R_GL_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right gastrocnemius lateralis','FontWeight','bold')
% title('EMG Right Gastrocnemius Lateralis:gait cycle','FontSize',15)
% plot(time_acq, R_GL_cen , 'Color','b');
% hold on
% plot(time_acq, R_GL_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_GL_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_GL_env , 'Color' , 'k');
% hold on
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% hold on; 
% plot(time_acq,R_GL_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% 
% xlabel('Time','FontWeight','bold')
% ylabel('Right gastrocnemius lateralis','FontWeight','bold')
% title('EMG Right Gastrocnemius Lateralis:gait cycle','FontSize',15)
% plot(time_acq, R_GL_cen , 'Color','b');
% hold on
% plot(time_acq, R_GL_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_GL_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_GL_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope', 'Location','best')
% hold off

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

%%plot
% figure(13); 
% hold on; 
% plot(time_acq,R_GMA_raw , 'Color', 'g');
% xlim([0 12])
% xlabel('Time','FontWeight','bold')
% ylabel('Right gluteus maximus','FontWeight','bold')
% title('EMG Right Gluteus Maximus:gait cycle','FontSize',15)
% plot(time_acq, R_GMA_cen , 'Color','b');
% hold on
% plot(time_acq, R_GMA_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_GMA_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_GMA_env , 'Color' , 'k');
% hold on
% 
% %Gait cycle: plot
% line([time_foot_strike_1(i,1), time_foot_strike_1(i,1)], ylim, 'Color', 'g', 'LineStyle', '--'); %right foot strike start
% line([time_foot_strike_2(i,1), time_foot_strike_2(i,1)], ylim, 'Color', 'm', 'LineStyle', '--'); %right foot strike stop
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope','Start cycle','End cycle','Location','best')
% hold off
% 
% %Time window: gait cycle
% figure; 
% hold on; 
% plot(time_acq,R_GMA_raw , 'Color', 'g');
% xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
% xlabel('Time','FontWeight','bold')
% ylabel('Right gluteus maximus','FontWeight','bold')
% title('EMG Right Gluteus Maximus:gait cycle','FontSize',15)
% plot(time_acq, R_GMA_cen , 'Color','b');
% hold on
% plot(time_acq, R_GMA_filt , 'Color', 'r');
% hold on
% plot(time_acq, R_GMA_abs , 'Color', 'y');
% hold on
% plot(time_acq,R_GMA_env , 'Color' , 'k');
% legend('EMG raw','EMG centered','EMG filtered','EMG rectified','EMG envelope','Location','best')
% hold off

%% SAMPLE NUMBER CALCULATION

time_window_length(i,1)=(time_foot_strike_2(i,1)-time_foot_strike_1(i,1));
num_samples(i,1)=time_window_length(i,1)*frequency;

end

%Save data for the mean envelope between all the subjects

name_file_excel="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\NUM_SAMPLES_NORMATIVI.xlsx";
if exist(name_file_excel,'file')==2
    existent_data=readmatrix(name_file_excel);
    dati_to_save=[existent_data;num_samples];
else
    dati_to_save=num_samples;
end
writematrix(dati_to_save,name_file_excel);

%% 

max_percentage=0.2;

%% RIGHT RECTUS FEMORIS

% Offset,Normalization and threshold
max_global_RF = -Inf;  

for i = 1:length(A)
    max_local_RF = max(R_RF{i});
    if max_local_RF > max_global_RF
        max_global_RF = max_local_RF;
    end
end

% min_global_RF=1;
% for i = 1:length(A)
%     min_local_RF = min(R_RF{i});
%     if min_local_RF < min_global_RF
%         min_global_RF = min_local_RF;
%     end
% end

for i=1:length(A)
normalized_R_RF{i}=(R_RF{i})./(max_global_RF);
end

threshold_R_RF=max_percentage*max_global_RF;
threshold_R_RF_norm = (threshold_R_RF) / (max_global_RF);


for i=1:length(A)

figure
plot(time_vector{i},normalized_R_RF{i},'Color','b');
xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
xlabel('Time','FontWeight','bold')
ylabel('Right Rectus Femoris','FontWeight','bold')
title('EMG Right Rectus Femoris','FontSize',15)
hold on
line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_RF_norm,threshold_R_RF_norm], 'Color', 'r', 'LineStyle', '--'); 
hold off
end


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
% Offset,Normalization and threshold

min_global_VL=1;
for i = 1:length(A)
    min_local_VL = min(R_VL{i});
    if min_local_VL < min_global_VL
        min_global_VL = min_local_VL;
    end
end

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
threshold_R_VL_norm = (threshold_R_VL) / (max_global_VL);

for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_VL{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Vastus Lateralis','FontWeight','bold')
    title('EMG Right Vastus Lateralis','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_VL_norm,threshold_R_VL_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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
% Offset,Normalization and threshold

min_global_VM=1;
for i = 1:length(A)
    min_local_VM = min(R_VM{i});
    if min_local_VM < min_global_VM
        min_global_VM = min_local_VM;
    end
end

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


for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_VM{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Vastus Medialis','FontWeight','bold')
    title('EMG Right Vastus Medialis','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_VM_norm,threshold_R_VM_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end
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

% Offset,Normalization and threshold

min_global_G=1;
for i = 1:length(A)
    min_local_G = min(R_G{i});
    if min_local_G < min_global_G
        min_global_G = min_local_G;
    end
end

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
threshold_R_G_norm = (threshold_R_G ) / (max_global_G);

for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_G{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Gluteus Medius','FontWeight','bold')
    title('EMG Right Gluteus Medius','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_G_norm,threshold_R_G_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_BFCL=1;
for i = 1:length(A)
    min_local_BFCL = min(R_BFCL{i});
    if min_local_BFCL < min_global_BFCL
        min_global_BFCL = min_local_BFCL;
    end
end

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


for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_BFCL{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Biceps Femoris Caput Longus','FontWeight','bold')
    title('EMG Right Biceps Femoris Caput Longus','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_BFCL_norm,threshold_R_BFCL_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_TA=1;
for i = 1:length(A)
    min_local_TA = min(R_TA{i});
    if min_local_TA < min_global_TA
        min_global_TA = min_local_TA;
    end
end

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

for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_TA{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Tibialis Anterior','FontWeight','bold')
    title('EMG Right Tibialis Anterior','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_TA_norm,threshold_R_TA_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_PL=1;
for i = 1:length(A)
    min_local_PL = min(R_PL{i});
    if min_local_PL < min_global_PL
        min_global_PL = min_local_PL;
    end
end

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

for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_PL{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Peroneus Longus','FontWeight','bold')
    title('EMG Right Peroneus Longus','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_PL_norm,threshold_R_PL_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_SO=1;
for i = 1:length(A)
    min_local_SO = min(R_SO{i});
    if min_local_SO < min_global_SO
        min_global_SO = min_local_SO;
    end
end

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


for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_SO{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Soleus','FontWeight','bold')
    title('EMG Right Soleus','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_SO_norm,threshold_R_SO_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_GL=1;
for i = 1:length(A)
    min_local_GL = min(R_GL{i});
    if min_local_GL < min_global_GL
        min_global_GL = min_local_GL;
    end
end

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


for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_GL{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Gastrocnemius Lateralis','FontWeight','bold')
    title('EMG Right Gastrocnemius Lateralis','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_GL_norm,threshold_R_GL_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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

% Offset,Normalization and threshold

min_global_GMA=1;
for i = 1:length(A)
    min_local_GMA = min(R_GMA{i});
    if min_local_GMA < min_global_GMA
        min_global_GMA = min_local_GMA;
    end
end

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
threshold_R_GMA_norm = (threshold_R_GMA) / (max_global_GMA);


for i=1:length(A)
    figure
    plot(time_vector{i},normalized_R_GMA{i},'Color','b');
    xlim([time_foot_strike_1(i,1) time_foot_strike_2(i,1)])
    xlabel('Time','FontWeight','bold')
    ylabel('Right Gluteus Maximus','FontWeight','bold')
    title('EMG Right Gluteus Maximus','FontSize',15)
    hold on
    line([time_foot_strike_1(i,1) ,time_foot_strike_2(i,1)],[threshold_R_GMA_norm,threshold_R_GMA_norm], 'Color', 'r', 'LineStyle', '--'); 
    hold off
end

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
    % figure
    % plot(time_vector{i},R_RF_above_threshold{i},'o', xval{i}, R_RF_unif_uniformed{i})
    % title('Interpolation RF','FontSize',15)
    % legend('Original Data RF', 'Interpolating Polynomial'); 
    % 

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
    % figure
    % plot(time_vector{i},R_VL_above_threshold{i},'o', xval{i}, R_VL_unif_uniformed{i})
    % title('Interpolation VL','FontSize',15)
    % legend('Original Data VL', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_VM_above_threshold{i},'o', xval{i}, R_VM_unif_uniformed{i})
    % title('Interpolation VM','FontSize',15)
    % legend('Original Data VM', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_G_above_threshold{i},'o', xval{i}, R_G_unif_uniformed{i})
    % title('Interpolation G','FontSize',15)
    % legend('Original Data G', 'Interpolating Polynomial'); 

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
    R_BFCL_unif_uniformed{i} = max(0, R_BFCL_unif_uniformed{i});
    R_BFCL_unif_uniformed{i} = min(1, R_BFCL_unif_uniformed{i});
    % figure
    % plot(time_vector{i},R_BFCL_above_threshold{i},'o', xval{i}, R_BFCL_unif_uniformed{i})
    % title('Interpolation BFCL','FontSize',15)
    % legend('Original Data BFCL', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_TA_above_threshold{i},'o', xval{i}, R_TA_unif_uniformed{i})
    % title('Interpolation TA','FontSize',15)
    % legend('Original Data TA', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_PL_above_threshold{i},'o', xval{i}, R_PL_unif_uniformed{i})
    % title('Interpolation PL','FontSize',15)
    % legend('Original Data PL', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_SO_above_threshold{i},'o', xval{i}, R_SO_unif_uniformed{i})
    % title('Interpolation SO','FontSize',15)
    % legend('Original Data SO', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_GL_above_threshold{i},'o', xval{i}, R_GL_unif_uniformed{i})
    % title('Interpolation GL','FontSize',15)
    % legend('Original Data GL', 'Interpolating Polynomial'); 

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
    % figure
    % plot(time_vector{i},R_GMA_above_threshold{i},'o', xval{i}, R_GMA_unif_uniformed{i})
    % title('Interpolation GMA','FontSize',15)
    % legend('Original Data GMA', 'Interpolating Polynomial'); 

end

%Vector column with GMA results
results_GMA=[];
for j=1:numel(R_GMA_unif_uniformed)
    vector_GMA=R_GMA_unif_uniformed{j}';
    results_GMA=[results_GMA;vector_GMA];
end

%% SUBPLOT INTERPOLATION

for i=1:length(A)
figure
subplot(5,2,1)
generalTitle = [ID ' ' Velocity];
sgtitle(generalTitle)
title('Interpolation RF')
plot(time_vector{i}, R_RF_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_RF_unif_uniformed{i}, 'LineWidth', 1);
hold off;


subplot(5,2,2)
title('Interpolation VL')
plot(time_vector{i}, R_VL_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_VL_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,3)
title('Interpolation VM')
plot(time_vector{i}, R_VM_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_VM_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,4)
title('Interpolation G')
plot(time_vector{i}, R_G_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_G_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,5)
title('Interpolation BFCL')
plot(time_vector{i}, R_BFCL_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_BFCL_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,6)
title('Interpolation TA')
plot(time_vector{i}, R_TA_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_TA_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,7)
title('Interpolation PL')
plot(time_vector{i}, R_PL_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_PL_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,8)
title('Interpolation SO')
plot(time_vector{i}, R_SO_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_SO_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,9)
title('Interpolation GL')
plot(time_vector{i}, R_GL_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_GL_unif_uniformed{i}, 'LineWidth', 1);
hold off


subplot(5,2,10)
title('Interpolation GMA')
plot(time_vector{i}, R_GMA_above_threshold{i}, 'o', 'MarkerSize', 3, 'LineWidth', 0.5);
hold on;
plot(xval{i}, R_GMA_unif_uniformed{i}, 'LineWidth', 1);
hold off

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

R_VL_transpose=transpose(R_VL_tot);
non_empty_cells=~cellfun('isempty',R_VL_transpose);
R_VL_TOT=R_VL_transpose(non_empty_cells);

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

R_VM_transpose=transpose(R_VM_tot);
non_empty_cells=~cellfun('isempty',R_VM_transpose);
R_VM_TOT=R_VM_transpose(non_empty_cells);

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

R_G_transpose=transpose(R_G_tot);
non_empty_cells=~cellfun('isempty',R_G_transpose);
R_G_TOT=R_G_transpose(non_empty_cells);



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

R_BFCL_transpose=transpose(R_BFCL_tot);
non_empty_cells=~cellfun('isempty',R_BFCL_transpose);
R_BFCL_TOT=R_BFCL_transpose(non_empty_cells);

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

R_TA_transpose=transpose(R_TA_tot);
non_empty_cells=~cellfun('isempty',R_TA_transpose);
R_TA_TOT=R_TA_transpose(non_empty_cells);

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

R_PL_transpose=transpose(R_PL_tot);
non_empty_cells=~cellfun('isempty',R_PL_transpose);
R_PL_TOT=R_PL_transpose(non_empty_cells);


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

R_SO_transpose=transpose(R_SO_tot);
non_empty_cells=~cellfun('isempty',R_SO_transpose);
R_SO_TOT=R_SO_transpose(non_empty_cells);

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

R_GL_transpose=transpose(R_GL_tot);
non_empty_cells=~cellfun('isempty',R_GL_transpose);
R_GL_TOT=R_GL_transpose(non_empty_cells);


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

R_GMA_transpose=transpose(R_GMA_tot);
non_empty_cells=~cellfun('isempty',R_GMA_transpose);
R_GMA_TOT=R_GMA_transpose(non_empty_cells);


%% MEAN ENVELOPE FOR EACH MUSCLE WITH ROOT MEAN SQUARE


% RF
R_RF_mean_envelope=zeros(size(R_RF_unif_uniformed{1}));
for i=1:numel(R_RF_unif_uniformed)
    R_RF_mean_envelope=R_RF_mean_envelope+R_RF_unif_uniformed{i}.^2;
end
R_RF_mean_envelope=sqrt(R_RF_mean_envelope/numel(R_RF_unif_uniformed));

% RF STD
R_RF_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_RF_values_at_time_i = cellfun(@(x) x(i), R_RF_unif_uniformed);
    R_RF_std(i) = std(R_RF_values_at_time_i);
end

% VL
R_VL_mean_envelope=zeros(size(R_VL_unif_uniformed{1}));
for i=1:numel(R_VL_unif_uniformed)
    R_VL_mean_envelope=R_VL_mean_envelope+R_VL_unif_uniformed{i}.^2;
end
R_VL_mean_envelope=sqrt(R_VL_mean_envelope/numel(R_VL_unif_uniformed));

% VL STD
R_VL_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_VL_values_at_time_i = cellfun(@(x) x(i), R_VL_unif_uniformed);
    R_VL_std(i) = std(R_VL_values_at_time_i);
end

% VM
R_VM_mean_envelope=zeros(size(R_VM_unif_uniformed{1}));
for i=1:numel(R_VM_unif_uniformed)
    R_VM_mean_envelope=R_VM_mean_envelope+R_VM_unif_uniformed{i}.^2;
end
R_VM_mean_envelope=sqrt(R_VM_mean_envelope/numel(R_VM_unif_uniformed));

% VM STD
R_VM_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_VM_values_at_time_i = cellfun(@(x) x(i), R_VM_unif_uniformed);
    R_VM_std(i) = std(R_VM_values_at_time_i);
end

% G
R_G_mean_envelope=zeros(size(R_G_unif_uniformed{1}));
for i=1:numel(R_G_unif_uniformed)
    R_G_mean_envelope=R_G_mean_envelope+R_G_unif_uniformed{i}.^2;
end
R_G_mean_envelope=sqrt(R_G_mean_envelope/numel(R_G_unif_uniformed));

% G STD
R_G_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_G_values_at_time_i = cellfun(@(x) x(i), R_G_unif_uniformed);
    R_G_std(i) = std(R_G_values_at_time_i);
end

% BFCL
R_BFCL_mean_envelope=zeros(size(R_BFCL_unif_uniformed{1}));
for i=1:numel(R_BFCL_unif_uniformed)
    R_BFCL_mean_envelope=R_BFCL_mean_envelope+R_BFCL_unif_uniformed{i}.^2;
end
R_BFCL_mean_envelope=sqrt(R_BFCL_mean_envelope/numel(R_BFCL_unif_uniformed));

% BFCL STD
R_BFCL_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_BFCL_values_at_time_i = cellfun(@(x) x(i), R_BFCL_unif_uniformed);
    R_BFCL_std(i) = std(R_BFCL_values_at_time_i);
end

% TA
R_TA_mean_envelope=zeros(size(R_TA_unif_uniformed{1}));
for i=1:numel(R_TA_unif_uniformed)
    R_TA_mean_envelope=R_TA_mean_envelope+R_TA_unif_uniformed{i}.^2;
end
R_TA_mean_envelope=sqrt(R_TA_mean_envelope/numel(R_TA_unif_uniformed));

% TA STD
R_TA_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_TA_values_at_time_i = cellfun(@(x) x(i), R_TA_unif_uniformed);
    R_TA_std(i) = std(R_TA_values_at_time_i);
end

% PL
R_PL_mean_envelope=zeros(size(R_PL_unif_uniformed{1}));
for i=1:numel(R_PL_unif_uniformed)
    R_PL_mean_envelope=R_PL_mean_envelope+R_PL_unif_uniformed{i}.^2;
end
R_PL_mean_envelope=sqrt(R_PL_mean_envelope/numel(R_PL_unif_uniformed));

% PL STD
R_PL_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_PL_values_at_time_i = cellfun(@(x) x(i), R_PL_unif_uniformed);
    R_PL_std(i) = std(R_PL_values_at_time_i);
end

% SO
R_SO_mean_envelope=zeros(size(R_SO_unif_uniformed{1}));
for i=1:numel(R_SO_unif_uniformed)
    R_SO_mean_envelope=R_SO_mean_envelope+R_SO_unif_uniformed{i}.^2;
end
R_SO_mean_envelope=sqrt(R_SO_mean_envelope/numel(R_SO_unif_uniformed));

% SO STD
R_SO_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_SO_values_at_time_i = cellfun(@(x) x(i), R_SO_unif_uniformed);
    R_SO_std(i) = std(R_SO_values_at_time_i);
end

% GL
R_GL_mean_envelope=zeros(size(R_GL_unif_uniformed{1}));
for i=1:numel(R_GL_unif_uniformed)
    R_GL_mean_envelope=R_GL_mean_envelope+R_GL_unif_uniformed{i}.^2;
end
R_GL_mean_envelope=sqrt(R_GL_mean_envelope/numel(R_GL_unif_uniformed));

% GL STD
R_GL_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_GL_values_at_time_i = cellfun(@(x) x(i), R_GL_unif_uniformed);
    R_GL_std(i) = std(R_GL_values_at_time_i);
end

% GMA
R_GMA_mean_envelope=zeros(size(R_GMA_unif_uniformed{1}));
for i=1:numel(R_GMA_unif_uniformed)
    R_GMA_mean_envelope=R_GMA_mean_envelope+R_GMA_unif_uniformed{i}.^2;
end
R_GMA_mean_envelope=sqrt(R_GMA_mean_envelope/numel(R_GMA_unif_uniformed));

% GMA STD
R_GMA_std = zeros(1, floor(mean_n_samples));
for i = 1:floor(mean_n_samples)
    R_GMA_values_at_time_i = cellfun(@(x) x(i), R_GMA_unif_uniformed);
    R_GMA_std(i) = std(R_GMA_values_at_time_i);
end

mean_envelope_matrix=[R_RF_mean_envelope',R_VL_mean_envelope', R_VM_mean_envelope',R_G_mean_envelope',R_BFCL_mean_envelope',R_TA_mean_envelope',R_PL_mean_envelope',R_SO_mean_envelope',R_GL_mean_envelope',R_GMA_mean_envelope'];

output_folder=uigetdir('','Select destination folder EXCEL'); 
if output_folder==0
    disp('Operation canceled.')%if the folder is empty
else
    input_user=inputdlg('Enter a name for the EXCEL file to save:','Save as');

    if isempty(input_user)
    disp('Operation canceled.')%if the user clicks on 'Annulla'
else
    file_name=input_user{1};
    file_path=fullfile(output_folder,[file_name,'.xlsx']);
    header={'R_RF','R_VL','R_VM','R_G','R_BFCL','R_TA','R_PL','R_SO','R_GL','R_GMA'};
    data_table=cell2table(num2cell(mean_envelope_matrix),'VariableNames', header);
    writetable(data_table, file_path);

    end
end



%% PLOT MEAN ENVELOPE MUSCLES SINGLE WALK

%Vector of percent times:first column of the final matrix
min_time=zeros(length(A),1);
max_time=zeros(length(A),1);
new_time_vector_percentage=cell(length(A),1);
results_time=[];


for i=1:length(A)
min_time(i,1)=min(new_time_vector{i}); 
max_time(i,1)=max(new_time_vector{i});

new_time_vector_percentage{i}=((new_time_vector{i}-min_time(i,1))/(max_time(i,1)-min_time(i,1)))*100;
end

for j=1:numel(new_time_vector_percentage)
    vector_time=new_time_vector_percentage{j}';
    results_time=[results_time;vector_time];
end

%% SUBPLOT MEAN ENVELOPE

figure
subplot(5,2,1)
generalTitle = [ID ' ' Velocity];
sgtitle(generalTitle)
xlim([0 100])
title('EMG Right Rectus Femoris: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_RF_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_RF_mean_envelope + R_RF_std, fliplr(R_RF_mean_envelope - R_RF_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,2)
xlim([0 100])
title('EMG Right Vastus Lateralis: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_VL_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_VL_mean_envelope + R_VL_std, fliplr(R_VL_mean_envelope - R_VL_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,3)
xlim([0 100])
title('EMG Right Vastus Medialis: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_VM_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_VM_mean_envelope + R_VM_std, fliplr(R_VM_mean_envelope - R_VM_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,4)
xlim([0 100])
title('EMG Right Gluteus Medius: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_G_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_G_mean_envelope + R_G_std, fliplr(R_G_mean_envelope - R_G_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,5)
xlim([0 100])
title('EMG Right Biceps femoris caput longus: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_BFCL_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_BFCL_mean_envelope + R_BFCL_std, fliplr(R_BFCL_mean_envelope - R_BFCL_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,6)
xlim([0 100])
title('EMG Tibialis anterior: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_TA_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_TA_mean_envelope + R_TA_std, fliplr(R_TA_mean_envelope - R_TA_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,7)
xlim([0 100])
title('EMG Peroneus longus: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_PL_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_PL_mean_envelope + R_PL_std, fliplr(R_PL_mean_envelope - R_PL_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,8)
xlim([0 100])
title('EMG Soleus: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_SO_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_SO_mean_envelope + R_SO_std, fliplr(R_SO_mean_envelope - R_SO_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,9)
xlim([0 100])
title('EMG Gastrocnemius latralis: gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_GL_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_GL_mean_envelope + R_GL_std, fliplr(R_GL_mean_envelope - R_GL_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,10)
xlim([0 100])
title('EMG Gluteus maximus:gait cycle','Color','b')
hold on
plot(new_time_vector_percentage{1},R_GMA_mean_envelope,'k')
fill([new_time_vector_percentage{1}, fliplr(new_time_vector_percentage{1})],...
    [R_GMA_mean_envelope + R_GMA_std, fliplr(R_GMA_mean_envelope - R_GMA_std)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


%ID column: second column of the final matrix
n_repetitions=((floor(mean_n_samples))*length(A));
results_ID=strings(n_repetitions,1);
for j=1:n_repetitions
    results_ID{j}=ID;
end

%Walk column:third column of the final matrix
results_Velocity=strings(n_repetitions,1);
for j=1:n_repetitions
    results_Velocity{j}=Velocity;
end

if strcmp(Patology,"NORMATIVO")

 values_WALK=repelem(1:length(A),(floor(mean_n_samples)));
 results_WALK=strcat(string(values_WALK(:)),results_Velocity);

elseif ~strcmp(Patology, "NORMATIVO" ) %when patology is ATASSIA or PARAPARESI, Velocity is a blank space
 values_WALK=repelem(1:length(A),(floor(mean_n_samples)));
 results_WALK=strcat(string(values_WALK(:)),results_Velocity);
 
end

%Patology column
results_Patology=strings(n_repetitions,1);
for j=1:n_repetitions
    results_Patology{j}=Patology; 
end

%Phases column

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

values_phases=horzcat(results_heel_contact,results_single_stance,results_terminal_stance,results_swing);
values_phases_transpose=transpose(values_phases);

non_empty_cells = ~cellfun('isempty', values_phases_transpose);%find the non-empty cells

results_phases = values_phases_transpose(non_empty_cells);




%% FINAL MATRIX 

results_matrix=[results_time, results_ID, results_WALK, results_RF,results_VL, results_VM,results_G,results_BFCL,results_TA,results_PL,results_SO,results_GL,results_GMA,results_Patology,results_phases,R_RF_TOT,R_VL_TOT,R_VM_TOT,R_G_TOT,R_BFCL_TOT,R_TA_TOT,R_PL_TOT,R_SO_TOT,R_GL_TOT,R_GMA_TOT];


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










