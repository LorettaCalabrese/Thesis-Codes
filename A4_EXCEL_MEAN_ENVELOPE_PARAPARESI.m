
clear all
clc
close all

%% NUM_SAMPLES
name_file_excel_NUM_SAMPLES="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\NUM_SAMPLES_PARAPARESI.xlsx";
table_excel_NUM_SAMPLES=readtable(name_file_excel_NUM_SAMPLES);

num_samples=min(table_excel_NUM_SAMPLES);
min_num_samples=table2array(num_samples);


%% PAZIENTI

% IVERSA CARLO M
filename = "C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\IVERSA CARLO M.xlsx";
dati_IVERSA_CARLO_M=readtable(filename);

% IVERSA SARA M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\IVERSA SARA M.xlsx";
dati_IVERSA_SARA_M=readtable(filename);

% BENEDETTI FABIANA M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\BENEDETTI FABIANA M.xlsx";
dati_BENEDETTI_FABIANA_M=readtable(filename);

% BENEDETTI ANTONIO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\BENEDETTI ANTONIO M.xlsx";
dati_BENEDETTI_ANTONIO_M=readtable(filename);

% MORETTINI M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\MORETTINI M.xlsx";
dati_MORETTINI_M=readtable(filename);

% MARINI ENNIO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\MARINI ENNIO M.xlsx";
dati_MARINI_ENNIO_M=readtable(filename);

% MARINI ALESSANDRO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\MARINI ALESSANDRO M.xlsx";
dati_MARINI_ALESSANDRO_M=readtable(filename);

% RIZZO FILIPPO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\RIZZO FILIPPO M.xlsx";
dati_RIZZO_FILIPPO_M=readtable(filename);

% RIZZO ENRICA M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\RIZZO ENRICA M.xlsx";
dati_RIZZO_ENRICA_M=readtable(filename);

% D'ANSELMI M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\D'ANSELMI M.xlsx";
dati_DANSELMI_M=readtable(filename);

% TRAIETTI O M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\TRAIETTI O M.xlsx";
dati_TRAIETTI_O_M=readtable(filename);

% TRAIETTI MT
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\TRAIETTI MT M.xlsx";
dati_TRAIETTI_MT_M=readtable(filename);

% MEOLI M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\MEOLI M.xlsx";
dati_MEOLI_M=readtable(filename);

% MELE M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\MELE M.xlsx";
dati_MELE_M=readtable(filename);

% CORMONS M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\PATOLOGICI\PARAPARESI\EXCEL_PARAPARESI\CORMONS M.xlsx";
dati_CORMONS_M=readtable(filename);

%% FINAL TABLE

new_time_vector=linspace(0,100,floor(min_num_samples));
num_subjects=15;

%% Table M

tabelle_M=cell(num_subjects,1);
tabelle_M{1}=dati_IVERSA_CARLO_M;
tabelle_M{2}=dati_IVERSA_SARA_M;
tabelle_M{3}=dati_BENEDETTI_FABIANA_M;
tabelle_M{4}=dati_BENEDETTI_ANTONIO_M;
tabelle_M{5}=dati_MORETTINI_M;
tabelle_M{6}=dati_MARINI_ENNIO_M;
tabelle_M{7}=dati_MARINI_ALESSANDRO_M;
tabelle_M{8}=dati_RIZZO_FILIPPO_M;
tabelle_M{9}=dati_RIZZO_ENRICA_M;
tabelle_M{10}=dati_DANSELMI_M;
tabelle_M{11}=dati_TRAIETTI_O_M;
tabelle_M{12}=dati_TRAIETTI_MT_M;
tabelle_M{13}=dati_MEOLI_M;
tabelle_M{14}=dati_MELE_M;
tabelle_M{15}=dati_CORMONS_M;

%% INTERPOLATION - M

% RF-M

cell_total_RF_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_RF_M{i}=table2array(tabelle_M{i}(:,1));
end


for i=1:length(cell_total_RF_M)

    x_original = linspace(1, length(cell_total_RF_M{i}), length(cell_total_RF_M{i}));
    x_interp = linspace(1, length(cell_total_RF_M{i}), min_num_samples);
    p_RF_M= polyfit(x_original, cell_total_RF_M{i}, 5);
    interpolated_signal_RF_M = polyval(p_RF_M, x_interp);

interpolated_RF_M{i} = interpolated_signal_RF_M;
end

% VL-M
cell_total_VL_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_VL_M{i}=table2array(tabelle_M{i}(:,2));
end


for i=1:length(cell_total_VL_M)

    x_original = linspace(1, length(cell_total_VL_M{i}), length(cell_total_VL_M{i}));
    x_interp = linspace(1, length(cell_total_VL_M{i}), min_num_samples);
    p_VL_M= polyfit(x_original, cell_total_VL_M{i}, 5);
    interpolated_signal_VL_M = polyval(p_VL_M, x_interp);

interpolated_VL_M{i} = interpolated_signal_VL_M;
end

% VM-M

cell_total_VM_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_VM_M{i}=table2array(tabelle_M{i}(:,3));
end


for i=1:length(cell_total_VM_M)

    x_original = linspace(1, length(cell_total_VM_M{i}), length(cell_total_VM_M{i}));
    x_interp = linspace(1, length(cell_total_VM_M{i}), min_num_samples);
    p_VM_M= polyfit(x_original, cell_total_VM_M{i}, 5);
    interpolated_signal_VM_M = polyval(p_VM_M, x_interp);

interpolated_VM_M{i} = interpolated_signal_VM_M;
end

% G-M

cell_total_G_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_G_M{i}=table2array(tabelle_M{i}(:,4));
end


for i=1:length(cell_total_G_M)

    x_original = linspace(1, length(cell_total_G_M{i}), length(cell_total_G_M{i}));
    x_interp = linspace(1, length(cell_total_G_M{i}), min_num_samples);
    p_G_M= polyfit(x_original, cell_total_G_M{i}, 5);
    interpolated_signal_G_M = polyval(p_G_M, x_interp);

interpolated_G_M{i} = interpolated_signal_G_M;
end

% BFCL-M

cell_total_BFCL_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_BFCL_M{i}=table2array(tabelle_M{i}(:,5));
end


for i=1:length(cell_total_BFCL_M)

    x_original = linspace(1, length(cell_total_BFCL_M{i}), length(cell_total_BFCL_M{i}));
    x_interp = linspace(1, length(cell_total_BFCL_M{i}), min_num_samples);
    p_BFCL_M= polyfit(x_original, cell_total_BFCL_M{i}, 5);
    interpolated_signal_BFCL_M = polyval(p_BFCL_M, x_interp);

interpolated_BFCL_M{i} = interpolated_signal_BFCL_M;
end

% TA-M

cell_total_TA_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_TA_M{i}=table2array(tabelle_M{i}(:,6));
end


for i=1:length(cell_total_TA_M)

    x_original = linspace(1, length(cell_total_TA_M{i}), length(cell_total_TA_M{i}));
    x_interp = linspace(1, length(cell_total_TA_M{i}), min_num_samples);
    p_TA_M= polyfit(x_original, cell_total_TA_M{i}, 5);
    interpolated_signal_TA_M = polyval(p_TA_M, x_interp);

interpolated_TA_M{i} = interpolated_signal_TA_M;
end

% PL-M

cell_total_PL_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_PL_M{i}=table2array(tabelle_M{i}(:,7));
end


for i=1:length(cell_total_PL_M)

    x_original = linspace(1, length(cell_total_PL_M{i}), length(cell_total_PL_M{i}));
    x_interp = linspace(1, length(cell_total_PL_M{i}), min_num_samples);
    p_PL_M= polyfit(x_original, cell_total_PL_M{i}, 5);
    interpolated_signal_PL_M = polyval(p_PL_M, x_interp);

interpolated_PL_M{i} = interpolated_signal_PL_M;
end

% SO-M

cell_total_SO_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_SO_M{i}=table2array(tabelle_M{i}(:,8));
end


for i=1:length(cell_total_SO_M)

    x_original = linspace(1, length(cell_total_SO_M{i}), length(cell_total_SO_M{i}));
    x_interp = linspace(1, length(cell_total_SO_M{i}), min_num_samples);
    p_SO_M= polyfit(x_original, cell_total_SO_M{i}, 5);
    interpolated_signal_SO_M = polyval(p_SO_M, x_interp);

interpolated_SO_M{i} = interpolated_signal_SO_M;
end

% GL-M

cell_total_GL_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_GL_M{i}=table2array(tabelle_M{i}(:,9));
end


for i=1:length(cell_total_GL_M)

    x_original = linspace(1, length(cell_total_GL_M{i}), length(cell_total_GL_M{i}));
    x_interp = linspace(1, length(cell_total_GL_M{i}), min_num_samples);
    p_GL_M= polyfit(x_original, cell_total_GL_M{i}, 5);
    interpolated_signal_GL_M = polyval(p_GL_M, x_interp);

interpolated_GL_M{i} = interpolated_signal_GL_M;
end

% GMA-M

cell_total_GMA_M=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_GMA_M{i}=table2array(tabelle_M{i}(:,10));
end


for i=1:length(cell_total_GMA_M)

    x_original = linspace(1, length(cell_total_GMA_M{i}), length(cell_total_GMA_M{i}));
    x_interp = linspace(1, length(cell_total_GMA_M{i}), min_num_samples);
    p_GMA_M= polyfit(x_original, cell_total_GMA_M{i}, 5);
    interpolated_signal_GMA_M = polyval(p_GMA_M, x_interp);

interpolated_GMA_M{i} = interpolated_signal_GMA_M;
end


%% MEAN ENVELOPE & STD -M

% RF-M
R_RF_mean_envelope_M=zeros(size(interpolated_RF_M{1}));
for i=1:numel(interpolated_RF_M)
    R_RF_mean_envelope_M=R_RF_mean_envelope_M+interpolated_RF_M{i}.^2;
end
R_RF_mean_envelope_M=sqrt(R_RF_mean_envelope_M/numel(interpolated_RF_M));

% RF STD-M
R_RF_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_RF_values_at_time_i_M = cellfun(@(x) x(i), interpolated_RF_M);
    R_RF_std_M(i) = std(R_RF_values_at_time_i_M);
end

% VL-M
R_VL_mean_envelope_M=zeros(size(interpolated_VL_M{1}));
for i=1:numel(interpolated_VL_M)
    R_VL_mean_envelope_M=R_VL_mean_envelope_M+interpolated_VL_M{i}.^2;
end
R_VL_mean_envelope_M=sqrt(R_VL_mean_envelope_M/numel(interpolated_VL_M));

% VL STD-M
R_VL_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_VL_values_at_time_i_M = cellfun(@(x) x(i), interpolated_VL_M);
    R_VL_std_M(i) = std(R_VL_values_at_time_i_M);
end

% VM_M
R_VM_mean_envelope_M=zeros(size(interpolated_VM_M{1}));
for i=1:numel(interpolated_VM_M)
    R_VM_mean_envelope_M=R_VM_mean_envelope_M+interpolated_VM_M{i}.^2;
end
R_VM_mean_envelope_M=sqrt(R_VM_mean_envelope_M/numel(interpolated_VM_M));

% VM STD-M
R_VM_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_VM_values_at_time_i_M = cellfun(@(x) x(i), interpolated_VM_M);
    R_VM_std_M(i) = std(R_VM_values_at_time_i_M);
end

% G_M
R_G_mean_envelope_M=zeros(size(interpolated_G_M{1}));
for i=1:numel(interpolated_G_M)
    R_G_mean_envelope_M=R_G_mean_envelope_M+interpolated_G_M{i}.^2;
end
R_G_mean_envelope_M=sqrt(R_G_mean_envelope_M/numel(interpolated_G_M));

% G STD-M
R_G_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_G_values_at_time_i_M = cellfun(@(x) x(i), interpolated_G_M);
    R_G_std_M(i) = std(R_G_values_at_time_i_M);
end

% BFCL-M
R_BFCL_mean_envelope_M=zeros(size(interpolated_BFCL_M{1}));
for i=1:numel(interpolated_BFCL_M)
    R_BFCL_mean_envelope_M=R_BFCL_mean_envelope_M+interpolated_BFCL_M{i}.^2;
end
R_BFCL_mean_envelope_M=sqrt(R_BFCL_mean_envelope_M/numel(interpolated_BFCL_M));

% BFCL STD-M
R_BFCL_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_BFCL_values_at_time_i_M = cellfun(@(x) x(i), interpolated_BFCL_M);
    R_BFCL_std_M(i) = std(R_BFCL_values_at_time_i_M);
end

% TA-M
R_TA_mean_envelope_M=zeros(size(interpolated_TA_M{1}));
for i=1:numel(interpolated_TA_M)
    R_TA_mean_envelope_M=R_TA_mean_envelope_M+interpolated_TA_M{i}.^2;
end
R_TA_mean_envelope_M=sqrt(R_TA_mean_envelope_M/numel(interpolated_TA_M));

% TA STD-M
R_TA_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_TA_values_at_time_i_M = cellfun(@(x) x(i), interpolated_TA_M);
    R_TA_std_M(i) = std(R_TA_values_at_time_i_M);
end

% PL-M
R_PL_mean_envelope_M=zeros(size(interpolated_PL_M{1}));
for i=1:numel(interpolated_PL_M)
    R_PL_mean_envelope_M=R_PL_mean_envelope_M+interpolated_PL_M{i}.^2;
end
R_PL_mean_envelope_M=sqrt(R_PL_mean_envelope_M/numel(interpolated_PL_M));

% PL STD-M
R_PL_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_PL_values_at_time_i_M = cellfun(@(x) x(i), interpolated_PL_M);
    R_PL_std_M(i) = std(R_PL_values_at_time_i_M);
end

% SO-M
R_SO_mean_envelope_M=zeros(size(interpolated_SO_M{1}));
for i=1:numel(interpolated_SO_M)
    R_SO_mean_envelope_M=R_SO_mean_envelope_M+interpolated_SO_M{i}.^2;
end
R_SO_mean_envelope_M=sqrt(R_SO_mean_envelope_M/numel(interpolated_SO_M));

% SO STD-M
R_SO_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_SO_values_at_time_i_M = cellfun(@(x) x(i), interpolated_SO_M);
    R_SO_std_M(i) = std(R_SO_values_at_time_i_M);
end

% GL-M
R_GL_mean_envelope_M=zeros(size(interpolated_GL_M{1}));
for i=1:numel(interpolated_GL_M)
    R_GL_mean_envelope_M=R_GL_mean_envelope_M+interpolated_GL_M{i}.^2;
end
R_GL_mean_envelope_M=sqrt(R_GL_mean_envelope_M/numel(interpolated_GL_M));

% GL STD-M
R_GL_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_GL_values_at_time_i_M = cellfun(@(x) x(i), interpolated_GL_M);
    R_GL_std_M(i) = std(R_GL_values_at_time_i_M);
end

% GMA-M
R_GMA_mean_envelope_M=zeros(size(interpolated_GMA_M{1}));
for i=1:numel(interpolated_GMA_M)
    R_GMA_mean_envelope_M=R_GMA_mean_envelope_M+interpolated_GMA_M{i}.^2;
end
R_GMA_mean_envelope_M=sqrt(R_GMA_mean_envelope_M/numel(interpolated_GMA_M));

% GMA STD-M
R_GMA_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_GMA_values_at_time_i_M = cellfun(@(x) x(i), interpolated_GMA_M);
    R_GMA_std_M(i) = std(R_GMA_values_at_time_i_M);
end



%% SUBPLOT-M

figure
subplot(5,2,1)
generalTitle = 'SPG4 Medium velocity';
sgtitle(generalTitle)
xlim([0 100])
ylim([0 1])
% title('EMG Right Rectus Femoris: gait cycle','Color','b')
hold on
plot(new_time_vector,R_RF_mean_envelope_M,'r')
fill([new_time_vector,fliplr(new_time_vector)],...
    [R_RF_mean_envelope_M + R_RF_std_M, fliplr(R_RF_mean_envelope_M - R_RF_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off

subplot(5,2,2)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Lateralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VL_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VL_mean_envelope_M + R_VL_std_M, fliplr(R_VL_mean_envelope_M - R_VL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,3)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Medialis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VM_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VM_mean_envelope_M + R_VM_std_M, fliplr(R_VM_mean_envelope_M - R_VM_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,4)
xlim([0 100])
ylim([0 1])
% title('EMG Right Gluteus Medius: gait cycle','Color','b')
hold on
plot(new_time_vector,R_G_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_G_mean_envelope_M + R_G_std_M, fliplr(R_G_mean_envelope_M - R_G_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,5)
xlim([0 100])
ylim([0 1])
% title('EMG Right Biceps femoris caput longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_BFCL_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_BFCL_mean_envelope_M + R_BFCL_std_M, fliplr(R_BFCL_mean_envelope_M - R_BFCL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,6)
xlim([0 100])
ylim([0 1])
% title('EMG Tibialis anterior: gait cycle','Color','b')
hold on
plot(new_time_vector,R_TA_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_TA_mean_envelope_M + R_TA_std_M, fliplr(R_TA_mean_envelope_M - R_TA_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,7)
xlim([0 100])
ylim([0 1])
% title('EMG Peroneus longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_PL_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_PL_mean_envelope_M + R_PL_std_M, fliplr(R_PL_mean_envelope_M - R_PL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,8)
xlim([0 100])
ylim([0 1])
% title('EMG Soleus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_SO_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_SO_mean_envelope_M + R_SO_std_M, fliplr(R_SO_mean_envelope_M - R_SO_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,9)
xlim([0 100])
ylim([0 1])
% title('EMG Gastrocnemius latralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_GL_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GL_mean_envelope_M + R_GL_std_M, fliplr(R_GL_mean_envelope_M - R_GL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,10)
xlim([0 100])
ylim([0 1])
% title('EMG Gluteus maximus:gait cycle','Color','b')
hold on
plot(new_time_vector,R_GMA_mean_envelope_M,'r')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GMA_mean_envelope_M + R_GMA_std_M, fliplr(R_GMA_mean_envelope_M - R_GMA_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off

