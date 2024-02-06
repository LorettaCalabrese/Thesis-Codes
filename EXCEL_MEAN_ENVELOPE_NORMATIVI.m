
clear all
clc
close all

%% NUM_SAMPLES
name_file_excel_NUM_SAMPLES="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\NUM_SAMPLES_NORMATIVI.xlsx";
table_excel_NUM_SAMPLES=readtable(name_file_excel_NUM_SAMPLES);

num_samples=min(table_excel_NUM_SAMPLES);
min_num_samples=table2array(num_samples);


%% PAZIENTI

% ALIFANO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\ALIFANO M.xlsx";
dati_ALIFANO_M=readtable(filename);

% ALIFANO L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\ALIFANO L.xlsx";
dati_ALIFANO_L=readtable(filename);

% MARINI M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\MARINI M.xlsx";
dati_MARINI_M=readtable(filename);

% MARINI L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\MARINI L.xlsx";
dati_MARINI_L=readtable(filename);

% LOSICCO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\LOSICCO M.xlsx";
dati_LOSICCO_M=readtable(filename);

% LOSICCO L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\LOSICCO L.xlsx";
dati_LOSICCO_L=readtable(filename);

% ALESSANDRELLI M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\ALESSANDRELLI M.xlsx";
dati_ALESSANDRELLI_M=readtable(filename);

% ALESSANDRELLI L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\ALESSANDRELLI L.xlsx";
dati_ALESSANDRELLI_L=readtable(filename);

% AMER M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\AMER M.xlsx";
dati_AMER_M=readtable(filename);

% AMER L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\AMER L.xlsx";
dati_AMER_L=readtable(filename);

% CARFORA M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\CARFORA M.xlsx";
dati_CARFORA_M=readtable(filename);

% CARFORA L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\CARFORA L.xlsx";
dati_CARFORA_L=readtable(filename);

% RUBBE M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\RUBBE M.xlsx";
dati_RUBBE_M=readtable(filename);

% RUBBE L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\RUBBE L.xlsx";
dati_RUBBE_L=readtable(filename);

% FERRO M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\FERRO M.xlsx";
dati_FERRO_M=readtable(filename);

% FERRO L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\FERRO L.xlsx";
dati_FERRO_L=readtable(filename);

% EL NAGAR M
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\EL NAGAR M.xlsx";
dati_ELNAGAR_M=readtable(filename);

% EL NAGAR L
filename ="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\NORMATIVI\EXCEL_NORMATIVI\EL NAGAR L.xlsx";
dati_ELNAGAR_L=readtable(filename);

%% FINAL TABLE

new_time_vector=linspace(0,100,floor(min_num_samples));
num_subjects=9;

%% Table M

tabelle_M=cell(num_subjects,1);
tabelle_M{1}=dati_ALIFANO_M;
tabelle_M{2}=dati_MARINI_M;
tabelle_M{3}=dati_LOSICCO_M;
tabelle_M{4}=dati_ALESSANDRELLI_M;
tabelle_M{5}=dati_AMER_M;
tabelle_M{6}=dati_CARFORA_M;
tabelle_M{7}=dati_RUBBE_M;
tabelle_M{8}=dati_FERRO_M;
tabelle_M{9}=dati_ELNAGAR_M;


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
generalTitle = 'HS Medium velocity';
sgtitle(generalTitle)
xlim([0 100])
ylim([0 1])
% title('EMG Right Rectus Femoris: gait cycle','Color','b')
hold on
plot(new_time_vector,R_RF_mean_envelope_M,'k')
fill([new_time_vector,fliplr(new_time_vector)],...
    [R_RF_mean_envelope_M + R_RF_std_M, fliplr(R_RF_mean_envelope_M - R_RF_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off

subplot(5,2,2)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Lateralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VL_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VL_mean_envelope_M + R_VL_std_M, fliplr(R_VL_mean_envelope_M - R_VL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,3)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Medialis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VM_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VM_mean_envelope_M + R_VM_std_M, fliplr(R_VM_mean_envelope_M - R_VM_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,4)
xlim([0 100])
ylim([0 1])
% title('EMG Right Gluteus Medius: gait cycle','Color','b')
hold on
plot(new_time_vector,R_G_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_G_mean_envelope_M + R_G_std_M, fliplr(R_G_mean_envelope_M - R_G_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,5)
xlim([0 100])
ylim([0 1])
% title('EMG Right Biceps femoris caput longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_BFCL_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_BFCL_mean_envelope_M + R_BFCL_std_M, fliplr(R_BFCL_mean_envelope_M - R_BFCL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,6)
xlim([0 100])
ylim([0 1])
% title('EMG Tibialis anterior: gait cycle','Color','b')
hold on
plot(new_time_vector,R_TA_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_TA_mean_envelope_M + R_TA_std_M, fliplr(R_TA_mean_envelope_M - R_TA_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,7)
xlim([0 100])
ylim([0 1])
% title('EMG Peroneus longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_PL_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_PL_mean_envelope_M + R_PL_std_M, fliplr(R_PL_mean_envelope_M - R_PL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,8)
xlim([0 100])
ylim([0 1])
% title('EMG Soleus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_SO_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_SO_mean_envelope_M + R_SO_std_M, fliplr(R_SO_mean_envelope_M - R_SO_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,9)
xlim([0 100])
ylim([0 1])
% title('EMG Gastrocnemius latralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_GL_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GL_mean_envelope_M + R_GL_std_M, fliplr(R_GL_mean_envelope_M - R_GL_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,10)
xlim([0 100])
ylim([0 1])
% title('EMG Gluteus maximus:gait cycle','Color','b')
hold on
plot(new_time_vector,R_GMA_mean_envelope_M,'k')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GMA_mean_envelope_M + R_GMA_std_M, fliplr(R_GMA_mean_envelope_M - R_GMA_std_M)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off




%% Table L

tabelle_L=cell(num_subjects,1);
tabelle_L{1}=dati_ALIFANO_L;
tabelle_L{2}=dati_MARINI_L;
tabelle_L{3}=dati_LOSICCO_L;
tabelle_L{4}=dati_ALESSANDRELLI_L;
tabelle_L{5}=dati_AMER_L;
tabelle_L{6}=dati_CARFORA_L;
tabelle_L{7}=dati_RUBBE_L;
tabelle_L{8}=dati_FERRO_L;
tabelle_L{9}=dati_ELNAGAR_L;


%% INTERPOLATION - L

% RF-L

cell_total_RF_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_RF_L{i}=table2array(tabelle_L{i}(:,1));
end


for i=1:length(cell_total_RF_L)

    x_original = linspace(1, length(cell_total_RF_L{i}), length(cell_total_RF_L{i}));
    x_interp = linspace(1, length(cell_total_RF_L{i}), min_num_samples);
    p_RF_L= polyfit(x_original, cell_total_RF_L{i}, 5);
    interpolated_signal_RF_L = polyval(p_RF_L, x_interp);

interpolated_RF_L{i} = interpolated_signal_RF_L;
end

% VL-L
cell_total_VL_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_VL_L{i}=table2array(tabelle_L{i}(:,2));
end


for i=1:length(cell_total_VL_L)

    x_original = linspace(1, length(cell_total_VL_L{i}), length(cell_total_VL_L{i}));
    x_interp = linspace(1, length(cell_total_VL_L{i}), min_num_samples);
    p_VL_L= polyfit(x_original, cell_total_VL_L{i}, 5);
    interpolated_signal_VL_L = polyval(p_VL_L, x_interp);

interpolated_VL_L{i} = interpolated_signal_VL_L;
end

% VM-L

cell_total_VM_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_VM_L{i}=table2array(tabelle_L{i}(:,3));
end


for i=1:length(cell_total_VM_L)

    x_original = linspace(1, length(cell_total_VM_L{i}), length(cell_total_VM_L{i}));
    x_interp = linspace(1, length(cell_total_VM_L{i}), min_num_samples);
    p_VM_L= polyfit(x_original, cell_total_VM_L{i}, 5);
    interpolated_signal_VM_L = polyval(p_VM_L, x_interp);

interpolated_VM_L{i} = interpolated_signal_VM_L;
end

% G-L

cell_total_G_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_G_L{i}=table2array(tabelle_L{i}(:,4));
end


for i=1:length(cell_total_G_L)

    x_original = linspace(1, length(cell_total_G_L{i}), length(cell_total_G_L{i}));
    x_interp = linspace(1, length(cell_total_G_L{i}), min_num_samples);
    p_G_L= polyfit(x_original, cell_total_G_L{i}, 5);
    interpolated_signal_G_L = polyval(p_G_L, x_interp);

interpolated_G_L{i} = interpolated_signal_G_L;
end

% BFCL-L

cell_total_BFCL_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_BFCL_L{i}=table2array(tabelle_L{i}(:,5));
end


for i=1:length(cell_total_BFCL_L)

    x_original = linspace(1, length(cell_total_BFCL_L{i}), length(cell_total_BFCL_L{i}));
    x_interp = linspace(1, length(cell_total_BFCL_L{i}), min_num_samples);
    p_BFCL_L= polyfit(x_original, cell_total_BFCL_L{i}, 5);
    interpolated_signal_BFCL_L = polyval(p_BFCL_L, x_interp);

interpolated_BFCL_L{i} = interpolated_signal_BFCL_L;
end

% TA-L

cell_total_TA_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_TA_L{i}=table2array(tabelle_L{i}(:,6));
end


for i=1:length(cell_total_TA_L)

    x_original = linspace(1, length(cell_total_TA_L{i}), length(cell_total_TA_L{i}));
    x_interp = linspace(1, length(cell_total_TA_L{i}), min_num_samples);
    p_TA_L= polyfit(x_original, cell_total_TA_L{i}, 5);
    interpolated_signal_TA_L = polyval(p_TA_L, x_interp);

interpolated_TA_L{i} = interpolated_signal_TA_L;
end

% PL-L

cell_total_PL_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_PL_L{i}=table2array(tabelle_L{i}(:,7));
end


for i=1:length(cell_total_PL_L)

    x_original = linspace(1, length(cell_total_PL_L{i}), length(cell_total_PL_L{i}));
    x_interp = linspace(1, length(cell_total_PL_L{i}), min_num_samples);
    p_PL_L= polyfit(x_original, cell_total_PL_L{i}, 5);
    interpolated_signal_PL_L = polyval(p_PL_L, x_interp);

interpolated_PL_L{i} = interpolated_signal_PL_L;
end

% SO-L

cell_total_SO_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_SO_L{i}=table2array(tabelle_L{i}(:,8));
end


for i=1:length(cell_total_SO_L)

    x_original = linspace(1, length(cell_total_SO_L{i}), length(cell_total_SO_L{i}));
    x_interp = linspace(1, length(cell_total_SO_L{i}), min_num_samples);
    p_SO_L= polyfit(x_original, cell_total_SO_L{i}, 5);
    interpolated_signal_SO_L = polyval(p_SO_L, x_interp);

interpolated_SO_L{i} = interpolated_signal_SO_L;
end

% GL-L

cell_total_GL_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_GL_L{i}=table2array(tabelle_L{i}(:,9));
end


for i=1:length(cell_total_GL_L)

    x_original = linspace(1, length(cell_total_GL_L{i}), length(cell_total_GL_L{i}));
    x_interp = linspace(1, length(cell_total_GL_L{i}), min_num_samples);
    p_GL_L= polyfit(x_original, cell_total_GL_L{i}, 5);
    interpolated_signal_GL_L = polyval(p_GL_L, x_interp);

interpolated_GL_L{i} = interpolated_signal_GL_L;
end

% GMA-L

cell_total_GMA_L=cell(num_subjects,1);

for i=1:num_subjects

    cell_total_GMA_L{i}=table2array(tabelle_L{i}(:,10));
end


for i=1:length(cell_total_GMA_L)

    x_original = linspace(1, length(cell_total_GMA_L{i}), length(cell_total_GMA_L{i}));
    x_interp = linspace(1, length(cell_total_GMA_L{i}), min_num_samples);
    p_GMA_L= polyfit(x_original, cell_total_GMA_L{i}, 5);
    interpolated_signal_GMA_L = polyval(p_GMA_L, x_interp);

interpolated_GMA_L{i} = interpolated_signal_GMA_L;
end


%% MEAN ENVELOPE & STD - L

% RF-L
R_RF_mean_envelope_L=zeros(size(interpolated_RF_L{1}));
for i=1:numel(interpolated_RF_L)
    R_RF_mean_envelope_L=R_RF_mean_envelope_L+interpolated_RF_L{i}.^2;
end
R_RF_mean_envelope_L=sqrt(R_RF_mean_envelope_L/numel(interpolated_RF_L));

% RF STD-L
R_RF_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_RF_values_at_time_i_L = cellfun(@(x) x(i), interpolated_RF_L);
    R_RF_std_L(i) = std(R_RF_values_at_time_i_L);
end

% VL-L
R_VL_mean_envelope_L =zeros(size(interpolated_VL_L{1}));
for i=1:numel(interpolated_VL_L)
    R_VL_mean_envelope_L=R_VL_mean_envelope_L+interpolated_VL_L{i}.^2;
end
R_VL_mean_envelope_L=sqrt(R_VL_mean_envelope_L/numel(interpolated_VL_L));

% VL STD-L
R_VL_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_VL_values_at_time_i_L = cellfun(@(x) x(i), interpolated_VL_L);
    R_VL_std_L(i) = std(R_VL_values_at_time_i_L);
end

% VM_L
R_VM_mean_envelope_L=zeros(size(interpolated_VM_L{1}));
for i=1:numel(interpolated_VM_L)
    R_VM_mean_envelope_L=R_VM_mean_envelope_L+interpolated_VM_L{i}.^2;
end
R_VM_mean_envelope_L=sqrt(R_VM_mean_envelope_L/numel(interpolated_VM_L));

% VM STD-L
R_VM_std_M = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_VM_values_at_time_i_L = cellfun(@(x) x(i), interpolated_VM_L);
    R_VM_std_L(i) = std(R_VM_values_at_time_i_L);
end

% G-L
R_G_mean_envelope_L=zeros(size(interpolated_G_L{1}));
for i=1:numel(interpolated_G_L)
    R_G_mean_envelope_L=R_G_mean_envelope_L+interpolated_G_L{i}.^2;
end
R_G_mean_envelope_L=sqrt(R_G_mean_envelope_L/numel(interpolated_G_L));

% G STD-L
R_G_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_G_values_at_time_i_L = cellfun(@(x) x(i), interpolated_G_L);
    R_G_std_L(i) = std(R_G_values_at_time_i_L);
end

% BFCL-L
R_BFCL_mean_envelope_L=zeros(size(interpolated_BFCL_L{1}));
for i=1:numel(interpolated_BFCL_L)
    R_BFCL_mean_envelope_L=R_BFCL_mean_envelope_L+interpolated_BFCL_L{i}.^2;
end
R_BFCL_mean_envelope_L=sqrt(R_BFCL_mean_envelope_L/numel(interpolated_BFCL_L));

% BFCL STD-L
R_BFCL_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_BFCL_values_at_time_i_L = cellfun(@(x) x(i), interpolated_BFCL_L);
    R_BFCL_std_L(i) = std(R_BFCL_values_at_time_i_L);
end

% TA-L
R_TA_mean_envelope_L=zeros(size(interpolated_TA_L{1}));
for i=1:numel(interpolated_TA_L)
    R_TA_mean_envelope_L=R_TA_mean_envelope_L+interpolated_TA_L{i}.^2;
end
R_TA_mean_envelope_L=sqrt(R_TA_mean_envelope_L/numel(interpolated_TA_L));

% TA STD-L
R_TA_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_TA_values_at_time_i_L = cellfun(@(x) x(i), interpolated_TA_L);
    R_TA_std_L(i) = std(R_TA_values_at_time_i_L);
end

% PL-L
R_PL_mean_envelope_L=zeros(size(interpolated_PL_L{1}));
for i=1:numel(interpolated_PL_L)
    R_PL_mean_envelope_L=R_PL_mean_envelope_L+interpolated_PL_L{i}.^2;
end
R_PL_mean_envelope_L=sqrt(R_PL_mean_envelope_L/numel(interpolated_PL_L));

% PL STD-L
R_PL_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_PL_values_at_time_i_L = cellfun(@(x) x(i), interpolated_PL_L);
    R_PL_std_L(i) = std(R_PL_values_at_time_i_L);
end

% SO-L
R_SO_mean_envelope_L=zeros(size(interpolated_SO_L{1}));
for i=1:numel(interpolated_SO_L)
    R_SO_mean_envelope_L=R_SO_mean_envelope_L+interpolated_SO_L{i}.^2;
end
R_SO_mean_envelope_L=sqrt(R_SO_mean_envelope_L/numel(interpolated_SO_L));

% SO STD-L
R_SO_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_SO_values_at_time_i_L = cellfun(@(x) x(i), interpolated_SO_L);
    R_SO_std_L(i) = std(R_SO_values_at_time_i_L);
end

% GL-L
R_GL_mean_envelope_L=zeros(size(interpolated_GL_L{1}));
for i=1:numel(interpolated_GL_L)
    R_GL_mean_envelope_L=R_GL_mean_envelope_L+interpolated_GL_L{i}.^2;
end
R_GL_mean_envelope_L=sqrt(R_GL_mean_envelope_L/numel(interpolated_GL_L));

% GL STD-L
R_GL_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_GL_values_at_time_i_L = cellfun(@(x) x(i), interpolated_GL_L);
    R_GL_std_L(i) = std(R_GL_values_at_time_i_L);
end

% GMA-L
R_GMA_mean_envelope_L=zeros(size(interpolated_GMA_L{1}));
for i=1:numel(interpolated_GMA_L)
    R_GMA_mean_envelope_L=R_GMA_mean_envelope_L+interpolated_GMA_L{i}.^2;
end
R_GMA_mean_envelope_L=sqrt(R_GMA_mean_envelope_L/numel(interpolated_GMA_L));

% GMA STD-L
R_GMA_std_L = zeros(1,floor(min_num_samples));
for i = 1:floor(min_num_samples)
    R_GMA_values_at_time_i_L = cellfun(@(x) x(i), interpolated_GMA_L);
    R_GMA_std_L(i) = std(R_GMA_values_at_time_i_L);
end



%% SUBPLOT-L

figure
subplot(5,2,1)
generalTitle = 'HS Low velocity';
sgtitle(generalTitle)
xlim([0 100])
ylim([0 1])
% title('EMG Right Rectus Femoris: gait cycle','Color','b')
hold on
plot(new_time_vector,R_RF_mean_envelope_L,'b')
fill([new_time_vector,fliplr(new_time_vector)],...
    [R_RF_mean_envelope_L + R_RF_std_L, fliplr(R_RF_mean_envelope_L - R_RF_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off

subplot(5,2,2)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Lateralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VL_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VL_mean_envelope_L + R_VL_std_L, fliplr(R_VL_mean_envelope_L - R_VL_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,3)
xlim([0 100])
ylim([0 1])
% title('EMG Right Vastus Medialis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_VM_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_VM_mean_envelope_L + R_VM_std_L, fliplr(R_VM_mean_envelope_L - R_VM_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,4)
xlim([0 100])
ylim([0 1])
% title('EMG Right Gluteus Medius: gait cycle','Color','b')
hold on
plot(new_time_vector,R_G_mean_envelope_M,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_G_mean_envelope_L + R_G_std_L, fliplr(R_G_mean_envelope_L - R_G_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,5)
xlim([0 100])
ylim([0 1])
% title('EMG Right Biceps femoris caput longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_BFCL_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_BFCL_mean_envelope_L + R_BFCL_std_L, fliplr(R_BFCL_mean_envelope_L - R_BFCL_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,6)
xlim([0 100])
ylim([0 1])
% title('EMG Tibialis anterior: gait cycle','Color','b')
hold on
plot(new_time_vector,R_TA_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_TA_mean_envelope_L + R_TA_std_L, fliplr(R_TA_mean_envelope_L - R_TA_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,7)
xlim([0 100])
ylim([0 1])
% title('EMG Peroneus longus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_PL_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_PL_mean_envelope_L + R_PL_std_L, fliplr(R_PL_mean_envelope_L - R_PL_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,8)
xlim([0 100])
ylim([0 1])
% title('EMG Soleus: gait cycle','Color','b')
hold on
plot(new_time_vector,R_SO_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_SO_mean_envelope_L + R_SO_std_L, fliplr(R_SO_mean_envelope_L - R_SO_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,9)
xlim([0 100])
ylim([0 1])
% title('EMG Gastrocnemius latralis: gait cycle','Color','b')
hold on
plot(new_time_vector,R_GL_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GL_mean_envelope_L + R_GL_std_L, fliplr(R_GL_mean_envelope_L - R_GL_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


subplot(5,2,10)
xlim([0 100])
ylim([0 1])
% title('EMG Gluteus maximus:gait cycle','Color','b')
hold on
plot(new_time_vector,R_GMA_mean_envelope_L,'b')
fill([new_time_vector, fliplr(new_time_vector)],...
    [R_GMA_mean_envelope_L + R_GMA_std_L, fliplr(R_GMA_mean_envelope_L - R_GMA_std_L)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Deviazione Standard');
hold off


