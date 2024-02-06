%% DATABASE
folder_path="C:\Users\Utente\Desktop\TESI MAGISTRALE\DATI EMG\DATABESE_MEAN_ENVELOPE_FINAL"; %This folder DATABASE_FINAL only contains the CSV files created using 
% the previous code WALK_TOT.

file_list=dir(fullfile(folder_path,'*.csv'));

table_Database=cell(1,numel(file_list));

for i=1:numel(file_list)
    file_name=fullfile(file_list(i).folder,file_list(i).name);
    table_Database{i}=readtable(file_name);
end

Database=[];

for i=1:numel(table_Database)
    Database=[Database;table_Database{i}];
end

writetable(Database, 'DATABASE_MEAN_ENVELOPE_FINAL_PARAPARESI.csv');
header={'Percentage Gait Cycle','ID','TRIAL','R_RF','R_VL','R_VM','R_G','R_BFCL','R_TA','R_PL','R_SO','R_GL','R_GMA','Patology','Phases','ON/OFF R_RF','ON/OFF R_VL','ON/OFF R_VM','ON/OFF R_G','ON/OFF R_BFCL','ON/OFF R_TA','ON/OFF R_PL','ON/OFF R_SO','ON/OFF R_GL','ON/OFF R_GMA'};
fileID=fopen('DATABASE_MEAN_ENVELOPE_FINAL_PARAPARESI.csv','r');
data=fread(fileID,'*char');
fclose(fileID);
fileID=fopen('DATABASE_MEAN_ENVELOPE_FINAL_PARAPARESI.csv','w');
fprintf(fileID, '%s\n', strjoin(header, ','));
fprintf(fileID, '%s', data);
fclose(fileID);

