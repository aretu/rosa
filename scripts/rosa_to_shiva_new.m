% questo script serve ad adattare i dati di ROSA allo script valido per
% SHIVA

clc
clear all
close all
% cd('C:\Users\Aretusini\Dropbox\Dottorato\ROSA')
cd('C:\LHVRFA\Data\Stefano\')
%% read file
cd csv

[FileName,PathName] = uigetfile('*.CSV');

file = [PathName FileName];

f = fopen(file,'r');  % Open text file
param1=textscan(f, '%*s %f %f %f %f %f %f %f',1,'Delimiter',',','HeaderLines',11);
fclose(f);
f = fopen(file,'r');  % Open text file
param2=textscan(f, '%*s %f %f %f %f %f %f %f',1,'Delimiter',',','HeaderLines',12);
fclose(f);
f = fopen(file,'r');  % Open text file
data=textscan(f, '%f %f %f %f %f %f %f %f', 'Delimiter',',','HeaderLines',14); %read the file considering header the first 14 rows
fclose(f);
f = fopen(file,'r');  % Open text file
freq=textscan(f, '%*s %f',1, 'Delimiter',',','HeaderLines',5);
fclose(f);

cd ../
%% output
% a = size(numeric,1);

Timestamp = repmat(1000/freq{1,1},size(data{1,1},1),1); %crea il vettore dalla freq campionamento (ms)
Axial = data{1,2}; %spinta ass
LVDT = data{1,3}; %lvdt
Torque = data{1,4}; %torque
Angle = data{1,5}; %angle
Encoder2 = data{1,6}; %tachometer
Encoder = data{1,7}; %puls counter
FG = data{1,8}; %function generator

%% supplementary: toggle calibration


h={'Timestamp' 'Torque' 'Axial' 'LVDT' 'Encoder' 'Encoder2' 'FG' 'Angle'};

output=[Timestamp,Torque,Axial,LVDT,Encoder,Encoder2,FG,Angle]; %ho scambiato un con in

output1=output';

%% salvataggio del file


exp=[FileName];
exp=num2str(exp);
ext='.txt';
exp2=strcat(exp,ext);

%% fprintf
% open a file for writing
cd txt
fid = fopen(exp2, 'w');
 
% Table Header

fprintf(fid,'Timestamp\tTorque\tAxial\tLVDT\tEncoder\tEncoder2\tFG\tAngle\n');
 
% print values in column order
fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', output1);
fclose(fid);

cd ../