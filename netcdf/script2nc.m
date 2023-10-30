function script2nc

clear;
close all;

path='/home/javier/Models/Kori_new/hercules_ensemble/';

experiments={'Bassis0N3A-16Enh0.1','Bassis1N3A-16Enh0.1','Bassis0N4A-16Enh0.1','Bassis1N4A-16Enh0.1','Bassis0N3A-19Enh0.1','Bassis1N3A-19Enh0.1','Bassis0N4A-19Enh0.1','Bassis1N4A-19Enh0.1'};
l=strlength(experiments);

for i = 1:l(2)
    disp(strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA.mat','"'));
    mat2nc(strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA.mat','"'),strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA.nc','"'),2); 
    mat2nc_toto(strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA_toto.mat','"'),strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA_toto.nc','"'),2); 
    mat2nc_1d(strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA_toto.mat','"'),strcat('"',path,string(experiments(i)),'/SIA/ASE2km_SIA_1D.nc','"'),2); 
end
