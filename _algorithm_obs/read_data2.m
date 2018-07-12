%% read_data2
% Initial script to read the newly generated line scan file 

%% 
nlines = 4;
angles_rad = [0,40,80,120]*pi/180;


data_xls = xlsread('Ca_trial2_originaltranslation.xlsx',...
                   'Compilation of all data');

for I = 1:nlines
    data.distance(:,I) = data_xls(:,3+(I-1)*5);
    data.current(:,I)  = data_xls(:,4+(I-1)*5);
    plot(data.distance(:,I),data.current(:,I)); hold on
end
%% 