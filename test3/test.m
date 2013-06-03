clear all;
load Exp1_LASSO.mat

load Exp1_Factorized.mat


for i=1:length(Data)
    gNorm(i) = Data{i}.gNorm;
    rNorm(i) = Data{i}.rNorm;
end

deri_clas = gNorm./rNorm

for i =1:size(Data_Fact,1)
    for j = 1:size(Data_Fact,2)
        deri_
