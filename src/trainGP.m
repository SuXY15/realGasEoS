% Loading data
dataName = 'c12h26';

dataFile = ['../mech/Alpha/' dataName '.csv'];
paraFile = ['../mech/Alpha/' dataName '_para.csv'];

Data = load(dataFile);
N = size(Data,1);
dim = size(Data,2)-1;
idx = randperm(N);
X = Data(idx,1:dim);
Y = Data(idx,end);

% Training
Ntrain = floor(N*1.0);
gprMdl = fitrgp(X(1:Ntrain,:), Y(1:Ntrain), 'BasisFunction', 'Linear', ...
               'KernelFunction','ardsquaredexponential', 'FitMethod','exact', 'PredictMethod','exact');

Xtest = X(Ntrain+1:end,:);
Ytest = Y(Ntrain+1:end);

Ypred = predict(gprMdl,Xtest);

plot(Y(1:Ntrain), Y(1:Ntrain), 'k--'); hold on;
plot(Y(1:Ntrain), resubPredict(gprMdl), 'rv'); hold on;
plot(Ytest, Ypred, 'g.'); hold on;

disp(gprMdl.Beta)
disp(gprMdl.KernelInformation.KernelParameterNames)
disp(gprMdl.KernelInformation.KernelParameters)

para = [gprMdl.KernelInformation.KernelParameters'; gprMdl.Beta'];

writematrix(para, paraFile);
