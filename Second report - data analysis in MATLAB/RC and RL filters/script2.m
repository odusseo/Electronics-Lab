%Importing data
%creo delle liste di oggetti che possono essere gruppi di vettori, nel
%nostro caso abbiamo per ogni singolo elemento della lista un gruppo di 3
%vettori: tempi, Vin, Vout. Lo script comprende tre cicli for: uno per
%riempire la lista, uno per creare una matrice Nx45 e l'ultimo per creare
%gli errorbar.
M = 24;
N = 15000;
fine = 19000;
inizio = 1000;
numeri = ones(M,N);
dati = {};
dati_log = {};
incertezze_s = {};
incertezze={};
incertezze_log = zeros(M,N);
dati_selezionati=zeros(M,N);
tempi_0 = [];
tempi_selezionati=zeros(M,N);
dati_selezionati_log=zeros(M,N);
A = [];
B = [];
C = [];
D = [];
E = [];
F = [];
dC = [];
dD = [];
dE = [];
dF = [];
dA=[];
dB=[];
tmpA=[];
tmpB=[];
R_l = [20.003 59.946 109.94 160.31 200.10];
%queste cose serviranno alla fine per fare tutti i grafici


for ii= 0:M
    
filename = strcat('./dati/scope_',int2str(ii),'.csv');
delimiter = ',';
startRow = 30;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
dati{ii} = [dataArray{0:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
dati_log{ii} = [dati{ii}(:,1) log(dati{ii}(:,2)) dati{ii}(:,3)];

end


%il ciclo for serve per prendere tutti i dataArray e posizionarli nella
%nuova lista creata.

for ii=0:24
    
filename = strcat('E:\Esperienza2\dati_', int2str(ii), '.txt');
delimiter = ',';
startRow = 2;
endRow = 2;
formatSpec = '%s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
incertezze_s{ii} = [dataArray{1:end-1}];
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
end

for ii=1:M
        str = char(incertezze_s{ii});
  
    if ii > 40 
        incertezze{ii}= 0.3*0.06*1e-3*str2double(str(12:14));
    else
       incertezze{ii}= 0.3*0.06*str2double(str(12:15));

    end
    
end



% il seguente ciclo crea una matrice 100x45 di numeri casuali compresi tra
% l'inizio e la fine del grafico (parametri scelti precentemente). 
for ii = 1:M

    for jj = 1:N
       
        %bisogna prendere la parte intera perch� i numeri casuali devono
        %essere numeri interi
        numeri(ii,jj) = fix(rand(1)*((fine - inizio)) + inizio);
        
    end
   
    
end

for jj = 1:M
    ii = 1;
while dati{jj}(ii, 1) > 6
   ii=ii+1;
end
tempi_0 = [tempi_0; dati{jj}(ii,1)];
end

for jj=1:M
   
    tempi_selezionati(jj,:)= dati{jj}(numeri(jj,1:end),1)-tempi_0(jj);
    dati_selezionati(jj,:)= dati{jj}(numeri(jj,1:end),2);
    dati_selezionati_log(jj,:)=log(dati{jj}(numeri(jj,1:end),2));
end


clear ii jj numeri str incertezze_s inizio fine dati dati_log;
% questo � il ciclo per creare le figure
for ii = 1:M
 
 immagini{ii} = figure();
 
 %con il comando dati{ii}(numeri(ii,1:end),2) gli dico di prendere il mio
 %ii-esimo elemento della lista, aprirmi la seconda colonna considerando
 %solamente gli elementi della ii-esima colonna della matrice 100x45
 
 incertezze_log(ii,:) = (dati_selezionati(ii,:).^(-1)).*(incertezze{ii}*ones(1,N));
 plots{ii} = errorbar(tempi_selezionati(ii,:), dati_selezionati_log(ii, :),incertezze_log(ii,:), incertezze_log(ii,:), 0*ones(N,1), 0*ones(N,1));
 set(plots{ii}, 'LineStyle', 'none', 'Marker', '.', 'CapSize', 0);
 grid on;
end

for ii=1:M/5
    
    for jj =1:5
        
X = tempi_selezionati((ii-1)*5+jj, :);       
Y = dati_selezionati((ii-1)*5+jj, :);
w_i = 1./(incertezze_log((ii-1)*5+jj, :)).^2;



Delta = sum(w_i)*sum(w_i.*X.^2) - (sum(w_i.*X))^2;
tmpA =[tmpA; (sum(w_i.*X.^2)*sum(w_i.*Y) - sum(w_i.*X)*sum(w_i.*X.*Y))/Delta];
tmpB =[tmpB; (sum(w_i)*sum(w_i.*X.*Y) - sum(w_i.*Y)*sum(w_i.*X))/Delta];

clear X Y w_i;
    end
    
    A = [A; mean(tmpA)];
    B = [B; mean(tmpB)];
    dA = [dA; std(tmpA)];
    dB = [dB; std(tmpB)];
    tmpA = [];
    tmpB = [];
end
close all;
    
    X =  (1./R_l)';    
    Y = -B(1:5);
    w_i = 1./dB(1:5).^(2);



    Delta = sum(w_i)*sum(w_i.*X.^2) - (sum(w_i.*X))^2;
    C = (sum(w_i.*X.^2)*sum(w_i.*Y) - sum(w_i.*X)*sum(w_i.*X.*Y))/Delta;
    D = (sum(w_i)*sum(w_i.*X.*Y) - sum(w_i.*Y)*sum(w_i.*X))/Delta;
    dC = sqrt(sum(w_i.*X.^2)/Delta);
    dD = sqrt(sum(w_i)/Delta);


    clear X Y w_i;
    X = (1./R_scond)';       
    Y = -B(6:9);
    w_i = 1./(dB(6:9)).^2;



    Delta = sum(w_i)*sum(w_i.*X.^2) - (sum(w_i.*X))^2;
    E = (sum(w_i.*X.^2)*sum(w_i.*Y) - sum(w_i.*X)*sum(w_i.*X.*Y))/Delta;
    F = (sum(w_i)*sum(w_i.*X.*Y) - sum(w_i.*Y)*sum(w_i.*X))/Delta;
    dE = sqrt(sum(w_i.*X.^2)/Delta);
    dF = sqrt(sum(w_i)/Delta);

    
   immagini{46} = figure();
   plots{46} = errorbar(1./R_l, -B(1:5), dB(1:5),dB(1:5));
   hold on;
   grid on;
   set(plots{46}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10);
   t = linspace(0,0.15*1e-3,10000);
   y = C + D*t;
   plot(t,y,'r');
   
   immagini{47} = figure();
   plots{47} = errorbar(1./R_scond, -B(6:9), dB(6:9),dB(6:9));
   hold on;
   grid on;
   set(plots{47}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10);
   t = linspace(0,.17*1e-4,10000);
   y = E + F*t;
   plot(t,y,'r');
