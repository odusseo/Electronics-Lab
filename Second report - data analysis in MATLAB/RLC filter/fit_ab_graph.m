function f = fit_ab_graph(x, y, ddx, ddy, a, b, xunits, yunits, yrunits, cfrx, cfry, cfrresx, cfrresy)
% function f = fit_ab_graph(x, y, ddx, ddy, a, b)
% grafica i dati con il fit lineare y = a + bx trovato da fit_ab
% doppio asse x per grafico dati-modello e residui
% residui con segno y - modello
% cfr da 0 a 5
% da fit_ab, a = fit_ab(2,1) b = fit_ab(2,3)
% oppure, da fit_ab, mettere la stessa matrice in a e b


if numel(ddx) == 1
    
    ddx = ddx * ones(size(x));
    
end

if numel(ddy) == 1
    
    ddy = ddy * ones(size(y));
    
end

if numel(a) == 5 && numel(b) == 5
    
    a = a;
    b = b;
    
end


fig1 = figure();
%grid on;
%hold on;

ax_top = axes('position',[.15 .45 .75 .45]);

h = errorbar(x, y, ddy, ddy, ddx, ddx);
%xlim([5600 7600])

set(h, 'linestyle', 'none');
set(h, 'marker', '.', 'markersize', 10);
set(h, 'markerfacecolor', 'b');
set(h, 'markeredgecolor', 'b');
set(h, 'color', [1 0 0]);
set(h, 'CapSize', 0.7);

hold on;
grid on;

title('Relazione tra resistenza e B','FontSize',13, 'FontName', 'David Libre');

model = b * x + a;

p = plot(x, model);
set(p, 'color', [0.5843 0.8157 0.9882])

yl = ylabel(yunits);
set(yl, 'FontSize', 12);



if (cfrx == 0)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.0f')); 
end
if (cfrx == 1)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.1f')); 
end
if (cfrx == 2)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.2f'));
end
if (cfrx == 3)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.3f')); 
end
if (cfrx == 4)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.4f'));
end
if (cfrx == 5)
    tix=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tix,'%.5f'));
end
if (cfry == 0)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy,'%.0f')); 
end
if (cfry == 1)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy, '%.1f')); 
end
if (cfry == 2)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy,'%.2f')); 
end
if (cfry == 3)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy,'%.3f'));
end
if (cfry == 4)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy,'%.4f')); 
end
if (cfry == 5)
    tiy=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiy,'%.5f'));
end



ax_bot = axes('position',[.15 .15 .75 .23]);

%hold on;

zeta = y - model;
r = errorbar(x, zeta, ddy, ddy, ddx, ddx);
set(r, 'linestyle', 'none')
set(r, 'marker', '.', 'markersize', 3);
set(r, 'markeredgecolor', 'k');
set(r, 'markerfacecolor', 'k');
set(r, 'color', [1 0 0]);
set(r, 'CapSize', 1);

grid on;

xll = xlabel(xunits);
set(xll, 'FontSize', 12);

yll = ylabel(yrunits);
set(yll, 'FontSize', 12);


if (cfrresx == 0)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr,'%.0f')); 
end
if (cfrresx == 1)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr,'%.1f')); 
end
if (cfrresx == 2)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr,'%.2f')); 
end
if (cfrresx == 3)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr,'%.3f')); 
end
if (cfrresx == 4)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr, '%.4f')); 
end
if (cfrresx == 5)
    tixr=get(gca,'xtick')';
    set(gca,'xticklabel',num2str(tixr,'%.5f')); 
end

if (cfrresy == 0)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr,'%.0f')); 
end
if (cfrresy == 1)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr, '%.1f')); 
end
if (cfrresy == 2)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr,'%.2f')); 
end
if (cfrresy == 3)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr,'%.3f')); 
end
if (cfrresy == 4)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr, '%.4f')); 
end
if (cfrresy == 5)
    tiyr=get(gca,'ytick')';
    set(gca,'yticklabel',num2str(tiyr,'%.5f')); 
end

%xlim([0 205])  %fin dove arriva la x(per non avere punti sul bordo)


f = [fig1 ax_top ax_bot];

end

