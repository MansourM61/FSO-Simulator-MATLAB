function MakeitPretty(FigureHandle, FigureProperties, AxisProperties, PlotProperties, OutputFile)
%MAKEITPRETTY Graphical Modification Function
%
% Make it Pretty function
% Written by Mojtaba Mansour Abadi
%
% Using this MATLAB function, shit gets in, a butterfly goes out!
% If you don't believe, it's OK, no one cares.
%
% To use the function, plot your figure and type:
% MakeitPretty(FH, [FW, FL], [XS, YS], [FS, LW, MS, AS], FN)
%
% FH = Figure handle, you can enter 'gcf' if you are dealing with one figure only.
% FS = desired Text font size.
% LW = desired line width.
% MS = desired marker size.
% AS = desired axes label size
% FW = desired figure width.
% FL = desired figure height.
% XS = desired X axis style: N = linear, G = logarithmic
% YS = desired Y axis style: N = linear, G = logarithmic
% FN = desired file name
%
% Example:
% X = 1:0.1*pi:2*pi;
% hold on;
% plot(X, sin(X), 'b-');
% plot(X, cos(X), 'k-*');
% xlabel('t');
% ylabel('f(t)');
% legend('sin', 'cos');
% MakeitPretty(gcf, [300, 400], 'LNLG', [16, 2.5, 16, 10], 'Figure');
%
% That's it.

Const_TS = 0.75;
Const_TL = [2, 1]*0.01;

Def_FS = 16.0;
Def_AS = 10.0;
Def_XS = 'linear';
Def_YS = 'log';
Def_LW = 1.5;
Def_MS = 5.0;
Def_FW = 10.0;
Def_FL = 9.0;
Def_FN = 'Figure';

switch nargin
    case 0
        FH = gcf;
        FS = Def_FS;
        LW = Def_LW;
        MS = Def_MS;
        AS = Def_AS;
        FW = Def_FW;
        FL = Def_FL;
        XS = Def_XS;
        YS = Def_YS;
        FN = Def_FN;
    case 1
        FH = FigureHandle;
        FS = Def_FS;
        LW = Def_LW;
        MS = Def_MS;
        AS = Def_AS;
        FW = Def_FW;
        FL = Def_FL;
        XS = Def_XS;
        YS = Def_YS;
        FN = Def_FN;
    case 2
        FH = FigureHandle;
        T1 = FigureProperties;
        FW = T1(1);
        FL = T1(2);
        XS = Def_XS;
        YS = Def_YS;
        FS = Def_FS;
        LW = Def_LW;
        MS = Def_MS;
        AS = Def_AS;        
        FN = Def_FN;
    case 3
        FH = FigureHandle;
        T1 = FigureProperties;
        FW = T1(1);
        FL = T1(2);
        T2 = AxisProperties;
        if T2(1) == 'G'
            XS = 'log';
        else
            XS = 'linear';
        end
        if T2(2) == 'G'
            YS = 'log';
        else
            YS = 'linear';
        end        
        FS = Def_FS;
        LW = Def_LW;
        MS = Def_MS;
        AS = Def_AS;        
        FN = Def_FN;
    case 4
        FH = FigureHandle;
        T1 = FigureProperties;
        FW = T1(1);
        FL = T1(2);
        T2 = AxisProperties;
        if T2(1) == 'G'
            XS = 'log';
        else
            XS = 'linear';
        end
        if T2(2) == 'G'
            YS = 'log';
        else
            YS = 'linear';
        end        
        T3 = PlotProperties;
        FS = T3(1);
        LW = T3(2);
        MS = T3(3);
        AS = T3(4);
        FN = Def_FN;
    case 5
        FH = FigureHandle;
        T1 = FigureProperties;
        FW = T1(1);
        FL = T1(2);
        T2 = AxisProperties;
        if T2(1) == 'G'
            XS = 'log';
        else
            XS = 'linear';
        end
        if T2(2) == 'G'
            YS = 'log';
        else
            YS = 'linear';
        end        
        T3 = PlotProperties;
        FS = T3(1);
        LW = T3(2);
        MS = T3(3);
        AS = T3(4);
        FN = OutputFile;
end

set(FH, 'Color', [1, 1, 1])
% set(FH, 'Resize', 'off');
set(FH, 'RendererMode', 'auto');
set(FH, 'Renderer', 'painters');

AX1 = get(FH, 'CurrentAxes');

set(AX1, 'XMinorTick', 'on');
set(AX1, 'YMinorTick', 'on');
set(AX1, 'ZMinorTick', 'on');

set(AX1, 'FontName','Times New Roman');
set(AX1, 'FontSize', FS);
set(AX1, 'Box', 'on');

set(AX1, 'XScale', XS);
set(AX1, 'YScale', YS);

set(AX1, 'LineWidth', Const_TS);
set(AX1, 'TickLength', Const_TL);

% OBJ = get(AX, 'Children');
OBJ = findobj('type', 'line');

% for child = OBJ
%     set(child, 'LineWidth', LW);
%     set(child, 'MarkerSize', MS);
% end

XL = get(AX1, 'XLabel');
set(XL, 'FontName','Times New Roman');
set(XL, 'FontSize', FS);

YL = get(AX1, 'YLabel');
set(YL, 'FontName','Times New Roman');
set(YL, 'FontSize', FS);

ZL = get(AX1, 'ZLabel');
set(ZL, 'FontName','Times New Roman');
set(ZL, 'FontSize', FS);

root = get(FH, 'Parent');

set(root, 'ShowHiddenHandles', 'on');

ANN = findobj('type', 'hggroup');

for annot = ANN
    if(~isprop(annot, 'LineWidth'))
        continue;
    end
    set(annot, 'LineWidth', LW);
end

for annot = ANN
    if(~isprop(annot, 'FontName'))
        continue;
    end
    set(annot, 'FontName','Times New Roman');
    set(annot, 'FontSize', FS);
end

OBJ = findobj('type', 'text');

for child = OBJ
    set(child, 'FontName','Times New Roman');
    set(child, 'FontSize', FS);
end

set(root, 'ShowHiddenHandles', 'off');

set(FH, 'paperunits', 'centimeters', 'paperposition', [0, 0, FW, FL]);

set(AX1, 'FontSize', AS);

print('-dtiff', '-r600', FN);

return;