function Handle = MarkerPlot(X, Y, Color, Style, Marker, NOM)

% Num_X = length(X);

% Step = floor(Num_X/NOM);

DelX = (X(end) - X(1))/NOM;

X_M = zeros(1, NOM + 1);
Y_M = zeros(1, NOM + 1);

X_M(1) = X(1);
Y_M(1) = Y(1);

for Index = 2:NOM
    X_c = (Index - 1)*DelX + X(1);
    Diff = X - X_c;
    
    Index_M = find(Diff(1:(end - 1)) <= 0);
    
    X_M(Index) = X(Index_M(end));
    Y_M(Index) = Y(Index_M(end));
end

X_M(end) = X(end);
Y_M(end) = Y(end);

HL = ishold;

plot(X, Y, [Color, Style]);

hold on;
if (strcmp(Marker, '') == 0)
%     Handle = plot(X(1:Step:end), Y(1:Step:end), [Color, Marker]);
    plot(X_M, Y_M, [Color, Marker]);
end

if(HL == 0)
    hold off;
end

Handle = plot(X_M(1), Y_M(1), [Color, Style, Marker]);
return;
end