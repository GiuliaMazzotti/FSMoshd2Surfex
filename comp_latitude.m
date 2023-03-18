function [lat] = comp_latitude(ch_y,ch_x)

% Convert Swiss coordinates to WGS84 latitude
% (Rewritten from R function which is published on swisstopo website)

%% Attention: SWISS coordinates have to be integer numbers!

% Convert CH y/x to WGS lat

% Converts militar to civil and  to unit = 1000km
% Axiliary values (% Bern)
y_aux = (ch_y - 600000)./1000000;
x_aux = (ch_x - 200000)./1000000;

% Process lat
lat = 16.9023892+x_aux.*3.238272-(y_aux.^2).*0.270978-(x_aux.^2).*0.002528-(y_aux.^2).*x_aux.*0.0447-(x_aux.^3).*0.0140;

% Unit 10000" to 1 " and converts seconds to degrees (dec)
lat = lat.*100/36;

end
