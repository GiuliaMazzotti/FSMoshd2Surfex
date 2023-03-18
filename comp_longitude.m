function [lon] = comp_longitude(ch_y,ch_x)

% Convert Swiss coordinates to WGS84 latitude
% (Rewritten from R function which is published on swisstopo website)

%% Attention: SWISS coordinates have to be integer numbers!

% Convert CH y/x to WGS long

% Converts militar to civil and  to unit = 1000km
% Axiliary values (% Bern)
y_aux = (ch_y - 600000)./1000000;
x_aux = (ch_x - 200000)./1000000;

% Process long
lon = 2.6779094 +y_aux.*4.728982+y_aux.*x_aux.*0.791484+y_aux.*(x_aux.^2).*0.1306-(y_aux.^3).*0.0436;

% Unit 10000" to 1 " and converts seconds to degrees (dec)
lon = lon.*100/36;

end