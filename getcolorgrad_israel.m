function [gradient_colors_red, gradient_colors_exo, gradient_colors_green, ...
    gradient_colors_eDot, gradient_colors_exp] = getcolorgrad_israel(N_simulations)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% GET COLOR
N = N_simulations; % Number of samples

% FOR LEVELS OF ASSISTANCE
initial_color = '#ffe5e5';  
final_color_red   = '#7f0000';   % red
gradient_colors_red = createColorGradient(initial_color, final_color_red, N);

initial_color = '#ffe5e5';  
final_color_red   = '#D70040';   % red
gradient_colors_exo = createColorGradient(initial_color, final_color_red, N);

initial_color = '#cbd394'; 
final_color_red   = '#6f8817';   % green
gradient_colors_green = createColorGradient(initial_color, final_color_red, N);

initial_color = '#ADD8E6'; 
final_color_blue   = '#00008B';   % Blue
gradient_colors_eDot = createColorGradient(initial_color, final_color_blue, N);

initial_color = '#000000'; 
final_color_black   = '#D3D3D3';   % black
gradient_colors_exp = createColorGradient(initial_color, final_color_black, N);

% FOR MUSCLES
%                ={'soleus_r' 'med_gas_r' 'tib_ant_r' 'bifemsh_r'};
final_color_array={'#0000FF' '#50C878' '#FFBF00' '#7E2F8E'}; %7E2F8E {'#D95319' '#50C878' '#FFBF00' '#7E2F8E'};
alpha= 0.4;
lighter_colors = lighten_colors(final_color_array, alpha);
gradient_colors_mus=double.empty;
for i=1:length(final_color_array)
    initial_color = lighter_colors{i}; % Red ffb2b2
    final_color   = final_color_array{i};   % Blue
    gradient_colors_mus(i,:,:) = createColorGradient(initial_color, final_color, N);
end

% Functions for the colors
function gradient_colors = createColorGradient(initial_color, final_color, N)
    % Convert hexadecimal colors to RGB values
    initial_rgb = hex2rgb(initial_color);
    final_rgb = hex2rgb(final_color);

    % Calculate the step size for each color channel
    step_size = (final_rgb - initial_rgb) / (N - 1);

    % Initialize the gradient_colors array
    gradient_colors = zeros(N, 3);

    % Generate the gradient colors
    for i = 1:N
        gradient_colors(i, :) = initial_rgb + (i - 1) * step_size;
    end

    % Convert the gradient colors back to uint8
    gradient_colors = uint8(gradient_colors * 255);
end

function rgb = hex2rgb(hex)
    % Convert hexadecimal color code to RGB values
    hex = hex(1, 2:end); % Remove '#' from the beginning of the hex code
    rgb = reshape(sscanf(hex, '%2x') / 255, 1, 3);
end

function lighter_colors = lighten_colors(final_color_array, alpha)
    % Initialize the output cell array
    lighter_colors = cell(size(final_color_array));
    
    % Loop through each color in the input array
    for i = 1:numel(final_color_array)
        % Get the current color
        current_color = final_color_array{i};
        
        % Convert the color from hex to RGB
        rgb_color = sscanf(current_color(2:end),'%2x%2x%2x',[1 3])/255;
        
        % Lighten the color based on the alpha value
        lighter_rgb_color = rgb_color + alpha * (1 - rgb_color);
        
        % Convert the RGB color back to hex
        lighter_hex_color = sprintf('#%02X%02X%02X', round(lighter_rgb_color * 255));
        
        % Store the lighter color in the output array
        lighter_colors{i} = lighter_hex_color;
    end
end

end