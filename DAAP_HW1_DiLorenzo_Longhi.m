%% DAAP - Homework A.A. 2024/25
% Di Lorenzo Giuliano
% Longhi Filippo

clc; clear; close all;

addpath('functions')

set(groot,  'DefaultFigureWindowStyle', 'docked', ...       % all figures as tabs in single window
            'DefaultTextInterpreter', 'latex', ...          % interpreter Latex - text and annotations
            'DefaultAxesTickLabelInterpreter', 'latex', ... % interpreter Latex - tick labels
            'DefaultLegendInterpreter', 'latex', ...        % interpreter Latex - legends
            'DefaultLineLineWidth', 1.5, ...                % functions
            'DefaultAxesFontSize', 12, ...                  % axis and title
            'DefaultTextFontSize', 14, ...                  % sgtitle
            'DefaultAxesFontName', 'Times New Roman', ...   % axis and title
            'DefaultTextFontName', 'Times New Roman', ...   % sgtitle
            'DefaultAxesLineWidth', 1, ...                  % axis
            'DefaultConstantLineLineWidth', 1.2, ...        % xline and yline
            'DefaultAxesTitleFontSizeMultiplier', 1.2, ...  % title
            'DefaultFigureColor', 'w', ...                  % background color
            'DefaultAxesBox', 'on', ...                     % plot box
            'DefaultLegendLocation', 'best' ...             % legend position
            );

%% ...
