
% clc; clear; close all; warning off

load('2D-Ring.mat');

prompt = {'Grid Accuracy:','Sampling Rate:'};
dlg_title = 'Input Params';
num_lines = 1;
defaultans = {'.05','.5'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

Eps = str2double(answer{1});
eta = str2double(answer{2});
[cncvStat,~] = concvAnals(ring,Eps,eta);

if ~cncvStat
    fprintf('The cluster is convex!\n');
else
    fprintf('The cluster is concave!\n');
end

