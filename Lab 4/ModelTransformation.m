function [ro] = ModelTransformation(Scale,M,t,ri)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ro = (Scale*M*ri)+t;

end

