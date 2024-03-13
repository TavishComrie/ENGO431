function [ro_update] = ModelTransformation(Scale,M,t,r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ri = r(:,[2,3,4]);


ro = (Scale*M* transpose(ri))+t;
ro_update = transpose(ro);
end

