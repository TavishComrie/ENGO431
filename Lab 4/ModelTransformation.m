function [ro_update] = ModelTransformation(Scale,M,t,r)
    % Summary: the ModelTransformation function determines the object space 
    %          coordinates of a point
    % 
    % Input
    %       Scale:  the scale factor from model to object space
    %       M:      the rotation matrix from model to object space
    %       t:      the translation vector 
    %       r:      the model space coordinates of the point
    % Output
    %       ro_update: the object space coordinates of the point

    % Extracting the model space coordinates 
    ri = r(:,[2,3,4]);

    % Computing the object space coordinates
    ro = (Scale * M * transpose(ri)) + t;

    % Reformatting the object space coordinate vector for output
    ro_update = transpose(ro);
end

