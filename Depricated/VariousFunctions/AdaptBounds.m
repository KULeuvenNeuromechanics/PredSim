function [bounds] = AdaptBounds(bounds,S,mai)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% adapt bounds based on user input
if ~isempty(S.Bounds)
    bounds.a.lower = ones(size(bounds.a.lower)).*S.Bounds.ActLower;
    %  (1)hip_flex (2)hip_add (3) ip_rot (4)knee (5)ankle (6)sub (7)mtp
    % (8) lumbar ext (9)lumbar bend (10)lumbar rot
    if ~isempty(S.Bounds.ActLowerHip)
        IndexHip = [mai(1).mus.l mai(1).mus.r];
        bounds.a.lower(IndexHip) = S.Bounds.ActLowerHip;
    end
    if ~isempty(S.Bounds.ActLowerKnee)
        IndexKnee = [mai(4).mus.l mai(4).mus.r];
        bounds.a.lower(IndexKnee) = S.Bounds.ActLowerKnee;
    end
    if ~isempty(S.Bounds.ActLowerAnkle)
        IndexAnkle = [mai(5).mus.l mai(5).mus.r];
        bounds.a.lower(IndexAnkle) = S.Bounds.ActLowerAnkle;
    end
end

% adapt bound on final time if this is an input setting
if ~isempty(S.Bounds.tf) & ~isnan(S.Bounds.tf)
    if length(S.Bounds.tf) == 1
        bounds.tf.lower = S.Bounds.tf(1);
        bounds.tf.upper = S.Bounds.tf(1);
    else
        bounds.tf.lower = S.Bounds.tf(1);
        bounds.tf.upper = S.Bounds.tf(2);
    end
end

end

