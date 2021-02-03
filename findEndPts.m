function [P0, t, endPoints, idx_min, idx_max] = findEndPts(ptc)
% This function is to find end points of the line segment from a set of
% point

validateattributes(ptc,{'numeric'},{'real','finite','nonnan','ncols', 3})
P0 = mean(ptc,1);
[t,~,~] = eigenspace(ptc,2);
temp_ptc = ptc-repmat(P0,[size(ptc,1),1]);
P11 = repmat(P0,[size(ptc,1),1]) + sum(repmat(t,[size(ptc,1),1]).*temp_ptc,2)*t; %Project P1 onto the line
s = sum((P11-repmat(P0,[size(ptc,1),1])).*repmat(t,[size(ptc,1),1]),2);
[~,idx_max] = max(s);
[~,idx_min] = min(s);
endPoints = [P11(idx_min,:) P11(idx_max,:)];