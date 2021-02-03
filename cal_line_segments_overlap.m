function Length = cal_line_segments_overlap(LP1,LP2)
% This function is to compute an overlap length between two line section,
% where LP1 is source while LP2 is target
% LP1 = [P1x,P1y,P1z,P2x,P2y,P2z] : source line
% LP2 = [P1x,P1y,P1z,P2x,P2y,P2z] : target line
%Ref: http://stackoverflow.com/questions/22456517/algorithm-for-finding-the-segment-overlapping-two-collinear-segments
%% Demo
% LP1 = [P0_Coord,cur_exam_seg_PP_Coord]
% LP2 = cur_pos_connect_seg_P_Coord
% LP2 = [3,5,0,5,6,0];
%% Project end points of LP2 onto LP1
t_LP1 = (LP1(4:6) - LP1(1:3));%/norm(LP1(4:6) - LP1(1:3));
LP1_P1 = LP1(1:3);
LP2 = reshape(LP2,[],2)';
% Project end points of LP2 onto LP1
P1_LP2 = bsxfun(@minus,LP2,LP1_P1);
LP1_length = norm(t_LP1);
proj_LP2 = sum(bsxfun(@times,P1_LP2,t_LP1),2)./(LP1_length);
s = proj_LP2/LP1_length;
LP2 = bsxfun(@plus,LP1_P1,bsxfun(@times,s,t_LP1));
LP2 = reshape(LP2',[],1)';
% Compute min and max each line segment
% Line 1
min1 = [min(LP1(1),LP1(4)),min(LP1(2),LP1(5)),min(LP1(3),LP1(6))];
max1 = [max(LP1(1),LP1(4)),max(LP1(2),LP1(5)),max(LP1(3),LP1(6))];
% Line 2
min2 = [min(LP2(1),LP2(4)),min(LP2(2),LP2(5)),min(LP2(3),LP2(6))];
max2 = [max(LP2(1),LP2(4)),max(LP2(2),LP2(5)),max(LP2(3),LP2(6))];
% Determine intersection points
minIntersection = [max(min1(1),min2(1)),max(min1(2),min2(2)),max(min1(3),min2(3))];
maxIntersection = [min(max1(1),max2(1)),min(max1(2),max2(2)),min(max1(3),max2(3))];
% Check intersection
intersect = (minIntersection(1)<=maxIntersection(1))&&(minIntersection(2)<=maxIntersection(2))&&(minIntersection(3)<= maxIntersection(3));
% Compute intersection line
if intersect
    Length = norm(maxIntersection - minIntersection);
else
    Length = 0;
end
% Plot results
% figure(10)
% hold all
% plot([LP1(1);LP1(4)],[LP1(2);LP1(5)],'k-','Linewidth',2)
% plot([LP2(1);LP2(4)],[LP2(2);LP2(5)],'r-','Linewidth',3)
% plot(minIntersection(1),minIntersection(2),'r.','Markersize',20)
% plot(maxIntersection(1),maxIntersection(2),'b.','Markersize',20)
