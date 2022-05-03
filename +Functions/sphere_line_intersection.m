function [Te,inter0] = sphere_line_intersection(c,alt,p_dir,r)
%sphere_line_intersection Calculate the location of intersection of 
% an object at r pointing in the direction of p_dir with a sphere located
% at c at an altitude alt.

Re = 6378.137+alt; %km radius of the earth
rmag = norm(r);
cmag = norm(c);

B = 2*(dot(r,p_dir) - dot(c,p_dir));
C = rmag^2 + cmag^2 - Re^2;

case_check = B^2-4*C;
if case_check >= 0
    P1 = (-B + sqrt(case_check))/2;
    P2 = (-B - sqrt(case_check))/2;
    temp_arr = [P1 P2];
    [~,indx] = min(abs(temp_arr));
    inter0 = abs(temp_arr(indx));
    p = inter0.*p_dir;
    Te = r - p;
else
   Te = 0;
   warning('No intersection');
end


end

