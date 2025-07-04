% % First off, you need to create an events function telling Matlab to stop integrating when the swing
% % foot has reached the ground, that is, when the second component of the vector px crosses zero
% % 6
% % from a positive to a negative value.

function [position,isterminal,direction] = ground_impact(t,x)
     % create an events function telling Matlab to stop integrating when the swing
     % foot has reached the ground, that is, when the second component of the vector px crosses zero
     % from a positive to a negative value.

     q1 = x(1);
     q2 = x(2);
     q3 = x(3);
     q4 = x(4);
     % q5 = x(5);

     l1 = 0.5;
     l2 = 0.5;
     % px taken from impact map step
     px = [cos(q1 + q2 + q3)/2 + cos(q1 + q2 + q3 + q4)/2 + cos(q1 + q2)/2 + cos(q1)/2;
        sin(q1 + q2 + q3)/2 + sin(q1 + q2 + q3 + q4)/2 + sin(q1 + q2)/2 + sin(q1)/2];
     position = px(2);
     isterminal = 1;
     direction = -1;

 end