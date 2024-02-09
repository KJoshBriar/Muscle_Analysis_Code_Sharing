function [r1,r2,r3] = dcm2angle( dcm, varargin );

narginchk(1,3);

[r1, r2, r3] = threeaxisrot(dcm(3,1,:), dcm(3,3,:), -dcm(3,2,:), ...
    dcm(1,2,:), dcm(2,2,:), ...
    -dcm(1,3,:), dcm(1,1,:));

r1 = r1(:);
r2 = r2(:);
r3 = r3(:);

end

function [r1, r2, r3] = threeaxisrot(r11, r12, r21, r31, r32, r11a, r12a)
% find angles for rotations about X, Y, Z axes

    r1 = atan2( r11, r12);
    r2 = asin( r21 );
    r3 = atan2(r31, r32);
    %if strcmpi( lim, 'zeror3')
        for i = find(abs( r21 ) >= 1.0)
            r1(i) = atan2( r11a(i), r12a(i));
            r2(i) = asin( r21(i));
            r3(i) = 0;
        end
    %end
end



