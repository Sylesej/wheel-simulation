function del = delFunk(v)

del = [...
    cos(v) sin(v) 0 0      0      0;...
    sin(v) cos(v) 0 0      0      0;...
    0      0      1 0      0      0;...
    0      0      0 cos(v) sin(v) 0;...
    0      0      0 sin(v) cos(v) 0;...
    0      0      0 0      0      1];