app.GroupName = 'test_group';
app.ModelName = 'test_name';

dist       = '1';
vel        = '2';
llen       = '3';
mass       = '4';
shoulders  = '5';
arm_upper  = '6';
arm_lower  = '7';
leg_upper  = '8';
leg_lower  = '9';
foot       = '10';
strength   = '11';

cmd = ['autoFill.py ' app.GroupName ' ' app.ModelName ...
    ' ' dist ' ' vel ' ' llen ' ' mass ' ' shoulders ' ' arm_upper...
    ' ' arm_lower ' ' leg_upper ' ' leg_lower ' ' foot ' ' strength];

pyrunfile(cmd);