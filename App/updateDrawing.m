function [] = updateDrawing(usr_height,usr_fingertip_elbow,usr_elbow_shoulder,...
    usr_shoulder_width,usr_hip_knee,usr_knee_ground,usr_foot_length,ax1,...
    ink_colour,paper_colour)

%% set defaults for user inputs
default_height = 1.8; % [m]
default_ratio_fingertip_elbow = 1/4;
default_ratio_elbow_shoulder = 1/8;
default_ratio_shoulder_width = 1/4;
default_ratio_hip_knee = 1/4;
default_ratio_knee_ground = 1/4;
default_ratio_foot_length = 1/7;

%% get user input ratios
usr_ratio_fingertip_elbow = usr_fingertip_elbow/usr_height;
usr_ratio_elbow_shoulder = usr_elbow_shoulder/usr_height;
usr_ratio_shoulder_width = usr_shoulder_width/usr_height;
usr_ratio_hip_knee = usr_hip_knee/usr_height;
usr_ratio_knee_ground = usr_knee_ground/usr_height;
usr_ratio_foot_length = usr_foot_length/usr_height;


%% calculate scale factors
sf_foot = usr_ratio_foot_length/default_ratio_foot_length;
sf_low_leg = usr_ratio_knee_ground/default_ratio_knee_ground;
sf_upp_leg = usr_ratio_hip_knee/default_ratio_hip_knee;
sf_upp_arm = usr_ratio_elbow_shoulder/default_ratio_elbow_shoulder;
sf_low_arm = usr_ratio_fingertip_elbow/default_ratio_fingertip_elbow;
sf_shoulder = usr_ratio_shoulder_width/default_ratio_shoulder_width;


%% dimentions for drawing
% foot
dim_foot.length = default_ratio_foot_length*sf_foot;
dim_foot.height = dim_foot.length/4;
dim_foot.width = dim_foot.length/3;
% lower leg
dim_low_leg.length = default_ratio_knee_ground*sf_low_leg - dim_foot.height;
dim_low_leg.width = dim_low_leg.length/3;
% upper leg
dim_upp_leg.length = default_ratio_hip_knee*sf_upp_leg;
dim_upp_leg.width = dim_upp_leg.length/3;
% upper arm
dim_upp_arm.length = default_ratio_elbow_shoulder*sf_upp_arm;
dim_upp_arm.width = dim_upp_arm.length/3;
% hand
dim_hand.length = 1/10;
dim_hand.width = dim_hand.length*0.6;
% lower arm
dim_low_arm.length = (default_ratio_fingertip_elbow*sf_low_arm - dim_hand.length);
dim_low_arm.width = dim_upp_arm.length/3;
% pelvis
dim_pel.width = dim_upp_leg.width*2.2;
dim_pel.height = dim_pel.width*0.4;
% leg length
leg_length = dim_foot.height + dim_low_leg.length + dim_upp_leg.length;
% torso
dim_trs.height = 1-1/6-leg_length;
dim_trs.shwidth = default_ratio_shoulder_width*sf_shoulder;

%% reference points for drawing
pos_hip = [dim_pel.width/2-dim_upp_leg.width/2; leg_length];
pos_knee = [pos_hip(1); pos_hip(2)-dim_upp_leg.length];
pos_ankle = [pos_knee(1); pos_knee(2)-dim_low_leg.length];
pos_shoulder = [dim_trs.shwidth/2; 1-1/6-dim_upp_arm.width/2];
pos_elbow = [pos_shoulder(1); pos_shoulder(2)-dim_upp_arm.length];
pos_wrist = [pos_elbow(1); pos_elbow(2)-dim_low_arm.length];
pos_fingertip = [pos_wrist(1); pos_wrist(2)-dim_hand.length];

%% create figure
% ink_colour = [139,69,19]/256;
% paper_colour = [235,222,173]/256;
% scs = get(0,'ScreenSize');
% f1 = figure('Position',[10,50,scs(3)-20, (scs(4)-140)]);
% set(f1,'Color',paper_colour);
axes(ax1)

% ax1 = gca;
ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';
ax1.Color = paper_colour;
ax1.YLim = [-0.05,1.25];
ax1.DataAspectRatio = [1,1,1]; % axis equal
ax1.NextPlot = 'add'; % hold on

cla(ax1)

%% draw reference cirkle and square
rectangle(ax1,'Position',[-1/2,0,1,1],'EdgeColor',ink_colour)

diam = 1.2;
rectangle(ax1,'Position',[-diam/2,0,diam,diam],'Curvature',1,'EdgeColor',ink_colour)

%% draw man
% head
drawHead(ax1,ink_colour,paper_colour)

% torso
drawTorso(ax1,pos_hip,pos_shoulder,dim_upp_leg.width,dim_upp_arm.width,ink_colour)

% left upper leg
drawLimbSegment(ax1,pos_hip,pos_knee,dim_upp_leg.width,0,ink_colour);
% left lower leg
drawLimbSegment(ax1,pos_knee,pos_ankle,dim_low_leg.width,0,ink_colour);
% left foot
drawFoot(ax1,pos_ankle,dim_low_leg.width,dim_foot.width,dim_foot.length,dim_foot.height,0,ink_colour)

% left upper leg
knee1 = drawLimbSegment(ax1,pos_hip,pos_knee,dim_upp_leg.width,-pi/7,ink_colour,[0.02;0.5]);
% left lower leg
ankle1 = drawLimbSegment(ax1,knee1,pos_ankle-pos_knee+knee1,dim_low_leg.width,-pi/7,ink_colour);
% left foot
drawFoot(ax1,ankle1,dim_low_leg.width,dim_foot.width,dim_foot.length,dim_foot.height,-pi/7,ink_colour)

% left upper arm
elbow1 = drawLimbSegment(ax1,pos_shoulder,pos_elbow,dim_upp_arm.width,-pi/2,ink_colour);
% left lower arm
wrist1 = drawLimbSegment(ax1,elbow1,pos_wrist-pos_elbow+elbow1,dim_low_arm.width,-pi/2,ink_colour);
% left hand
drawLimbSegment(ax1,wrist1,pos_fingertip-pos_wrist+wrist1,dim_hand.width,-pi/2,ink_colour);

% left upper arm
elbow1 = drawLimbSegment(ax1,pos_shoulder,pos_elbow,dim_upp_arm.width,-pi*0.65,ink_colour);
% left lower arm
wrist1 = drawLimbSegment(ax1,elbow1,pos_wrist-pos_elbow+elbow1,dim_low_arm.width,-pi*0.65,ink_colour);
% left hand
drawLimbSegment(ax1,wrist1,pos_fingertip-pos_wrist+wrist1,dim_hand.width,-pi*0.65,ink_colour);

end

%% helper functions

% draw upper/lower arm/leg
function [varargout] = drawLimbSegment(ax,pt1,pt2,width,phi,colr,varargin)
if length(varargin) >= 1
    ptrot = varargin{1};
else
    ptrot = pt1;
end
rot = [cos(phi),sin(phi);-sin(phi),cos(phi)];
axes(ax)

p11 = [pt1(1)-width/2;pt1(2)];
p12 = [pt1(1)-width/2;pt2(2)];
p22 = [pt2(1)+width/2;pt2(2)];
p21 = [pt2(1)+width/2;pt1(2)];
pos = [p11,p12,p22,p21];
pos = pos - ptrot;
for ir=1:size(pos,2)
    tmp = rot*pos(:,ir);
    pos(:,ir) = tmp;
end
pos = pos + ptrot;

r1 = fill(ax,pos(1,:)',pos(2,:)',colr);
r1.FaceAlpha = 0;
r1.LineWidth = 1;
r1.EdgeColor = colr;

r2 = fill(ax,-pos(1,:)',pos(2,:)',colr);
r2.FaceAlpha = 0;
r2.LineWidth = 1;
r2.EdgeColor = colr;

if nargout>=1
    tmp = pt2-ptrot;
    tmp2 = rot*tmp;
    tmp3 = tmp2+ptrot;
    varargout{1} = tmp3;

end
end

% draw head
function [] = drawHead(ax,colr1,colr2)
axes(ax)
r1 = rectangle(ax,'Position',[-1/12,1-1/6,1/6,1/6],'Curvature',0.9);
r1.EdgeColor = colr1;
% r1.FaceColor = colr1;

face_h = 1/8-0.02;
face_w = face_h*0.9;

p1 = [face_w*0.4; 1-face_h];
p2 = [face_w*0.5; 1-1/6];
pos = [p1,p2,p2,p1];
pos(1,3:4) = -pos(1,3:4);
r3 = fill(ax,pos(1,:)',pos(2,:)',colr2);
r3.EdgeColor = colr1;

r2 = rectangle(ax,'Position',[-face_w/2,1-1/8,face_w,face_h],'Curvature',0.8);
r2.EdgeColor = colr1;
r2.FaceColor = colr2;

end

% draw torso
function [] = drawTorso(ax,pt1,pt2,w1,w2,colr)
axes(ax)

p1_l = [abs(pt1(1))+w1/2;pt1(2)]; % pelvis left
p2_l = [p1_l(1); p1_l(2)+(pt2(2)-pt1(2))*0.3];
p3_l = [abs(pt2(1));pt2(2)-w2/2]; % left shoulder
p4_l = [abs(pt2(1));pt2(2)+w2/2];

pos_l = [p1_l,p2_l,p3_l,p4_l];
pos_r = pos_l;
pos_r(1,:) = -pos_r(1,:);
pos = [pos_l,fliplr(pos_r)];

r1 = fill(ax,pos(1,:)',pos(2,:)',colr);
r1.FaceAlpha = 0;
r1.LineWidth = 1;
r1.EdgeColor = colr;

end

% draw foot
function [] = drawFoot(ax,pt1,w0,w,l,h,phi,colr)
axes(ax)
rot = [cos(phi),sin(phi);-sin(phi),cos(phi)];

% side view
xpt = [w0/2,l*0.7,l*0.7,-l*0.3,-w0/2];
ypt = [0,-h*0.7,-h,-h,0];
pos = [xpt;ypt];
for ir=1:size(pos,2)
    tmp = rot*pos(:,ir);
    pos(:,ir) = tmp;
end
pos = pos + pt1;

r1 = fill(ax,pos(1,:)',pos(2,:)',colr);
r1.FaceAlpha = 0;
r1.LineWidth = 1;
r1.EdgeColor = colr;

% front view
xpt = [w0/2,w0/2,-w0/2,-w0/2];
ypt = [0,-h,-h,0];
pos = [xpt;ypt];
for ir=1:size(pos,2)
    tmp = rot'*pos(:,ir);
    pos(:,ir) = tmp;
end
pos = pos + [-pt1(1);pt1(2)];

r2 = fill(ax,pos(1,:)',pos(2,:)',colr);
r2.FaceAlpha = 0;
r2.LineWidth = 1;
r2.EdgeColor = colr;

end
