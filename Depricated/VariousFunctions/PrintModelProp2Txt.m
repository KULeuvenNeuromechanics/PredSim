function [] = PrintModelProp2Txt(Modelname,OutFile,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

boolExo = true;
if ~isempty(varargin)
    ExoType = varargin{1};
    boolExo = true;
end

%% Read info from the model
import org.opensim.modeling.*

% make it a bit easier to adapt models in cpp by create tables in advance

m = Model(Modelname);
m.initSystem();
% get the names, mass and inertia of the bodies
nb = m.getBodySet.getSize();
BodyMass = nan(nb,1);
BodyCOM = nan(nb,3);
BodyInertia = nan(nb,3);
BodyNames = cell(nb,1);
Bodies = {'pelvis','femur_l','femur_r','tibia_l','tibia_r','talus_l','talus_r',...
    'calcn_l','calcn_r','toes_l','toes_r','torso','humerus_l','humerus_r','ulna_l',...
    'ulna_r','radius_l','radius_r','hand_l','hand_r'};

for i=1:length(Bodies)
    BSel = m.getBodySet.get(Bodies{i});
    BodyMass(i) = BSel.getMass();
    BodyNames{i} = char(BSel.getName());
    Isel = BSel.getInertia().getMoments();
    for j=1:3
        BodyInertia(i,j) = Isel.get(j-1);
    end
    COMsel = BSel.getMassCenter();
    for j=1:3
        BodyCOM(i,j) = COMsel.get(j-1);
    end
end

JointsModel = {'ground_pelvis','hip_l','hip_r','knee_l','knee_r','ankle_l',...
    'ankle_r','subtalar_l','subtalar_r','mtp_l','mtp_r','back','acromial_l','shoulder_r'...
    'elbow_l','elbow_r','radioulnar_l','radioulnar_r','radius_hand_l','radius_hand_r'};

% get the location in parent and child for the different joints
Joints = m.getJointSet();
nJoints = Joints.getSize();
JointNames = cell(nJoints,1);
LocationInParent = nan(nJoints,3);
LocationInChild = nan(nJoints,3); 
OrientInParent = nan(nJoints,3); 
OrientInChild = nan(nJoints,3); 
ParentNames = cell(nJoints,1);
ChildNames = cell(nJoints,1);
for i=1:length(JointsModel)
    Jsel = Joints.get(JointsModel{i});
    JointNames{i} = char(Jsel.getName());
    if strcmp(JointNames{i},'acromial_l')
        JointNames{i} = 'shoulder_l';
    end
    Child = Jsel.getChildFrame();
    ChildFrame = PhysicalOffsetFrame.safeDownCast(Child);
    Child_transl = ChildFrame.get_translation();
    Child_orient = ChildFrame.get_orientation();
    Parent = Jsel.getParentFrame();
    ParentFrame = PhysicalOffsetFrame.safeDownCast(Parent);
    Parent_transl = ParentFrame.get_translation();
    Parent_orient = ParentFrame.get_orientation();
    for j = 1:3
        LocationInChild(i,j) = Child_transl.get(j-1);
        OrientInChild(i,j) = Child_orient.get(j-1);
        LocationInParent(i,j) = Parent_transl.get(j-1);
        OrientInParent(i,j) = Parent_orient.get(j-1);
    end
    % get names of connected bodies
    ParentNames{i} = char(ParentFrame.getSocket('parent').getConnecteeAsObject.getName());
    ChildNames{i} = char(ChildFrame.getSocket('parent').getConnecteeAsObject.getName());
end
AxisNames = {'st_ground_pelvis','st_hip_l','st_hip_r','st_knee_l','st_knee_r','st_ankle_l',...
    'st_ankle_r','st_subtalar_l','st_subtalar_r','','','st_back','st_sho_l','st_sho_r','st_elb_l',...
    'st_elb_r','st_radioulnar_l','st_radioulnar_r','',''};

% Hard Coded - Joint defenitions
IndCustom = [1:9 12:18];
IndWeld = [19:20];
IndPin = [10:11];
JointDef = cell(20,1);
JointDef(IndCustom) = {'CustomJoint'};
JointDef(IndWeld) = {'WeldJoint'};
JointDef(IndPin) = {'PinJoint'};

%% adapt inertia with exoskeleton
if strcmp(ExoType,'AnkleFoot_Poggensee')
    lTibia = -LocationInParent(6,2);
    
    % exoskeleton properties
    mShank_Exo = 0.88*0.5;
    mFoot_Exo = 0.88*0.5;
    COMShank_Exo = [-0.0381, -lTibia + 0.1016, 0];
    COMFoot_Exo = [0, -0.0381, 0];
    IFootExo = [0.0021, 0.0068, 0.0050];
    IShankExo = [0.0073, 0.0027, 0.0066];
    
    % adapt tibia
    iTibia = find(strcmp(Bodies,'tibia_l') | strcmp(Bodies,'tibia_r'));
    for is = iTibia
        mBone = BodyMass(is);
        COMbone = BodyCOM(is,:);
        IBone = BodyInertia(is,:);
        mExo = mShank_Exo;
        COMexo = COMShank_Exo;
        IExo = IShankExo;
        [NewMass,NewCOM, INew] = AddExoToSegment(mBone,COMbone,IBone,mExo,COMexo,IExo);
        BodyMass(is) = NewMass;
        BodyCOM(is,:) = NewCOM;
        BodyInertia(is,:) = INew;        
    end
    
    % combined COM locations
    iTalus= find(strcmp(Bodies,'talus_l') | strcmp(Bodies,'talus_r'));
    for is = iTalus
        mBone = BodyMass(is);
        COMbone = BodyCOM(is,:);
        IBone = BodyInertia(is,:);
        mExo = mFoot_Exo;
        COMexo = COMFoot_Exo;
        IExo = IFootExo;
        [NewMass,NewCOM, INew] = AddExoToSegment(mBone,COMbone,IBone,mExo,COMexo,IExo);
        BodyMass(is) = NewMass;
        BodyCOM(is,:) = NewCOM;
        BodyInertia(is,:) = INew;        
    end
end

%% Print results to a text file to make the copy - paste a bit easier

fid = fopen(OutFile,'wt');
fprintf( fid, '%s\r\n', 'Export opensim model ');
fprintf( fid, '%s\r\n', ' -- Bodies: -- ');

for i=1:length(BodyNames)
    fprintf( fid, '%s', [BodyNames{i} ' = new OpenSim::Body(']);
    fprintf( fid, '%s', ['"' BodyNames{i} '",' ]);
    
    fprintf( fid, '%.4f%s', BodyMass(i), ', Vec3(');    
    fprintf( fid, '%.4f%s', BodyCOM(i,1), ', ');
    fprintf( fid, '%.4f%s', BodyCOM(i,2), ', ');
    fprintf( fid, '%.4f', BodyCOM(i,3));
    
    fprintf( fid, '%s', '), Inertia(');
    fprintf( fid, '%.4f%s', BodyInertia(i,1), ', ');
    fprintf( fid, '%.4f%s', BodyInertia(i,2), ', ');
    fprintf( fid, '%.4f', BodyInertia(i,3));
    
    fprintf( fid, '%s\n', ', 0, 0, 0));');
end

fprintf( fid, '\n\n');

fprintf( fid, '%s\r\n', ' -- Joints: -- ');

for i=1:length(JointNames)
    fprintf( fid, '%s', [JointNames{i} ' = new ' JointDef{i} ' (']);
    fprintf( fid, '%s', ['"' JointNames{i} '", ' ]);
    if i==1
        fprintf( fid, '%s', 'model->getGround(), Vec3(');    
    else
        fprintf( fid, '%s', ['*' ParentNames{i} ', Vec3(']);    
    end
    fprintf( fid, '%.4f%s', LocationInParent(i,1), ', ');
    fprintf( fid, '%.4f%s', LocationInParent(i,2), ', ');
    fprintf( fid, '%.4f%s', LocationInParent(i,3),'), Vec3(');
    fprintf( fid, '%.4f%s', OrientInParent(i,1), ', ');
    fprintf( fid, '%.4f%s', OrientInParent(i,2), ', ');
    fprintf( fid, '%.4f%s', OrientInParent(i,3),'), ');
    fprintf( fid, '%s', ['*' ChildNames{i} ' , Vec3(']);
    fprintf( fid, '%.4f%s', LocationInChild(i,1), ', ');
    fprintf( fid, '%.4f%s', LocationInChild(i,2), ', ');
    fprintf( fid, '%.4f%s', LocationInChild(i,3),'), Vec3(');
    fprintf( fid, '%.4f%s', OrientInChild(i,1), ', ');
    fprintf( fid, '%.4f%s', OrientInChild(i,2), ', ');
    fprintf( fid, '%.4f%s', OrientInChild(i,3), ')');
    if any(IndCustom == i)
        fprintf( fid, '%s\n',[', ' AxisNames{i} ');']); 
    else
        fprintf( fid, '%s\n',');'); 
    end
end

fclose(fid);

function [NewMass,NewCOM, INew] = AddExoToSegment(mBone,COMBone,IBone,mExo,COMExo,IExo)
    % adapted mass
    NewMass =  mBone + mExo;
    % adapted COM position
    NewCOM = (COMBone*mBone + COMExo*mExo)./(mBone+mExo);
    
    % adapted inertia
    dBone = sqrt(sum((NewCOM(2:3)-COMBone(2:3)).^2));
    dExo = sqrt(sum((NewCOM(2:3)-COMExo(2:3)).^2));
    Ix = IBone(1) + mBone*dBone.^2 + IExo(1) + mExo*dExo.^2;
    
    dBone = sqrt(sum((NewCOM([1 3])-COMBone([1 3])).^2));
    dExo = sqrt(sum((NewCOM([1 3])-COMExo([1 3])).^2));
    Iy = IBone(2) + mBone*dBone.^2 + IExo(2) + mExo*dExo.^2;
    
    dBone = sqrt(sum((NewCOM(1:2)-COMBone(1:2)).^2));
    dExo = sqrt(sum((NewCOM(1:2)-COMExo(1:2)).^2));
    Iz = IBone(3) + mBone*dBone.^2 + IExo(3) + mExo*dExo.^2;
    
    INew = [Ix Iy Iz];
        
end



end

