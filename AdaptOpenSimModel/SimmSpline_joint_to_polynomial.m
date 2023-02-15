function [] = SimmSpline_joint_to_polynomial(osim_path)
% --------------------------------------------------------------------------
% SimmSpline_joint_to_polynomial
%   opensimAD is not compatible bith spline-based descriptions of joints.
%   This function replaces them by polynomials
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
% 
% Original author: Lars D'Hondt
% Original date: 25/Nov/2022
%
% Last edit by: Dhruv Gupta
% Last edit date: 15/02/2023
% --------------------------------------------------------------------------

%% load model
import org.opensim.modeling.*;
model = Model(osim_path);

%% loop over joints
for j_j=1:model.getJointSet().getSize()
    joint_j = model.getJointSet().get(j_j-1);
    
    % not all joint classes can have SimmSpline
    if strcmp(joint_j.getConcreteClassName(),"CustomJoint")
        joint_j = CustomJoint.safeDownCast(joint_j);

        % note: might need to add more joint classes
    else
        continue
    end

    % loop over 6 components of spatial transform 
    sptr = joint_j.getSpatialTransform();
    for j_ta=1:6
    
        tr1 = sptr.getTransformAxis(j_ta-1);
        f1 = tr1.getFunction();

        % skip if it is not a SimmSpline
        if strcmp(f1.getConcreteClassName(),"SimmSpline")
            s1 = SimmSpline.safeDownCast(f1);
            sf = 1;
            % if it is a multiplierfunction, check the function that is inside
        elseif strcmp(f1.getConcreteClassName(),"MultiplierFunction")
            mf = MultiplierFunction.safeDownCast(f1);
            if strcmp(mf.getFunction().getConcreteClassName(),"SimmSpline")
                s1 = SimmSpline.safeDownCast(mf.getFunction());
                if mf.getNumProperties>1
                    sf = str2num(string(mf.getPropertyByName('scale')));
                else
                    sf=1;
                end
            else
                continue
            end
        else
            continue
        end


        % get points that define spline
        x1 = s1.getX();
        y1 = s1.getY();
        % ... as a matrix
        x = nan(1,x1.size());
        y = x;
        for i=1:x1.size()
            x(i) = x1.get(i-1);
            y(i) = y1.get(i-1);
        end
        
        y = sf*y;
        
        % get coordinate range
        xrange = [min(x),max(x)];
        
        coord_name = tr1.get_coordinates(0);
        coord = model.getCoordinateSet().get(coord_name);
        xrange(1) = max([xrange(1)-0.1*diff(xrange),coord.get_range(0)]);
        xrange(2) = min([xrange(2)+0.1*diff(xrange),coord.get_range(1)]);
        
        % generate sample points within range
        N = length(find(x>=xrange(1) & x<=xrange(2)))*3;
        xf = linspace(xrange(1),xrange(2),N);

        % evaluate function
        yf = nan(1,N);
        for i=1:N
            x_i = Vector.createFromMat(xf(i));
            yf(i) = double(f1.calcValue(x_i));
        end

%         % include spline points outside range
%         xf = [x(x<xrange(1)),xf,x(x>xrange(2))];
%         yf = [y(x<xrange(1)),yf,y(x>xrange(2))];

        % fit polynomial
        coeffs = polyfit(xf,yf,4);
        
        % plot figure to show quality of fit
        if exist('h1','var') && isa(h1,'matlab.ui.Figure')
            figure(h1);
        else
            h1 = figure;
            tiledlayout('flow');
        end
        nexttile
        p1=xline(coord.get_range(0),'--k','DisplayName','Range');
        hold on
        xline(coord.get_range(1),'--k','DisplayName','Range')
        xs = linspace(min(x),max(x),500);
        ys = polyval(coeffs,xs);
        p2=plot(xs,ys,'DisplayName','Polynomial fit');
        p3=plot(xf,yf,'.','DisplayName','Interpolation points');
        p4=plot(x,y,'*','DisplayName','SimmSpline points');
        title([char(joint_j.getName()) ' / ' char(tr1.getName())],Interpreter="none")
        xlabel(char(coord_name),Interpreter="none")
        ylabel(char(tr1.getName()),Interpreter="none")
        
        % overwrite SimmSpline with polynomial
        coeffs = Vector.createFromMat(coeffs);
        f2 = PolynomialFunction(coeffs);
        tr1.set_function(f2);

    end % end of loop over spatial transform
end % end of loop over joints

% add legend to figure
if exist('h1','var') && isa(h1,'matlab.ui.Figure')
    figure(h1);
    legend([p1,p2,p3,p4],'Location','best')
end

%% save model
model.finalizeConnections();
model.initSystem();
model.print(osim_path);



end