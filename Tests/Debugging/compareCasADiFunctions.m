% Compare a function to a reference CasADi-function. 
% If the relative error between each pair of corresponding output elements,
% obtained by evaluating both functions with the same input, is within
% tolerance, the functions are accepted as equal.
%
% mandatory inputs:
%   - f_ref:  handle of the reference function. Should be a CasADi-function
%   - f_test: handle of the function to be tested. Can be CasADi or matlab
%
% additional inputs for settings:
%   - relative_tolerance: (output_ref - output_test)/output_ref < tolerance
%                               default = eps(1) = 2.2204e-16
%   - n_test: repeat the test n times, default = 1
%   - order of magnitude of the inputs to use for testing. Argument values
%           are used as minimum and maximum power of 10, default = [-3,3]
%
%
% Author: Lars D'Hondt
%
% Date: January 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = compareCasADiFunctions(f_ref,f_test,varargin)

import casadi.*

%% settings

% set relative tolerance to accept function outputs as equal
if length(varargin)>=1 && ~isempty(varargin{1})
    relative_tolerance = varargin{1};
else
    relative_tolerance = eps;
end

% repeat the test multiple times
if length(varargin)>=2 && ~isempty(varargin{2})
    n_tests = varargin{2};
else
    n_tests = 1;
end

% order of magnitude of the inputs 10^(...)
if length(varargin)>=3 && ~isempty(varargin{3})
    order_pref = varargin{3};
    if length(order_pref)>=2
        order_min = min(order_pref,[],'all');
        order_max = max(order_pref,[],'all');
        
    else
        order_min = 0;
        order_max = order_pref;
    end
    order_range = abs(order_max - order_min);
else
    order_range = 6;
    order_min = -3;
end

%% prepare test
% get the number of inputs and outputs
if contains(class(f_ref),'casadi')
    % use the reference function if it is a CasADi-function
    n_in = f_ref.n_in;
    n_out = f_ref.n_out;
else
    error('The reference function should be a CasADi-function')
end

% no output is proven equal untill proven otherwise
out_is_eq = zeros(n_tests,n_out);

% prepare cell array to store output differences
out_diff = cell(n_tests,n_out);


%% run test
for j=1:n_tests

    % prepare empty input cell array
    test_input = cell(1,n_in);

    % loop over inputs
    for i=1:n_in
        % get size of the i th input
        n_i1 = f_ref.size1_in(i-1);
        n_i2 = f_ref.size2_in(i-1);
        % generate an input with the proper size, with random elements from
        % range [-1 1]*10^[-3 3]
        test_input{i} = (lhsdesign(n_i1,n_i2)*2-1).*10.^(lhsdesign(n_i1,n_i2)*order_range+order_min);
    end

    % prepare empty output cell array for reference function
    ref_output = cell(1,n_out);
    % prepare empty output cell array for tested function
    test_output = cell(1,n_out);
    
    % evaluate reference function
    [ref_output{:}] = f_ref(test_input{:});
    
    % evaluate tested function
    [test_output{:}] = f_test(test_input{:});
    
    % loop over evaluated outputs
    for i=1:n_out
        % transform output from DM to double
        ref_o = full(ref_output{i});
        test_o = full(test_output{i});
        % absolute difference
        diff_o = ref_o - test_o;
        out_diff{j,i} = diff_o;
        % test relative difference versus tolerance (implicit formulation to
        % avoid- dividing by 0)
        rel_diff = abs(diff_o) - abs(ref_o)*relative_tolerance;
        out_is_eq(j,i) = max(rel_diff < 0);
    end

end

% display outcome
for i=1:n_out
    if min(out_is_eq(:,i))<1
        disp(['Output ' num2str(i) ' is not equal to reference.'])
    end
end
if out_is_eq
    disp('Function outputs are equal to reference function')
end



if nargout>=1
    varargout{1} = out_diff;
end








