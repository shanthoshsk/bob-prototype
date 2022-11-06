
classdef BRISKPoints < vision.internal.FeaturePoints
    
   properties (Access='public', Dependent = true)
       %Scale Array of point scales
       Scale;
       %SignOfLaplacian Sign of Laplacian
       SignOfLaplacian;
       %Orientation Orientation of the feature
       Orientation;
   end

   % Internal properties that are accessible only indirectly through
   % dependent properties
   properties (Access='private')
       pScale           = ones(0,1,'single');
       pSignOfLaplacian = ones(0,1,'int8'  );
       pOrientation     = ones(0,1,'single');
   end
   
   methods % Accessors for Dependent properties
       %------------------------------------------------
       function this = set.Scale(this, in)
           this.checkForResizing(in);
           checkScale(in);
           this.pScale = single(in);
       end
       function out = get.Scale(this)
           out = this.pScale;
       end
       %-------------------------------------------------
       function this = set.SignOfLaplacian(this, in)
           this.checkForResizing(in);           
           checkSignOfLaplacian(in);
           this.pSignOfLaplacian = int8(in);
       end
       function out = get.SignOfLaplacian(this)
           out = this.pSignOfLaplacian;
       end
       %-------------------------------------------------
       function this = set.Orientation(this, in)
           this.checkForResizing(in);
           checkOrientation(in);
           this.pOrientation = single(in);
       end
       function out = get.Orientation(this)
           out = this.pOrientation;
       end
   end
   
   %-----------------------------------------------------------------------
   methods(Access=private, Static)
       function name = matlabCodegenRedirect(~)
         name = 'vision.internal.SURFPoints_cg';
       end
   end
   
   %-----------------------------------------------------------------------
   methods (Access='public')
       function this = SURFPoints(varargin)
           %SURFPoints constructor
           this = this@vision.internal.FeaturePoints();
           
           if nargin >= 1
               this = this.append(varargin{:});
           end
       end  
              
       %-------------------------------------------------------------------
       function varargout = plot(this, varargin)
           %plot Plot SURF points
           %
           %   surfPoints.plot plots SURF points in the current axis.
           %
           %   surfPoints.plot(AXES_HANDLE,...) plots using axes with
           %   the handle AXES_HANDLE instead of the current axes (gca).
           %
           %   surfPoints.plot(AXES_HANDLE, PARAM1, VAL1, PARAM2,
           %   VAL2, ...) controls additional plot parameters:
           %
           %      'showScale'   true or false.  When true, a circle 
           %                    proportional to the scale of the detected
           %                    feature is drawn around the point's
           %                    location
           %
           %                    Default: true
           %
           %      'showOrientation' true or false. When true, a line
           %                    corresponding to the point's orientation 
           %                    is drawn from the point's location to the
           %                    edge of the circle indicating the scale
           %
           %                    Default: false
           %
           %   Notes
           %   -----
           %   - Scale of the feature is represented by a circle of
           %     6*Scale radius, which is equivalent to the size of
           %     circular area used by the SURF algorithm to compute 
           %     orientation of the feature
           %
           %   Example
           %   -------
           %   % Extract SURF features
           %   I = imread('cameraman.tif');
           %   points = detectSURFFeatures(I);
           %   [features, valid_points] = extractFeatures(I, points);
           %   % Visualize 10 strongest SURF features, including their 
           %   % scales and orientation which were determined during the 
           %   % descriptor extraction process.
           %   imshow(I); hold on;
           %   strongestPoints = valid_points.selectStrongest(10);
           %   strongestPoints.plot('showOrientation',true);
              
           nargoutchk(0,1);
           
           [h, inputs] = parsePlotInputs(this,varargin{:});
           
           if ~inputs.showScale && ~inputs.showOrientation
               % invoke basic points plot
               plot(h, this.pLocation(:,1), this.pLocation(:,2),'g+');
           else
               phi = linspace(0,2*pi);
               x = cos(phi);  % will cause horizontal line +90 for vertical
               y = sin(phi);
               
               if inputs.showOrientation % Plot orientation
                   % the two zeros result in a horizontal line which
                   % will be rotated at a later stage
                   unitCircle = [x 0; y 0];
               else
                   unitCircle = [x; y];
               end
               
               wasHeld = ishold;
               
               for k = 1:this.Count
                   scale  = 6*this.pScale(k);
                   pt     = this.pLocation(k,:)';
                   % negate the orientation so that it's adjusted for 
                   % plotting in the HG's "image" type which assumes that Y
                   % axis is pointing downward
                   orient = -this.pOrientation(k); 
                   
                   rotationMat = [cos(orient) -sin(orient);...
                       sin(orient) cos(orient)];
                   
                   surfCircle = scale*rotationMat*unitCircle + ...
                       pt*ones(1,size(unitCircle,2));
                   plot(h, pt(1), pt(2), 'g+'); % + marking center
                   % turn the hold state on so that all points are 
                   % plotted; this is necessary since we are plotting in
                   % a loop, otherwise each plot command will overwrite the
                   % previous result
                   if k==1 && ~wasHeld
                       hold('on');
                   end
                   plot(h, surfCircle(1,:),surfCircle(2,:),'g-');
                   
               end
               
               if ~wasHeld
                   hold('off'); % restore original states of hold
               end
           end
          
           if nargout == 1
               varargout{1} = h;
           end
           
       end
       
       %-------------------------------------------------------------------
       function this = selectStrongest(this, N)
       %selectStrongest Return N points with strongest metrics
       %
       %   surfPoints = surfPoints.selectStrongest(N) keeps N points with
       %   strongest metrics.
       %
       %   Example
       %   -------
       %   % create object holding 50 points
       %   points = SURFPoints(ones(50,2), 'Metric', 1:50);
       %   % keep 2 strongest features
       %   points = points.selectStrongest(2)
       
       validateattributes(N, {'numeric'}, {'scalar', 'integer', ...
           'positive'}, mfilename);

       [this.pMetric, idx] = sort(this.pMetric,'descend');
       
       this.pLocation        = this.pLocation(idx,:);
       this.pScale           = this.pScale(idx);
       this.pSignOfLaplacian = this.pSignOfLaplacian(idx);
       this.pOrientation     = this.pOrientation(idx);
       
       this = this.remove(N+1:this.Count);
       
       end
       
   end
              
   methods (Access='public', Hidden=true)
       %-------------------------------------------------------------------
       function this = append(this,varargin)
           %append Appends additional SURF points
           
           indexS = this.Count + 1;
           inputs = parseInputs(this, varargin{:});
           indexE = indexS + size(inputs.Location,1) - 1;
           
           this.pLocation(indexS:indexE, 1:2)     = inputs.Location;
           this.pScale(indexS:indexE, 1)          = inputs.Scale;
           this.pMetric(indexS:indexE, 1)         = inputs.Metric;
           this.pSignOfLaplacian(indexS:indexE,1) = inputs.SignOfLaplacian;
           this.pOrientation(indexS:indexE, 1)    = inputs.Orientation;
       end
   end
   
   methods (Access='protected')
       %-------------------------------------------------------------------
       % Copy data for subsref. This method is used in subsref
       function this = subsref_data(this, option)
           this = subsref_data@vision.internal.FeaturePoints(this, option);
           
           % Scale, SignOfLaplacian, and Orientation are Mx1 matrices. When
           % the indices for sub-referencing is a 1-D array, we explicitly
           % specify the size for the second dimension.
           if length(option.subs) == 1
               option.subs{2} = 1;
           end
           
           this.pScale           = subsref(this.pScale,option);
           this.pSignOfLaplacian = subsref(this.pSignOfLaplacian,option);
           this.pOrientation     = subsref(this.pOrientation,option);
       end       
       
       %-------------------------------------------------------------------
       % Copy data for subsasgn. This method is used in subsasgn
       function this = subsasgn_data(this, option, in)
           this = subsasgn_data@vision.internal.FeaturePoints(this, option, in);
           
           if isempty(in)
               this.pScale = ...
                   subsasgn(this.pScale, option, in);
               this.pSignOfLaplacian = ...
                   subsasgn(this.pSignOfLaplacian, option, in);
               this.pOrientation = ...
                   subsasgn(this.pOrientation, option, in);
           else
               this.pScale = ...
                   subsasgn(this.pScale, option, in.pScale);
               this.pSignOfLaplacian = ...
                   subsasgn(this.pSignOfLaplacian, option, in.pSignOfLaplacian);
               this.pOrientation = ...
                   subsasgn(this.pOrientation, option, in.pOrientation);
           end
       end       

       %------------------------------------------------------------------
       % Removes points; this method is used by selectStrongest
       function this = remove(this, index)
           this.pScale(index,:)           = [];
           this.pSignOfLaplacian(index,:) = [];
           this.pOrientation(index,:)     = [];
           this = remove@vision.internal.FeaturePoints(this, index);
       end
       
   end
   
end

%--------------------------------------------------------------------------
% Plot input parser
%--------------------------------------------------------------------------
function [h, inputs] = parsePlotInputs(~,varargin)

% Parse the PV pairs
parser = inputParser;
parser.CaseSensitive = true;

parser.addOptional('AXIS_HANDLE', [], @vision.internal.validateAxesHandle)

parser.addParamValue('showScale',       true,  @checkFlag);

parser.addParamValue('showOrientation', false, @checkFlag);

% Parse input
parser.parse(varargin{:});

% Assign return values
h = parser.Results.AXIS_HANDLE;

if isempty(h)
    h = gca;
end

inputs.showScale        = parser.Results.showScale;
inputs.showOrientation  = parser.Results.showOrientation;

end

%--------------------------------------------------------------------------
function tf = checkFlag(in)

validateattributes(in, {'logical'}, {'scalar'}, mfilename);

tf = true;
end

%--------------------------------------------------------------------------
% Main parser for the class
%--------------------------------------------------------------------------
function [inputs, unused] = parseInputs(this, varargin)

% Parse the PV pairs
parser = inputParser;
parser.CaseSensitive = true;

parser.addRequired('Location', @(x)checkLocation(this,x));
%parser.addRequired('Location');

parser.addParamValue('Scale',          1.6, @checkScale);
parser.addParamValue('Metric',           0, @(x)checkMetric(this,x));
parser.addParamValue('SignOfLaplacian',  0, @checkSignOfLaplacian);
parser.addParamValue('Orientation',      0, @checkOrientation);

% Parse input
parser.parse(varargin{:});

% Populate the parameters to pass into OpenCV's icvfastHessianDetector()
inputs.Location =        single(parser.Results.Location);
inputs.Scale    =        single(parser.Results.Scale);
inputs.Metric   =        single(parser.Results.Metric);
inputs.SignOfLaplacian = int8(parser.Results.SignOfLaplacian);
inputs.Orientation     = single(parser.Results.Orientation);

numPts = size(inputs.Location,1);

% All parameters must have the same number of elements or be a scalar
validateParamLength(this, numel(inputs.Scale), 'Scale', numPts);
%validateParamLength(numel(inputs.Metric), 'Metric', numPts);
validateParamLength(this, numel(inputs.SignOfLaplacian), 'SignOfLaplacian',...
        numPts);
validateParamLength(this, numel(inputs.Orientation), 'Orientation', numPts);

% Pack the input vectors into column vectors. For convenience, the input
% vector can be Mx1 or 1xN
if ~isscalar(inputs.Scale)
    reshape(inputs.Scale, numPts, 1);
end
% if ~isscalar(inputs.Metric)
%     reshape(inputs.Metric, numPts, 1);
% end
if ~isscalar(inputs.SignOfLaplacian)
    reshape(inputs.SignOfLaplacian, numPts, 1);
end
if ~isscalar(inputs.Orientation)
    reshape(inputs.Orientation, numPts, 1);
end

unused = parser.UsingDefaults;

end

%--------------------------------------------------------------------------
function tf = checkScale(in)
validateattributes(in, {'numeric'}, {'nonnan', ...
    'finite', 'nonsparse', 'real', '>=', 1.6, 'vector'}, ...
    mfilename);

tf = true;
end

%--------------------------------------------------------------------------
function tf = checkSignOfLaplacian(in)
validateattributes(in, {'numeric'}, {'integer', 'real', 'vector'}, ...
    mfilename);

if ~(all(in >= -1) && all( in <= 1))
    error(message('vision:SURFPoints:invalidSignOfLaplacian'));
end

tf = true;
end

%--------------------------------------------------------------------------
function tf = checkOrientation(in)
validateattributes(in, {'numeric'}, {'nonnan', 'finite', 'nonsparse',...
    'real', 'vector'}, mfilename);
tf = true;
end

% LocalWords:  Laplacian
% LocalWords:  OpenCV
