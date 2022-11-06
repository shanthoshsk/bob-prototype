function Pts=detectBRISKFeatures(I, varargin)
if isSimMode()
    [Iu8, params] = parseInputs(I,varargin{:});
    PtsStruct=ocvFastHessianDetector(Iu8, params);
    Pts = Points(PtsStruct.Location, PtsStruct);
else
    [I_u8, params] = parseInputs_cg(I,varargin{:});
    
    % get original image size
    nRows = size(I_u8, 1);
    nCols = size(I_u8, 2);
    numInDims = 2;
    
    % column-major (matlab) to row-major (opencv)
    Iu8 = I_u8';
    
    % output variable size and it's size cannot be determined here; 
    % Inside OpenCV algorithm, vector is used to hold output; 
    % Vector is grown by pushing element into it; Once OpenCV computation is
    % done, output size is known, and we use that size to create output
    % memory using malloc; Then elements are copied from OpenCV Vector to EML
    % output buffer
    
    [PtsStruct_Location, PtsStruct_Scale, PtsStruct_Metric, PtsStruct_SignOfLaplacian] = ...
        vision.internal.buildable.fastHessianDetectorBuildable.fastHessianDetector_uint8(Iu8, ...
        int32(nRows), int32(nCols), int32(numInDims), ...
        int32(params.nOctaveLayers), int32(params.nOctaves), int32(params.hessianThreshold));  
    
    PtsStruct.Location        = PtsStruct_Location;
    PtsStruct.Scale           = PtsStruct_Scale;
    PtsStruct.Metric          = PtsStruct_Metric;
    PtsStruct.SignOfLaplacian = PtsStruct_SignOfLaplacian;
    
    Pts = SURFPoints(PtsStruct.Location, PtsStruct);
end

%========================================================================== 
function flag = isSimMode()

flag = isempty(coder.target);

%==========================================================================
% Parse and check inputs - simulation
%==========================================================================
function [Iu8, params] = parseInputs(I, varargin)

validateattributes(I,{'logical', 'uint8', 'int16', 'uint16', ...
    'single', 'double'}, {'2d', 'nonempty', 'nonsparse', 'real'},...
                   'detectSURFFeatures', 'I', 1); %#ok<EMCA>

if isa(I,'uint8')
    Iu8 = I;
else
    Iu8 = im2uint8(I);
end

% Parse the PV pairs
parser = inputParser;
parser.CaseSensitive = true;
parser.addParamValue('MetricThreshold', 1000, @checkMetricThreshold);
parser.addParamValue('NumOctaves',         3, @checkNumOctaves);
parser.addParamValue('NumScaleLevels',     4, @checkNumScaleLevels);

% Parse input
parser.parse(varargin{:});

% Populate the parameters to pass into OpenCV's icvfastHessianDetector()
params.nOctaveLayers    = parser.Results.NumScaleLevels-2;
params.nOctaves         = parser.Results.NumOctaves;
params.hessianThreshold = parser.Results.MetricThreshold;

%==========================================================================
% Parse and check inputs - code-generation
%==========================================================================
function [Iu8, params] = parseInputs_cg(I, varargin)

validateattributes(I,{'logical', 'uint8', 'int16', 'uint16', ...
    'single', 'double'}, {'2d', 'nonempty', 'nonsparse', 'real'},...
                   'detectSURFFeatures', 'I', 1); %#ok<EMCA>

if isa(I,'uint8')
    Iu8 = I;
else
    % im2uint8 does not support codegen
    h_idtc = getImageDataTypeConverter(class(I));
    Iu8 = step(h_idtc,I);
end

% varargin must be non-empty
defaultsVal   = getDefaultParametersVal();
defaultsNoVal = getDefaultParametersNoVal();
properties    = getEmlParserProperties();
optarg = eml_parse_parameter_inputs(defaultsNoVal, properties, varargin{:});
MetricThreshold = (eml_get_parameter_value( ...
        optarg.MetricThreshold, defaultsVal.MetricThreshold, varargin{:}));
NumOctaves = (eml_get_parameter_value( ...
        optarg.NumOctaves, defaultsVal.NumOctaves, varargin{:}));
NumScaleLevels = (eml_get_parameter_value( ...
        optarg.NumScaleLevels, defaultsVal.NumScaleLevels, varargin{:}));        

checkMetricThreshold(MetricThreshold);
checkNumOctaves(NumOctaves);
checkNumScaleLevels(NumScaleLevels);

params.nOctaveLayers    = uint32(NumScaleLevels)-uint32(2);
params.nOctaves         = uint32(NumOctaves);
params.hessianThreshold = uint32(MetricThreshold);

%==========================================================================
function h_idtc = getImageDataTypeConverter(inpClass)

persistent h1 h2 h3 h4 h5

inDTypeIdx = coder.internal.const(getDTypeIdx(inpClass));
% setup method is not required as ImageDataTypeConverter does not have
% tunable parameter

switch inDTypeIdx
      case 1 % double
          if isempty(h1) 
              h1 = vision.ImageDataTypeConverter('OutputDataType','uint8'); 
          end
          h_idtc = h1;   
      case 2 % single
          if isempty(h2) 
              h2 = vision.ImageDataTypeConverter('OutputDataType','uint8'); 
          end
          h_idtc = h2;    
      case 3 % uint16
          if isempty(h3) 
              h3 = vision.ImageDataTypeConverter('OutputDataType','uint8');
          end
          h_idtc = h3;  
      case 4 % int16
          if isempty(h4) 
              h4 = vision.ImageDataTypeConverter('OutputDataType','uint8'); 
          end
          h_idtc = h4;  
      case 5 % logical
          if isempty(h5)
              h5 = vision.ImageDataTypeConverter('OutputDataType','uint8'); 
          end
          h_idtc = h5;          
end

%==========================================================================
function dtIdx = getDTypeIdx(dtClass)

switch dtClass
    case 'double',
        dtIdx = 1;
    case 'single',
        dtIdx = 2;
    case 'uint16',
        dtIdx = 3;
    case 'int16',
        dtIdx = 4;
    case 'logical',
        dtIdx = 5;        
end

%==========================================================================
function defaultsVal = getDefaultParametersVal()

defaultsVal = struct(...
    'MetricThreshold', uint32(1000), ...
    'NumOctaves', uint32(3), ...
    'NumScaleLevels', uint32(4));

%==========================================================================
function defaultsNoVal = getDefaultParametersNoVal()

defaultsNoVal = struct(...
    'MetricThreshold', uint32(0), ... 
    'NumOctaves',      uint32(0), ... 
    'NumScaleLevels',  uint32(0));

%==========================================================================
function properties = getEmlParserProperties()

properties = struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    true, ...
    'PartialMatching', false);

%==========================================================================
function tf = checkMetricThreshold(threshold)
validateattributes(threshold, {'numeric'}, {'scalar','finite',...
    'nonsparse', 'real', 'nonnegative'}, 'detectSURFFeatures'); %#ok<EMCA>
tf = true;

%==========================================================================
function tf = checkNumOctaves(numOctaves)
validateattributes(numOctaves, {'numeric'}, {'integer',... 
    'nonsparse', 'real', 'scalar', 'positive'}, 'detectSURFFeatures'); %#ok<EMCA>
tf = true;

%==========================================================================
function tf = checkNumScaleLevels(scales)
validateattributes(scales, {'numeric'}, {'integer',...
    'nonsparse', 'real', 'scalar', '>=', 3}, 'detectSURFFeatures'); %#ok<EMCA>
tf = true;

