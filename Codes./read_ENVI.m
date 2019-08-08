function [imageStack,wavelengths] = read_ENVI(fileName,varargin)

% READ_ENVI    [imageStack, wavelengths] = read_ENVI(fileName {,rows,cols,bands})
%
%   Load ENVI spectral image data file.
%
%
%   INPUTS:     'fileName' =  Name of the file without file extension
%                             (e.g., 'test_image')
%
%               'rows'   =  (optional) Vector containing the indices of the
%                                      spectral image rows to be read
%                                      (if empty, all rows will be selected)
%
%               'cols'   =  (optional) Vector containing the indices of the
%                                      spectral image columns to be read
%                                      (if empty, all columns will be selected)
%
%               'bands'  =  (optional) Vector containing the indices of the
%                                      spectral image bands to be read
%                                      (if empty, all bands will be selected)
%
%
%   OUTPUTS:    'imageStack'  =  Spectral image cube
%
%               'wavelengths' =  Wavelengths of the spectral image bands
%
%
%
%   *Example 1*:   Read an ENVI file:
%
%                  read_ENVI('test_image')
%
%
%   *Example 2*:   Read image rows 51-100, all columns, and all bands:
%
%                  read_ENVI('test_image', 51:100)
%
%
%   *Example 3*:   Read image rows 51-100, all columns, and every
%                  second band between bands 1 and 30:
%
%                  read_ENVI('test_image', 51:100, [], 1:2:30)
%
%
% October 19, 2016
% Pauli Fält
%


% --------- ENVI file extensions: ---------

headerExtension = '.hdr';
binaryExtension = '.raw';


% --------- Read data from ENVI header file: ---------

fid = fopen([fileName headerExtension],'r');
headerText = fread(fid,'char');
headerText = char(headerText');
fclose(fid);


% --------- Data interleave types: ---------
% BSQ: Band Sequential (X[col,row,band])
% BIL: Band Interleave by Line (X[col,band,row])
% BIP: Band Interleave by Pixel (X[band,col,row])

interleaveType = regexpi(headerText, '\ninterleave *= *([a-zA-Z]+)\n', 'tokens');
interleaveType = interleaveType{1}{1};

% --------- Data types: ---------
%  1: 1-byte unsigned integer
%  2: 2-byte signed integer
%  3: 4-byte signed integer
%  4: 4-byte float
%  5: 8-byte double
%  9: 2x8-byte complex number made up from 2 doubles
% 12: 2-byte unsigned integer

dataType = regexpi(headerText, '\ndata type *= *([0-9]+)\n', 'tokens');
dataType = str2double(dataType{1});

if dataType == 12
    dataType = 'uint16';
elseif dataType == 8
    dataType = 'uint8';
else
    error('Unknown data type!')
end


% --------- Byte order: ---------
% 0: little-endian byte order
% 1: big-endian byte order

byteOrder = regexpi(headerText, '\nbyte order *= *([0-9]+)\n', 'tokens');
byteOrder = str2double(byteOrder{1});

if byteOrder == 0
    byteOrder = 'ieee-le';
elseif byteOrder == 1
    byteOrder = 'ieee-be';
end


% --------- Number of samples, lines and bands: ---------

samples = regexpi(headerText,'\nsamples *= *([0-9]+)\n', 'tokens');
samples = str2double(samples{1});

lines = regexpi(headerText, '\nlines *= *([0-9]+)\n', 'tokens');
lines = str2double(lines{1});

bands = regexpi(headerText, '\nbands *= *([0-9]+)\n', 'tokens');
bands = str2double(bands{1});


% --------- Wavelengths: ---------

dummyIndex1 = strfind(headerText,'wavelength');
dummyIndex2 = strfind(headerText,'{');
dummyIndex3 = strfind(headerText,'}');

% If there's more than one '{' (there is) select the one following 'wavelength'
if ismatrix(dummyIndex2)
    for I=1:length(dummyIndex2)
        if dummyIndex2(I) >= dummyIndex1
            dummyIndex2 = dummyIndex2(I);
            break;
        end
    end
end

% If there's more than one '}' (there is) select the one following the
% chosen '{'
if ismatrix(dummyIndex3)
    for I=1:length(dummyIndex3)
        if dummyIndex3(I) >= dummyIndex2
            dummyIndex3 = dummyIndex3(I);
        end
    end
end

% Make sure that 'wavelength' is followed by '{', which is followed by '}'
assert(dummyIndex1 < dummyIndex2)
assert(dummyIndex2 < dummyIndex3)

% +/-2 to exclude the braces and the (hopefully) following and preceding
% newline characters
wavelengths = headerText(dummyIndex2+2 : dummyIndex3-2);
wavelengths = textscan(wavelengths,'%f %s');
wavelengths = cell2mat(wavelengths(1));


% --------- Subsetting: ---------

if isempty(varargin)
    rowIndices  = 1:lines;
    colIndices  = 1:samples;
    bandIndices = 1:bands;

elseif length(varargin) == 1
    if isempty(varargin{1})
        rowIndices = 1:lines;
    else
        rowIndices = varargin{1};
    end
    colIndices  = 1:samples;
    bandIndices = 1:bands;

elseif length(varargin) == 2
    if isempty(varargin{1})
        rowIndices = 1:lines;
    else
        rowIndices = varargin{1};
    end
    if isempty(varargin{2})
        colIndices = 1:samples;
    else
        colIndices = varargin{2};
    end
    bandIndices = 1:bands;

elseif length(varargin) == 3
    if isempty(varargin{1})
        rowIndices = 1:lines;
    else
        rowIndices = varargin{1};
    end
    if isempty(varargin{2})
        colIndices = 1:samples;
    else
        colIndices = varargin{2};
    end
    if isempty(varargin{3})
        bandIndices = 1:bands;
    else
        bandIndices = varargin{3};
    end

end

subsetRows  = {'Row','Direct',rowIndices};
subsetCols  = {'Column','Direct',colIndices};
subsetBands = {'Band','Direct',bandIndices};


% --------- Read data from ENVI binary file: ---------

imageStack = multibandread([fileName binaryExtension], ...
    [lines, samples, bands], dataType, 0, interleaveType, byteOrder, ...
    subsetRows, subsetCols, subsetBands);


% --------- Flip the wavelength order from  ---------
% --------- 'long->short' to 'short->long': ---------

if wavelengths(1) > wavelengths(end)
    imageStack  = flip(imageStack,3);
    wavelengths = flip(wavelengths);
end

