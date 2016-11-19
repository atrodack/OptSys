classdef OptSysArray  < matlab.mixin.Copyable
    %OptSysArray The fundamental Class from which all others will derive
    % This class is meant to be the basic structure for holding data and
    % doing simple, mathematical calculations. It is the basis for the real
    % space coordinate systems each element will live in.
    %
    % NOTES:
    % "x" is dim 2
    % "y" is dim 1
    %
    % Authors: Alexander Rodack and JLCodona
    % Date: 11/14/2016
    %
    % Special Thanks to JLCodona. This will bear resemblance to the AOGrid
    % class within AOSim2.
    
    
    %% Properties
    
    % Public Properties
    properties(GetAccess='public',SetAccess='public')
        
        % Properties
        name; % name the object internally
        default_array_size = [1,1]; % set a default size for an array to use if no input arg is given to contstructor
        FFTsize = [0,0]; % default size of the fftgrid
        Offset = [0 0]; % moves the object.
        verbose = false; % prints info and plots
        
        % GPU properties
        NGPUs = 0; % default to 0, but is initialized in constructor
        GPU_device; % which GPU is being used
        useGPU = false; % flag to use a detected GPU or use the CPU
                        % set to true if GPU is detected, can be turned off
                        % if desired.
        
        
        % Miscellaneous Properties
        interpolation_method = []; % selects a method of interpolation
        seed = []; % Set seed value for repeatable random numbers
        
    end
    
    
    properties(GetAccess='public',SetAccess='protected')
        
        % Storage arrays
        array_;
        fftarray_ = [];
        
        % Storage arrays on the GPU
        GPUarray_;
        GPUfftarray_ = [];
        
        % Array element spacing in meters
        spacing_ = [0.01,0.01];

        % Default Data type to store in array_
        default_data_type = 'single';

        % Coordinate Properties
        origin_ = [0,0]; % real world coordinates of the "center" of the array
        X_ = []; % cached coordinates for X
        Y_ = [];% cached coordinates for Y
        Xextremes = []; % properties to help COORDS() decide if it needs to
        Yextremes = []; % recompute, or if the cache is still vaild
        Nx_=nan;
        Ny_=nan;
        
        % Timing property
        time_ = [];
    end
    
    %% Methods
    methods
        %% Constructor
        function OSA = OptSysArray(nxy)
            % OSA = OptSysArray(nxy)
            % This is the contructor function for the OptSysArray class. The
            % input argument can be of several types: A scalar value, a
            % two-vector, or another OptSysArray.
            if nargin == 0
                OSA.array_ = zeros(OSA.default_array_size);
            else
                if(isa(nxy,'OptSysArray')) % copy parameters from input
                    OSA.name = ['copy of ' nxy.name];
                    OSA.default_array_size = nxy.default_array_size;
                    OSA.FFTsize = nxy.FFTsize;
                    OSA.verbose = nxy.verbose;
                    OSA.interpolation_method = nxy.interpolation_method;
                    OSA.seed = nxy.seed;
                    
                    OSA.array_ = nxy.array_;
                    OSA.fftarray_ = nxy.fftarray_;
                    OSA.spacing_ = nxy.spacing_;
                    OSA.origin_ = nxy.origin_;
                    OSA.X_ = nxy.X_;
                    OSA.Y_ = nxy.Y_;
                else
                    if(isscalar(nxy)) % make a square array from input
                        OSA.array_ = zeros(nxy,nxy);
                    else % make rectangular array from input
                        OSA.array_ = zeros(nxy(1:2));
                    end
                end
            end
                                        
            % Check GPUs
            OSA.NGPUs = gpuDeviceCount;
            OSA.initGPU;
            
            
            % Set default data type
            OSA.setdatatype();
            
        end % of constructor function
        
                
        %% Read Out Utilities
        function description = describe(OSA)
            % description = describe(OSA)
            % Prints a summary of the important object parameters
            % name, class, [size of array], [coordinate spacing] data tyep,
            % GPU name, GPU device number
            if OSA.NGPUs < 1
                description = sprintf('%s {%s: [%dx%d] [%g,%g] %s}',OSA.name,class(OSA),OSA.nx,OSA.ny,OSA.dx,OSA.dy, OSA.default_data_type);
            else
                description = sprintf('%s {%s: [%dx%d] [%g,%g] %s} {%s : Device %d} ',OSA.name,class(OSA),OSA.nx,OSA.ny,OSA.dx,OSA.dy,OSA.default_data_type,OSA.GPU_device.Name,OSA.GPU_device.Index);
            end
        end % of describe
        
        function value = nx(OSA)
            % value = nx(OSA)
            % Returns number of pixels in "x"
            value = size(OSA.array_,2);
        end % of nx
        
        function value = ny(OSA)
            % value = ny(OSA)
            % Returns number of pixels in "y"
            value = size(OSA.array_,1);
        end % of ny
        
        function value = dx(OSA)
            % value = dx(OSA)
            % Returns the spacing between x samples
            value = OSA.spacing_(2);
        end % of dx
        
        function value = dy(OSA)
            % value = dy(OSA)
            % Returns the spacing between y samples
            value = OSA.spacing_(1);
        end % of ny
        
        function value = dk_(OSA)
            % value = dk_(OSA)
            % Calculates the k spacing
            value = 2*pi./(size(OSA.grid_) .* OSA.spacing_);
        end % of dk

        function sz = size(OSA,dim)  
            % sz = size(OSA,dim)
            % Returns the size of the data array
            % dim is used to specify which dimmension is desired. Providing
            % no input returns both dimmensions
            if nargin < 2
                sz = size(OSA.array);
            else
                sz = size(OSA.array,dim);
            end
        end % of size
        
        
        %% Utilitites for Setting Protected Properties
        
        function o = origin(OSA,varargin)
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    if(isscalar(arg))
                        OSA.origin_ = [1 1]*arg;
                    else
                        OSA.origin_ = [arg(1) arg(2)];
                        OSA.touch;
                    end
                    
                case 2
                    OSA.origin_ = [varargin{1} varargin{2}];
                    OSA.touch;
                    
                otherwise
                    error('too many arguments');
            end
            
            o = OSA.origin_;
        end % of origin
        
        function OSA = setdefaultdatatype(OSA,datatype)
           % OSA = setdefaultdatatype(OSA,datatype)
           % Sets the default_data_type property that is used in the
           % calling of OSA.array(nuarray). datatype must be one of the
           % following strings exactly:
           % 1) 'single'
           % 2) 'double'
           
           OSA.default_data_type = datatype;
           OSA.setdatatype;
        
        end % of setdefaultdatatype
        
        
        
        function a = array(OSA,nuarray,mask) % MAKE GPU FRIENDLY
            % a = array(OSA,nuarray,mask)
            % Multi-functional call
            % 1) Calling with no input arguments returns what is in array_
            % 2) Set the array to a new one, return the new OSA object
            
            if nargin == 1
                a = OSA.array_;
            else % Set the value to the input
                if nargin == 2
                    nuarray = squeeze(nuarray);
                    nuarray = nuarray(:,:,1);
                    nuarray = squeeze(nuarray);
                    if numel(OSA) == numel(nuarray)
                        OSA.array_ = nuarray(:);
                        OSA.setdatatype();
                        OSA.touch;
                    else
                        OSA.resize(size(nuarray));
                        OSA.array_ = nuarray;
                        OSA.setdatatype();
                        OSA.touch;
                    end
                else
                    OSA.array_(mask(:)) = nuarray(:);
                end
                a = OSA;
            end
        end % of array
        
        function OSA = starttiming(OSA)
            % OSA = starttiming(OSA)
            % Starts a stopwatch
            OSA.time_ = tic;
        end % of starttiming
        
        %% Micellaneous Utilities
            
        function OSA = touch(OSA)
            % OSA = touch(OSA)
            % Clears the fft cache
            OSA.fftarray_ = [];
            OSA.GPUfftarray_ = [];
        end % of touch
        
        
        function newsize = resize(OSA,varargin) % MAKE GPU FRIENDLY
            % newsize = OSA.resize(newsize);
            % Change the size of an OptSysArray object's data array.
            % (Assigning a new data array does this automatically.)
            % Note that resize leaves the array spacing unchanged.
            % Use OSA.spacing(dx) to set the new spacing.
            
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};                    
                    if(isscalar(arg))
                        OSA.array_ = zeros([1 1]*arg);
                    else
                        OSA.array_ = zeros(arg(1:2));
                    end
                                       
                case 2
                    OSA.array_ = zeros([varargin{1} varargin{2}]);
                    
                otherwise
                    error('too many arguments');
            end
            
            newsize = size(OSA);
            OSA.touch;
        end % of resize
            
        
        function OSA = setdatatype(OSA,default_data_type) % MAKE GPU FRIENDLY
            % OSA = datatype(default_data_type)
            % sets the datatype to use
            % Currently supported:
            % single
            % double
            % uint8
            
            if nargin < 2
                default_data_type = OSA.default_data_type;
            else
                OSA.default_data_type = default_data_type;
            end
            
            switch default_data_type
                case 'single'
                    OSA.array_ = single(OSA.array_);
                    OSA.clearGPU;
                    OSA.send2GPU;
                    
                case 'double'
                    OSA.array_ = double(OSA.array_);
                    OSA.clearGPU;
                    OSA.send2GPU;
                                        
                case 'uint8'
                    OSA.array_ = uint8(OSA.array_);
                    OSA.clearGPU;
                    OSA.send2GPU;
                                        
                otherwise
                    error('I do not understand that data type (yet)!');
                    
            end
            
        end % of setdatatype
        
        function t = stoptiming(OSA)
            % t = stoptiming(OSA)
            % Returns time since OSA.time_
            t = toc;
        end
        
        
            %% GPU Utilities
        function OSA = initGPU(OSA,device_num)
            % OSA = initGPU(OSA,device_num)
            % Detects if the computer has an CUDA capable GPU, and enables
            % it for usage.
            
            if OSA.NGPUs < 1
                warning('Computer does not have NVIDIA CUDA capable GPU! Not Initializing....');
                OSA.GPU_device = 0;
                return;
            end
            
            % Assume GPU 1 if no device is input
            if nargin < 2
                device_num = 1;
            end
            
            OSA.GPU_device = gpuDevice(device_num);
            OSA.useGPU = true;
        end % of initGPU
        
        function OSA = clearGPU(OSA)
            % OSA = clearGPU(OSA)
            % Clears GPU properties
            
            if OSA.useGPU == true
                OSA.GPUarray_ = [];
                OSA.GPUfftarray_ = [];
            else
                warning('GPU is not being used/not available!');
            end
        end % of clearGPU
        
        function OSA = send2GPU(OSA,nuarray)
            % OSA = send2GPU(OSA,nuarray)
            % Sends the array_ to GPUarray_ with no argument, sets nuarray
            % to array_, and then send that to the GPU
            if OSA.useGPU == true
                if nargin < 2
                    OSA.GPUarray_ = OSA.array_;
                else
                    OSA.array(nuarray);
                    OSA.GPUarray_ = OSA.array_;
                end
            else
                warning('GPU is not being used/not available!');
            end
        end % of send2GPU
        
        function OSA = gather(OSA,~)
            % a = gather(OSA,fftflag)
            % Gathers array from the GPU back to the CPU. If no second
            % argument is given, pull from GPUarray_. If second argument is
            % given, pull from GPUfftarray_. If called by the user, it will
            % probe what is in the arrays. If called by array method,
            % stores gathered matrix into corresponding CPU arrays.
            if OSA.useGPU == true
                if nargin < 2
                    OSA.array_ = gather(OSA.GPUarray_);
                else
                    OSA.fftarray_ = gather(OSA.GPUfftarray_);
                end
            else
                warning('GPU is not being used/not available!');
            end
        end % of gather
        %% Math Operators
        
        
        %% Coordinate System
        
        % See coords utitility in @OptSysArray directory
        function [X,Y] = COORDS(OSA,local)
            % [X,Y] = COORDS(OSA,local)
            % Returns the 2D coordinates of the Array OSA
            % COORDS method by JLCodona in AOGrid
            if(nargin>1)
                [x,y] = coords(OSA,local);
            else
                [x,y] = coords(OSA);
            end
            
             % Okay.  Try to not have to compute this.
            
            % see if anything has changed since the caching...
            if(OSA.Nx_==length(x) && OSA.Ny_==length(y))
                
                %A.Xextremes==x([1 end])
                %A.Yextremes==y([1 end])
                
                % yeah, I know this is ugly.  What I really want is a call
                % like...
                % if(&&(Boolean_Array)) ...
                % If you are going to complain about elegance, provide a
                % better solution.
                if(prod(double([(OSA.Xextremes==x([1 end])) (OSA.Yextremes==y([1 end]))]))~=0)
                    % Looks like nothing has changed.  Use the cache!
                    X = OSA.X_;
                    Y = OSA.Y_;
                    fprintf('<Using COORDS cached values.>\n');
                    return;
                else
                    % fprintf('<COORDS Cache failed values test. x(%g %g) y(%g %g)>\n',...
                    % A.Xextremes,x([1 end]),A.Yextremes,y([1 end]));
                end
            else
                % fprintf('<COORDS Cache failed length test. x(%g %g) y(%g %g)>\n',...
                % A.Nx_,length(x),A.Ny_,length(y));
            end
            % Try to speed this line up...
            %[X,Y] = meshgrid(x,y);
            % Note that coords returns row vectors
            
            %fprintf('DEBUG: Computing COORDS for %s <%s>\n',class(A),A.name);
            
            % Go ahead and compute it.
            X = ones(length(y),1)*x;
            Y = y'*ones(1,length(x));
            
            % save the cached coordinates.
            OSA.X_ = X;
            OSA.Y_ = Y;
            OSA.Nx_ = length(x);
            OSA.Ny_ = length(y);
            
            OSA.Xextremes = x([1 end]);
            OSA.Yextremes = y([1 end]);
        end
            
        
        %% Interpolation
        
                
    end % of methods
    
    
    %% Static Methods
    % by JLCodona in AOGrid via AOSim2
    methods(Static=true)
        function yn = differ(mat1,mat2)
            % This is to work around an apparent bug in MATLAB.
            
            yn = prod(double(mat1(:)==mat2(:)))==0;
            
            %             if(prod(double(mat1(:)==mat2(:))))
            %                 yn = false;
            %             else
            %                 yn = true;
            %             end
            
        end
        
        function copy = copyobj(obj)
            % Create a shallow copy of the calling object.
            copy = eval(class(obj));
            meta = eval(['?',class(obj)]);
            for p = 1: size(meta.Properties,1)
                pname = meta.Properties{p}.Name;
                try
                    eval(['copy.',pname,' = obj.',pname,';']);
                catch this
                    fprintf(['\nCould not copy ',pname,' ',this ,'.\n']);
                end
            end
        end
        
        function org = middlePixel(n)
            % STATIC: org = AOGrid.middlePixel(n)
            % Figure out which pixel would be the origin in an FFT of
            % length n.
            if(isa(n,'OptSysArray'))
                n = n.size;
            end
            org = (n+2-mod(n,2))/2;
        end
        
    end % of static methods
    
end

