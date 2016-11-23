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
    % Authors: Alexander T. Rodack and Johanan L. Codona
    % Date: 11/14/2016 - 
    %
    % Special Thanks to JLCodona. This will bear significant resemblance
    % to the AOGrid class within AOSim2: 
    % https://github.com/jlcodona/AOSim2/tree/SCDA
    
    % TO DO:
    % 1) finish basic functionality
    % 2) make multiple GPUs usable
    %       a) Perhaps with cell arrays?
    
    
    %% Properties
    
    % Public Properties
    properties(GetAccess='public',SetAccess='public')
        
        % Properties
        name; % name the object internally
        default_array_size = [1,1]; % set a default size for an array to use if no input arg is given to contstructor
        FFTsize = [0,0]; % default size of the fftgrid
        Offset = [0 0]; % moves the object.
        verbose = false; % prints info and plots
        
        % Miscellaneous Properties
        interpolation_method = []; % selects a method of interpolation
        seed = []; % Set seed value for repeatable random numbers
        
    end
    
    
    properties(GetAccess='public',SetAccess='protected')
        
        % Storage arrays
        array_;
        persistedstorage_ = [];
        GPUarraystorage_ = [1,1,2];
        fftarray_ = [];
        
        % Storage arrays on the GPU
        GPUarray_;
        GPUfftarray_ = [];
        GPUpersisted_ = [];
        
        % GPU properties
        NGPUs_ = 0; % default to 0, but is initialized in constructor
        GPU_device_; % which GPU is being used
        useGPU_ = false; % flag to use a detected GPU or use the CPU
        % set to true if GPU is detected, can be turned off
        % if desired using disableGPU().
        
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
        R_ = []; % cached radial coordinate (no check, always recomputed)
        THETA_ = []; % cached angular coordinate (always recomputed)
        
        % Stop Watch Properties
        starttime_ = [];
        stoptime_ = [];
        exectime_ = [];
        
    end
    
    %% Methods
    methods
        %% Constructor
        function OSA = OptSysArray(nxy)
            % OSA = OptSysArray(nxy)
            % This is the contructor function for the OptSysArray class. The
            % input argument can be of several types: A scalar value, a
            % two-vector, or another OptSysArray. Automatically searches
            % for CUDA enabled GPUs, and initializes them for use if found.
            
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
            
            % Initialize GPUs
            OSA.initGPU(1,2); % give inputs as if 2 GPUs are available, so if 2 are, they are both used. If < 2 are available, the function handles it accordingly
            
            % Set default data type
            OSA.setdatatype();
            
        end % of Constructor
        
        
        %% Read Out Utilities
        function description = describe(OSA)
            % description = describe(OSA)
            % Prints a summary of the important object parameters:
            % name, {class [size of array] [coordinate spacing] data type}
            % {GPU name, GPU device number}
            
            if OSA.NGPUs_ < 1
                description = sprintf('%s {%s: [%dx%d] [%g,%g] %s}',OSA.name,class(OSA),OSA.nx,OSA.ny,OSA.dx,OSA.dy, OSA.default_data_type);
            elseif OSA.NGPUs_ == 1
                description = sprintf('%s {%s: [%dx%d] [%g,%g] %s} {%s : Device %d} ',OSA.name,class(OSA),OSA.nx,OSA.ny,OSA.dx,OSA.dy,OSA.default_data_type,OSA.GPU_device_{1}.Name,OSA.GPU_device_{1}.Index);
            elseif OSA.NGPUs_ == 2
                description = sprintf('%s {%s: [%dx%d] [%g,%g] %s} {%s : Device %d} {%s : Device %d} ',OSA.name,class(OSA),OSA.nx,OSA.ny,OSA.dx,OSA.dy,OSA.default_data_type,OSA.GPU_device_{1}.Name,OSA.GPU_device_{1}.Index,OSA.GPU_device_{2}.Name,OSA.GPU_device_{2}.Index);
            else
                
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
        end % of dy
        
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
        % Not inclusive of all protected properties. Most are set in
        % utilities that are better categorized in other sections
                
        function OSA = setdefaultdatatype(OSA,datatype)
            % OSA = setdefaultdatatype(OSA,datatype)
            % Sets the default_data_type property that is used in the
            % calling of OSA.array(nuarray). datatype must be one of the
            % following strings exactly:
            % 1) 'single'
            % 2) 'double'
            % 3) 'uint8'
            
            OSA.default_data_type = datatype;
            OSA.setdatatype;
            
        end % of setdefaultdatatype
        
        
        
        function a = array(OSA,nuarray,mask)
            % a = array(OSA,nuarray,mask)
            % Multi-functional call
            % 1) Calling with no input arguments returns what is in array_
            % 2) Set the array to a new one, return the new OSA object
            %
            % method by JLCodona in AOGrid with edits for data type and GPU
            % compatibility by ATRodack
            
            if nargin == 1
                a = OSA.array_;
            else % Set the array_ to the input
                if nargin == 2
                    nuarray = squeeze(nuarray);
                    nuarray = nuarray(:,:,1);
                    nuarray = squeeze(nuarray);
                    if numel(OSA) == numel(nuarray)
                        OSA.array_ = nuarray(:);
                        OSA.setdatatype();
                        OSA.send2GPU;
                        OSA.touch;
                    else
                        OSA.resize(size(nuarray));
                        OSA.array_ = nuarray;
                        OSA.setdatatype();
                        OSA.send2GPU;
                        OSA.touch;
                    end
                else
                    OSA.array_(mask(:)) = nuarray(:);
                end
                a = OSA;
            end
        end % of array
        
        function OSA = clearCPU(OSA)
            % OSA = clearCPU(OSA)
            % Clears the array_ property. Mostly here for testing GPU
            % gather method, but it could be useful in other respects
            
            OSA.array_ = [];
        end % of clearCPU
        
        
        %% Micellaneous Utilities
        
        function OSA = touch(OSA)
            % OSA = touch(OSA)
            % Clears the fft cache
            %
            % method by JLCodona in AOGrid
            
            OSA.fftarray_ = [];
            OSA.GPUfftarray_ = [];
        end % of touch
        
        
        function newsize = resize(OSA,varargin)
            % newsize = OSA.resize(newsize);
            % Change the size of an OptSysArray object's data array.
            % (Assigning a new data array does this automatically.)
            % Note that resize leaves the array spacing unchanged.
            % Use OSA.spacing(dx) to set the new spacing.
            %
            % method by JLCodona in AOGrid
            
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
        
        
        function OSA = setdatatype(OSA,default_data_type)
            % OSA = datatype(default_data_type)
            % Sets the datatype to use. Do not do anything if already of
            % the correct data type.
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
                    if ~isa(OSA.array_,'single')
                        OSA.array_ = single(OSA.array_);
                    end
                    
                case 'double'
                    if ~isa(OSA.array_,'double')
                        OSA.array_ = double(OSA.array_);
                    end
                    
                case 'uint8'
                    if ~isa(OSA.array_,'uint8')
                        OSA.array_ = uint8(OSA.array_);
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of setdatatype
        
        function t = starttime(OSA)
            % OSA = starttime(OSA)
            % Starts a stopwatch
            
            t = clock;
            OSA.starttime_ = t;
        end % of starttiming
        
        function t = stoptime(OSA)
            % t = stoptime(OSA)
            % Return the time the stop watch is stopped. Also sets the time
            % difference between this time and when the stop watch was
            % started in exectime_.
            
            if isempty(OSA.starttime_)
                warning('Stopwatch not started');
                return;
            end
            t = clock;
            OSA.stoptime_ = t;
            
            % Compute execution time
            OSA.exectime_ = etime(OSA.stoptime_,OSA.starttime_);
            % tic toc is more reliable for timing, but seems difficult to
            % implement using properties. This is a good work around, and
            % if you need to know the time more accurately/reliably than
            % this, use tic toc in outside calls or add it in to the method
            % you need to time.
            
        end
        
        
        %% GPU Utilities
        function OSA = initGPU(OSA,device_num_calc, device_num_persist)
            %  OSA = initGPU(OSA,device_num_calc, device_num_persist)
            % Detects if the computer has an CUDA capable GPU, and enables
            % it for use. If more than one GPU is detected, allow for a
            % second one to be used as well. In this case, the first GPU is
            % set to be the active one.
            
            % Count the number of GPUs
            OSA.NGPUs_ = gpuDeviceCount;
            
            % Return if no GPUs are detected
            if OSA.NGPUs_ < 1
                warning('GPU:noDevice','Computer does not have NVIDIA CUDA capable GPU! Not Initializing....');
                OSA.GPU_device_ = 0;
                return;
            end
            
            % Initialize cell array to hold GPU Device Structures
            OSA.GPU_device_ = cell(1,OSA.NGPUs_);
            
            % Initialize the GPUs
            if OSA.NGPUs_ == 1 % One GPU detected
                
                if nargin == 1 % no input, don't use a GPU (user surely knows what they are doing, right?)
                    OSA.useGPU_ = false;
                else % use the detected GPU, ignore any 3rd input
                    device_num_calc = 1; % Don't allow user to supply incorrect input if only 1 GPU is present
                    OSA.GPU_device_{1} = gpuDevice(device_num_calc);
                    OSA.useGPU_ = true;
                end
                
            elseif OSA.NGPUs_ > 1 % Multiple GPUs detected
                
                if nargin == 1 % don't initialize a GPU
                    OSA.useGPU_ = false;
                    
                elseif nargin == 2 % use 1 GPU
                    OSA.GPU_device_{1} = gpuDevice(device_num_calc);
                    OSA.useGPU_ = true;
                    
                elseif nargin == 3 % use 2 GPUs
                    OSA.GPU_device_{1} = gpuDevice(device_num_calc);
                    OSA.GPU_device_{2} = gpuDevice(device_num_persist);
                    
                    gpuDevice(OSA.GPU_device_{1}); % set calc GPU to active
                    OSA.useGPU_ = [true, true];
                end
                
            end
            
        end % of initGPU
        
        function OSA = clearGPU(OSA)
            % OSA = clearGPU(OSA)
            % Clears calculation GPU properties. Does not clear persistent
            % GPU.
            
            if prod(double(OSA.useGPU_)) == 1
                % cache current matrices
                OSA.GPUarraystorage_(:,:,1) = OSA.gather(0);
                OSA.GPUarraystorage_(:,:,2) = OSA.gather(1);
                
                % make sure active GPU is GPU 1
                
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
            
            if prod(double(OSA.useGPU_)) == 1
                if OSA.NGPUs_ < 2
                    
                    % check that GPU 1 is active
                    if gpuDevice ~= 1
                        gpuDevice(OSA.GPU_device_{1});
                        % NOTE: If this is done, even if device 1 was
                        % already active, it resets all the gpuArrays.
                        % clearGPU() is called after this step anyway, so
                        % it shouldn't matter here.
                    end
                    
                    OSA.clearGPU;
                    if nargin < 2
                        OSA.GPUarray_ = gpuArray(OSA.array_);
                    else
                        OSA.array(nuarray);
                        OSA.GPUarray_ = gpuArray(OSA.array_);
                    end
                elseif OSA.NGPUs_ >= 2
                    
                    
                    
                end
            else
                warning('GPU is not being used/not available!');
            end
        end % of send2GPU
        
        function OSA = send2persistentGPU(OSA,M)
            % OSA = send2persistentGPU(OSA,M)
            % Sends a Matrix M to the GPU Device being used for
            % persistence.
            
            if prod(double(OSA.useGPU_)) == 1
                % set active GPU to persisted if necessary
                
                
                
                OSA.persistedstorage_ = M;
                OSA.GPUpersisted_ = gpuArray(OSA.persistedstorage_);
            end
        end % of send2persistentGPU
        
        function OSA = gather(OSA,fftflag)
            % a = gather(OSA,fftflag)
            % Gathers array from the GPU back to the CPU. If no second
            % argument is given, pull from GPUarray_. If second argument is
            % given, pull from GPUfftarray_. If called by the user, it will
            % probe what is in the arrays. If called by array method,
            % stores gathered matrix into corresponding CPU arrays.
            
            if prod(double(OSA.useGPU_)) == 1
                if nargin < 2
                    OSA.array_ = gather(OSA.GPUarray_);
                else
                    if fftflag == 1
                        OSA.fftarray_ = gather(OSA.GPUfftarray_);
                    else
                        OSA.array_ = gather(OSA.GPUarray_);
                    end
                end
            else
                warning('GPU is not being used/not available!');
            end
        end % of gather
        
        function OSA = switchGPU(OSA,device_num)
            % OSA = switchGPUs(OSA,device_num)
            % If more than 1 GPU is present, switch which one is being used
            
            if(OSA.NGPUs_ > 1)
                if(device_num <= OSA.NGPUs_)
                    OSA.GPU_device_ = gpuDevice(device_num);
                else
                    warning('GPU:noDevice','No GPU Device %d available',device_num);
                    OSA.GPU_device_ = gpuDevice(1);
                    return;
                end
            else
                warning('GPU:noDevice','No GPU Device %d available',device_num);
                OSA.GPU_device_ = gpuDevice(1);
                return;
            end
        end % of switchGPU
        
        function OSA = disableGPU(OSA)
            % OSA = disableGPU(OSA)
            % Disables use of GPU
            
            OSA.clearGPU;
            OSA.useGPU_ = false;
        end % of disableGPU
        
        function OSA = enableGPU(OSA)
            % OSA = enableGPU(OSA)
            
            OSA.useGPU_ = true;
        end % of enableGPU
        
        %% Math Operators
        
        
        %% Coordinate System
        
        % See coords utitility in @OptSysArray directory
        function [X,Y] = COORDS(OSA,local)
            % [X,Y] = COORDS(OSA,local)
            % Returns the 2D coordinates of the Array OSA
            %
            % method by JLCodona in AOGrid
            
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
                    %fprintf('<Using COORDS cached values.>\n');
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
        end % of COORDS
        
        function [THETA,R] = POLARCOORDS(OSA)
            % [THETA,R] = POLARCOORDS(OSA)
            % Returns the polar coordinates computed from the cached
            % Cartesian Coordinates
            
            % Run COORDS....If it has already been done, just use the
            % cached coordinate system
            OSA.COORDS();
            [THETA,R] = cart2pol(OSA.X_,OSA.Y_);
            
            % save to cache
            OSA.R_ = R;
            OSA.THETA_ = THETA;
        end % of POLARCOORDS
        
        
        function o = origin(OSA,varargin)
            % o = origin(OSA,varargin)
            % Sets the real world "center" of the coordinate system defined
            % by coords().
            %
            % method by JLCodona in AOGrid
            
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
        
    %% Interpolation
        
        
    end % of methods
    
    
    %% Static Methods
    % All methods by JLCodona in AOGrid via AOSim2
    methods(Static=true)
        function yn = differ(mat1,mat2)
            % This is to work around an apparent bug in MATLAB.
            
            yn = prod(double(mat1(:)==mat2(:)))==0;
            
            %             if(prod(double(mat1(:)==mat2(:))))
            %                 yn = false;
            %             else
            %                 yn = true;
            %             end
            
        end % of differ
        
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
        end % of copy
        
        function org = middlePixel(n)
            % STATIC: org = AOGrid.middlePixel(n)
            % Figure out which pixel would be the origin in an FFT of
            % length n.
            if(isa(n,'OptSysArray'))
                n = n.size;
            end
            org = (n+2-mod(n,2))/2;
        end % of middlePixel
        
    end % of static methods
    
end

