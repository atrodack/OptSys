function result = TestOSArray(test_case)
% result = TestOSAClass(test_case)
% Function to test the functionality of the OptSysArray class as it is
% contructed. Input "test_case" determines what to test

%% Test OSA class
switch test_case
    
    case 1
        %% Test Constructor for single input
        % Use Describe to also test initializations, dx, dy, nx, ny
        fprintf('**************\n    Case 1 \n');
        nxy = 512;
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 1';
        result = OSA.describe;
        device = gpuDevice(1);
        teststring = sprintf('Test Case 1 {OptSysArray: [512x512] [0.01,0.01] single} {%s : Device 1}',device.Name);
        test = strncmp(teststring,result,length(result)-1);
        if test == 1
            fprintf('Case 1 passed\n');
        else
            fprintf('Case 1 failed\n');
        end
    case 2 
        %% Test Constructor for 2-vector input
        fprintf('**************\n    Case 2 \n');
        nxy = [512,1024];
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 2';
        result = OSA.describe;
        device = gpuDevice(1);
        teststring = sprintf('Test Case 2 {OptSysArray: [1024x512] [0.01,0.01] single} {%s : Device 1}',device.Name);
        test = strncmp(teststring,result,length(result)-1);
        if test == 1
            fprintf('Case 2 passed\n')
        else
            fprintf('Case 2 failed\n');
        end
    case 3 
        %% Test Constructor for OptSysArray input
        fprintf('**************\n    Case 3 \n');
        OSA = OptSysArray(512);
        OSA.name = 'Test Case 3';
        OSA2 = OptSysArray(OSA);
        result = OSA2.describe;
        device = gpuDevice(1);
        teststring = sprintf('copy of Test Case 3 {OptSysArray: [512x512] [0.01,0.01] single} {%s : Device 1}',device.Name);
        test = strncmp(teststring,result,length(result)-1);
        if test == 1
            fprintf('Case 3 passed\n');
        else
            fprintf('Case 2 failed\n');
        end
    case 4 
        %% Test array() and resize
        fprintf('**************\n    Case 4 \n');
        a = magic(6);
        OSA = OptSysArray(512);
        OSA.array(a);
        test = isequal(a(:),OSA.array_(:));
        if test == 1
            fprintf('Case 4 passed\n');
            result = true;
        else
            fprintf('Case 4 failed\n');
            result = false;
        end
    case 5 
        %% Test setdatatype
        fprintf('**************\n    Case 5 \n');
        a = magic(6);
        if isa(a,'double') == 1
            type1 = 'double';
        end
        OSA = OptSysArray(6);
        type2 = OSA.default_data_type;
        if strncmp(type1,type2,4) == 0
            OSA.array(a);
            test = isa(OSA.array_,type2);
        end
        if test == 1
            fprintf('Case 5 passed\n');
            result = true;
        else
            fprintf('Case 5 failed\n');
            result = false;
        end
    case 6 
        %% Test GPU utilities
        fprintf('**************\n    Case 6 \n');
            a = magic(6);
            OSA = OptSysArray(6);
            OSA.array(a);
        if(OSA.NGPUs_ > 0)
            % test 1: send2GPU
            OSA.send2GPU;
            test1 = isequal(OSA.array_,OSA.GPUarray_);
            if test1 == 0
                fprintf('send2GPU test failed\n');
            end
            
            % test 2: clearGPU
            OSA.clearGPU;
            test2 = isempty(OSA.GPUarray_);
            if test2 == 0
                fprintf('clearGPU test failed\n');
            end
            
            % test 3: gather
            OSA.send2GPU(a); % last test cleared it, send it back
            OSA.clearCPU; % clear array_
            OSA.gather; % write GPU back to CPU
            test3 = isequal(OSA.array_,a);
            if test3 == 0
                fprintf('gather test failed\n');
            end
            
            % test 4: switchGPU
            warning('off','GPU:noDevice')
            if OSA.NGPUs_ > 1
                OSA.switchGPU(2);
                test4 = isequal(OSA.GPU_device_.Index,2);
                if test4 == 0
                    fprintf('switchGPU test failed\n');
                end
            elseif OSA.NGPUs_ == 1
                OSA.switchGPU(2); % should throw warning that is suppressed,
                % and set the device back to 1
                test4 = isequal(OSA.GPU_device_.Index,1);
                if test4 == 0
                    fprintf('switchGPU test failed\n');
                end
            else
                fprintf('Cannot test\n');
                test4 = 1;
            end
            warning('on','GPU:noDevice')
            
            test = prod([test1, test2, test3, test4]);
            if test == 1
                fprintf('Case 6 passed\n');
                result = true;
            else
                fprintf('Case 6 failed\n');
                result = false;
            end
            
        else
            fprintf('Cannot test GPU, no GPU device available\n');
        end
        


end








end
