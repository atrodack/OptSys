function result = TestOSArray(test_case)
% result = TestOSAClass(test_case)
% Function to test the functionality of the OptSysArray class as it is
% contructed. Input "test_case" determines what to test
% Strncmp calls for my GPU configuration --- needs to be generalized

%% Test OSA class
switch test_case
    % Use Describe to also test initializations
    case 1 %Test Constructor for single input
        nxy = 512;
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 1';
        result = OSA.describe;
        test = strncmp('Test Case 1 {OptSysArray: [512x512] [0.01,0.01] single} {GeForce GT 750M : Device 1}',result,length(result)-1);
        if test == 1
            fprintf('Test 1 passed\n');
        else
            fprintf('Test Case 1 failed\n');
        end
    case 2 %Test Constructor for 2-vector input
        nxy = [512,1024];
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 2';
        result = OSA.describe;
        test = strncmp('Test Case 2 {OptSysArray: [1024x512] [0.01,0.01] single} {GeForce GT 750M : Device 1}',result,length(result)-1);
        if test == 1
            fprintf('Test 2 passed\n')
        else
            fprintf('Test Case 2 failed\n');
        end
    case 3 %Test Constructor for OptSysArray input
        OSA = OptSysArray(512);
        OSA.name = 'Test Case 3';
        OSA2 = OptSysArray(OSA);
        result = OSA2.describe;
        test = strncmp('copy of Test Case 3 {OptSysArray: [512x512] [0.01,0.01] single} {GeForce GT 750M : Device 1}',result,length(result)-1);
        if test == 1
            fprintf('Test 3 passed\n');
        else
            fprintf('Test Case 2 failed\n');
        end
    case 4 % Test array()
        a = magic(6);
        OSA = OptSysArray(512);
        OSA.array(a);
        test = isequal(a(:),OSA.array_(:));
        if test == 1
            fprintf('Test 4 passed\n');
            result = true;
        else
            fprintf('Test Case 4 failed\n');
            result = false;
        end
    case 5 %Test setdatatype
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
            fprintf('Test 5 passed\n');
            result = true;
        else
            fprintf('Test Case 5 failed\n');
            result = false;
        end
        
        

end









end
