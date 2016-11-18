function result = TestOSAClass(test_case)
% result = TestOSAClass(test_case)
% Function to test the functionality of the OptSysArray class as it is
% contructed. Input "test_case" determines what to test


%% Test OSA class
switch test_case
    case 1
        nxy = 512;
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 1';
        result = OSA.describe;
    case 2
        nxy = [512,1024];
        OSA = OptSysArray(nxy);
        OSA.name = 'Test Case 2';
        result = OSA.describe;
    case 3
        OSA = OptSysArray(512);
        OSA.name = 'Test Case 3';
        OSA2 = OptSysArray(OSA);
        result = OSA2.describe;

end









end
