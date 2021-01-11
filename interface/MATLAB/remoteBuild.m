clc; clear all;
cluster = parcluster();
RootPath = getenv('PROJDIR');
%% Generate defineDimer.m/mlx
job1 = batch(cluster, @ReadHeader, 0, ...
	{'Dimer', RootPath}, ...
	'AutoAddClientPath',false,...
	'AutoAttachFiles',false,...
	'Pool', 0);
%% Build defineDimer.m/mlx
%  Make sure to uncomment an adjust the MATLAB interface methods
job2 = batch(cluster, @BuildInterface, 0, ...
	{'Dimer', [RootPath '/interface/MATLAB'], ...
		"InstallPrefix", RootPath}, ...
	'AutoAddClientPath',false,...
	'AutoAttachFiles',false,...
	'Pool', 0);
%% Test the interface functions
%  The CMake project includes simple C++ native tests. Compare the results
%  of this test with the ones in /test/TestHamil.cpp
%  TODO: automate this process and compare the results
job3 = batch(cluster, 'TestInterface', ...
	'AdditionalPaths', string([RootPath '/lib']), ...
	'AutoAddClientPath',false,...
	'AutoAttachFiles',false,...
	'Pool', 0);
load(job3);