function BuildInterface(Name,RootPath,Args)
arguments
	Name
	RootPath			= './'
	Args.PkgName		= Name
	Args.InstallPrefix	= RootPath
end
copyfile([RootPath '/*' Args.PkgName '*'],'./');
build(eval(['define' Args.PkgName]));
if ~isempty(Args.InstallPrefix)
	copyfile(['./' Args.PkgName '/*'],[Args.InstallPrefix '/lib']);
end
end