-- Local variables
local myfn  = myFileName()
local mymfn = myModuleFullName()
local mydir = myfn:gsub(mymfn .. ".lua",""):gsub("/etc/lmod/modules/","")
local myshell = "bash" -- set this value to your default shell
if (os.getenv("SHELL") ~= nil) then
	myshell = os.getenv("SHELL"):gsub("(.*/)(.*)", "%2")
end
local ext = {bash = ".sh",zsh = ".zsh",tcsh = ".csh",csh = ".csh",sh = ".sh",fish = ".fish"}

-- Whatis
whatis("Name: " .. mymfn)
whatis("Path: " .. myfn)
whatis("Version: 0.1 ")
whatis("Category: data analysis ")
whatis("Description: Rotifer's shell setup.")
whatis("URL: https://github.com/robsonfsouza/rotifer ")

-- Dependencies
load("leep")

-- Environment variables
prepend_path("PYTHONPATH", mydir .. "/lib")
prepend_path("PERL5LIB", mydir .. "/perl/lib")
if (os.getenv("ROTIFER_DATA") == nil) then
	setenv("ROTIFER_DATA","/databases")
end
if (os.getenv("DATABASES") == nil) then
	setenv("DATABASES","/databases")
end

-- PATH
prepend_path("PATH", mydir .. "/bin")

-- Load shell functions
if (ext[myshell] ~= nil) then
	require "lfs"
	for script in lfs.dir(pathJoin(mydir,"etc","profile.d")) do
		if (script:sub(script:len() - ext[myshell]:len() + 1) == ext[myshell]) then
			local name = script:gsub(ext[myshell],"")
			local code = capture("cat " .. pathJoin(mydir,"etc","profile.d",script))
			if (code:match('^[ \t\n]*function *[^{]+{\n[\t ]+')) then
				code = code:gsub("\n$",""):gsub("\n}$",""):gsub('^[ \t\n]*function *[^{]+{\n[\t ]+',"")
				set_shell_function(name,code)
			end
		end
	end
end

-- Help message
local helpMsg = [[
ROTIFER
=======

Rapid Open-source Tools and Infrastructure For data Exploration and Research
----------------------------------------------------------------------------

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.

Current installation directory: ]] .. mydir .. [[
Shell: ]] .. myshell
