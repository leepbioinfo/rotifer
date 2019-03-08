import subprocess
s = subprocess.Popen(['man ./rlineage.man'], shell = True)
s.communicate()[0]
