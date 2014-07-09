import re
File=input("file name:\n")
handle=open(File,"r")
number=len(re.findall("\>",handle.read()))
print(number,"sequences.")
