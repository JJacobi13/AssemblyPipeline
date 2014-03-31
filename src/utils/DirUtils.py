import os, fnmatch

def fileRegexToList(regex):
    if "," in regex:
        return regex.split(",")
    elif "*" in regex:
        files = []
        for fileName in os.listdir(os.path.dirname(regex)):
            if fnmatch.fnmatch(fileName, os.path.basename(regex)):
                files.append(os.path.dirname(regex) + "/" + fileName)
        return files
    else:
        return [regex]
                