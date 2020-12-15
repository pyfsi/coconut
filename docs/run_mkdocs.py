import glob
import os
import shutil
from sys import argv
import re
import fileinput

"""
README
------

Execute run_mkdocs.py from this directory to update the documentation. 
Default behavior is to build the documentation but not deploy it.

Required software:
    mkdocs            (install with: pip install mkdocs)
    material theme    (install with: pip install mkdocs-material)

Optional arguments:
    --deploy        build documentation and deploy on pyfsi.github.io/coconut/, 
                    this requires admin rights on GitHub
    --preview X     build documentation and preview webpage of file X.md, 
                    if X is not supplied, README.md is shown
                    
Example:
    python run_mkdocs.py --preview mappers
"""

# check directory
if os.getcwd() != os.path.dirname(os.path.realpath(__file__)):
    raise SystemError('execute run_mkdcos.py from its directory')

# clean docs folder
shutil.rmtree('docs', ignore_errors=True)
shutil.rmtree('site', ignore_errors=True)

os.mkdir('docs')
os.mkdir('docs/images')
os.mkdir('docs/assets')
os.mkdir('docs/assets/images')


# find all MarkDown files in CoCoNuT
files = glob.glob('../**/*.md', recursive=True)

# check for duplicate MarkDown files
filenames = []
for file in files:
    filenames.append(file.split('/')[-1])

for i, filename in enumerate(filenames):
    if filenames.count(filename) > 1:
        tools.print_info(f'WARNING - duplicate filename "{files[i]}"')

# copy all MarkDown files to docs folder
for file in files:
    shutil.copy(file, 'docs/')

# check if all MarkDown files are mentioned in nav
unused = []
used = False
for filename in filenames:
    with open('mkdocs.yml', 'r') as file:
        for line in file:
            if filename in line:
                used = True
                break
    if not used:
        unused.append(filename)
    used = False
for file in unused:
    print(f'WARNING - file "{file}" is not used in mkdocs.yml')

# fix links to MarkDown files
for filename in filenames:
    with fileinput.FileInput('docs/' + filename, inplace=True) as file:
        for line in file:
            matches = []
            for match in re.finditer(r'\([a-zA-Z0-9_./\-]+\.md\)', line):
                matches.append(match)
            if len(matches) > 0:
                for match in matches[::-1]:
                    tmp = match[0][1:-4].split('/')[-1]
                    if tmp == 'README':
                        url = '(https://pyfsi.github.io/coconut/)'
                    else:
                        url = '(https://pyfsi.github.io/coconut/' + tmp + '/)'
                    line = line[:match.start()] + url + line[match.end():]
            print(line, end='')

# find all relevant images in CoCoNuT
extensions = ['png', 'jpg', 'jpeg']
images = []
for ext in extensions:
    images += glob.glob(f'../**/images/*.{ext}', recursive=True)

# check for duplicate images
imagenames = []
for image in images:
    imagenames.append(image.split('/')[-1])

for i, imagename in enumerate(imagenames):
    if imagenames.count(imagename) > 1:
        print(f'WARNING - duplicate image "{images[i]}"')

# copy all images to docs/images folder
for image in images:
    shutil.copy(image, 'docs/images/')

# add favicons
shutil.copy('logo.png', 'docs/images/logo.png')
shutil.copy('favicon.ico', 'docs/assets/images/favicon.ico')

# build and deploy website
print('\n---MkDocs output from here---')
if len(argv) > 1:
    if argv[1] == '--preview':
        os.system('mkdocs build --clean')
        cwd = os.getcwd()
        if len(argv) == 2:
            cmd = 'firefox ' + os.path.join(cwd, 'site', 'index.html &')
        else:
            if argv[2] == 'README':
                cmd = 'firefox ' + os.path.join(cwd, 'site', 'index.html &')
            else:
                cmd = 'firefox ' + os.path.join(cwd, 'site', argv[2], 'index.html &')
        os.system(cmd)
    elif argv[1] == '--deploy':
        os.system('mkdocs gh-deploy')
    else:
        os.system('mkdocs build')
else:
    os.system('mkdocs build')
