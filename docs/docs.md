# Documentation

This website is a collection of all MarkDown (MD) files found in the `coconut` package.
These files are converted to a static website using [MkDocs][1], based on the outline defined in the `mkdocs.yml` file.


## Documentation in MarkDown 

MD files are easy to write and read in raw form, and can easily be converted to HTML for nicer rendering. It is easy to learn, see e.g. this [MarkDown tutorial](https://www.markdowntutorial.com/). Several flavors of MD exist, which differ only in some details. GitHub has its own [GitHub Flavored Markdown (GFM)](https://github.github.com/gfm/), which is used to render MD documents found in repositories on the run. PyCharm comes with an automatic rendering tool built in, and I assume many other IDE's too. 

The MD documents for the CoCoNuT documentation must be written in the [Python-Markdown](https://python-markdown.github.io/) flavor, because they are processed by [MkDocs][1]. This flavor stays very close to the *reference implementation* of MD. Pay attention to lists: they should be indented with 4 spaces (contrary to GFM). 

An important rule for writing MD files for this documentation website is that the name of the file must be unique: use names like `fluent.md` or `examples.md`, not `readme.md`.


### Links to other MarkDown files
It is possible to use relative links to other MarkDown files in CoCoNuT, using the syntax
```markdown
[link description](relative_path)
```
where `relative_path` is the relative path to another MarkDown file, e.g. `../coupling_components/mappers/mappers.md`. 

These links can be used in rendered MarkDown, e.g. in PyCharm, but also on GitHub itself (see [this blogpost](https://github.blog/2013-01-31-relative-links-in-markup-files/)). 
These links also work on the documentation website, as they are automatically replaced by the correct URL. 
Take for example a look at the documentation of the [mappers](../coupling_components/mappers/mappers.md) or the [examples](../examples/examples.md). 

In addition, it is also possible to link to (sub)sections of a file. It is less straightforward to find the correct address in this case, but PyCharm gives you suggestions if you type a # behind the filename. Some examples: [this paragraph](docs.md#links-to-other-markdown-files), [a CoCoNuT tutorial](../README.md#getting-started), [setting up a Fluent case](../coupling_components/solver_wrappers/fluent/fluent.md#setting-up-a-new-case), [radial basis mappers](../coupling_components/mappers/mappers.md#mapperradialbasis). This feature has some limitations: it is not possible to link to a section name that appears several times in a file (non-unique name), or to a section name that is a link itself.



### Math
LaTeX notation can be used for writing mathematical expressions. Inline equations must be enclosed in single dollar signs (e.g. $E = m c^2$), block-style equations in double dolar signs, e.g.

$$
e^{i \pi} + 1 = 0.
$$

LaTeX expressions will *not* be rendered on GitHub, but only on the documentation website. For the latter, the MD extension [Arithmatex](https://facelessuser.github.io/pymdown-extensions/extensions/arithmatex/) is used to render the expressions with MathJax. Note that [MathJax syntax](https://math.meta.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference/) is a little more restrictive than a real LaTeX installation. 

### Images
External images can be included with their URL. Adding locally stored images is a bit more complicated: these images must be stored in a directory `images` next to the MD file. If another location is used, they will not be shown on the website, only on GitHub. Furthermore, images must have a unique name. A warning is shown when this is not the case. 

An image can be added with the MD command

```markdown
![alt](images/example.png "description")
```

with `alt` displayed when the image cannot be shown/viewed for some reason, and `description` appearing when hovering over the image with your mouse. For example:

![example image](images/lachend_kakske.png "diefstal in Oostende")

Only image formats specified in `run_mkdocs.py` (i.e. png, jpg, jpeg, gif, svg) are copied to the website; missing extensions can be added. 
Images from all `coconut` subdirectories called `images` are copied to the website, so care must be taken that `images` is not used in e.g. the output of the test examples.

### Style & layout guide

This section gives some guidelines about style and layout of MarkDown files, to keep the documentation consistent.

* Use code style for:

    * class and method names (and plurals): `Model`, `ModelParts`,`__init__`,`finalize`
    * JSON keys and values: `coupled_solver`,`delta_t`
    
* Use code style + italics for:

    * files: `run_simulation.py`, `parameters.json`
    * folders: `data_structure`
    * paths: `coupling_components/solver_wrappers/mapped.py`

* Use normal text for:

    * referring to abstract terms and concepts in CoCoNuT (i.e. not a specific class): solver wrappers, mappers, coupled solver, data structure
    
* Titles of MD files (e.g. `# Mappers`, the first line of the MarkDown file):

    * should be brief and not repeat information that can be deducted from the structure of the documentation, e.g. for the Fluent solver wrapper: just use `# Fluent` and not `# Fluent solver wrapper`, as it is beneath `Solver wrappers` on the website.
    * don't use class names (and hence no *camelcase*), e.g. not `# SolverWrapperOpenFOAM`
    
* If you refer to other MarkDown pages in the documentation, it can be useful to use a [relative link](docs.md#links-to-other-markdown-files).

* Recommendation for links: it is nice that the link text gives you some information about where the link goes, for example: 

    * good example: [coconut documentation][2]
    * bad example: [this link][3]
    
    
## Creating a static website with MkDocs

[MkDocs][1] can be installed using pip: 
```bash
pip install mkdocs-material
```
Apart from MkDocs, this also automatically installs the [Material theme](https://squidfunk.github.io/mkdocs-material/) and the [PyMdown](https://squidfunk.github.io/mkdocs-material/extensions/pymdown/) MD extensions (including Arithmatex). 

To install locally, add a relative or absolute path with the `-t` argument:
```bash
pip install mkdocs-material -t /some/absolute/path
```
In this case, two lines must be added to your `.bashrc` file: 
```bash
export PYTHONPATH=/some/absolute/path:$PYTHONPATH
export PATH=/some/absolute/path/bin:$PATH
```

The structure/outline of the website is dictated by the `nav` variable in `mkdocs.yml`. This is the only variable that must be adjusted when new MD files are added to the code.

The complete process to create the documentation website is automated by `run_mkdocs.py`. This does the following things:

*   Copy all MD files in `coconut` and its subfolders to folder `docs`. 
*   Check if there are duplicate filenames: these overwrite each other! If duplicates are found, a warning is given with the original paths.
*   Check if each MD file is mentioned in `mkdocs.yml`. If a file is not mentioned, a warning is given.
*   Build static HTML website using `mkdocs build`. 

You can run `run_mkdocs.py` with an extra command line argument: 

```bash
python run_mkdocs.py --preview example
```
opens a preview of the website in Firefox, showing the webpage corresponding to the file `example.md`. This can be used to check MD and LaTeX syntax. 


```bash
python run_mkdocs.py --deploy
```
deploys the website on GitHub Pages using `mkdocs gh-deploy`. It seems this is currently only possible if the remote is configured with SSH, not with HTTPS.  




[1]:    https://www.mkdocs.org/
[2]:    https://pyfsi.github.io/coconut/
[3]:    https://www.youtube.com/watch?v=w0AOGeqOnFY