# Software and Tools

## Widgets

Here is a list of interactive widgets that we'll go over during the lectures:

* [Slab line formation](https://tiagopereira.space/widgets/slab.html)


## Programming language

This course can be followed using the [Julia](https://julialang.org/) programming language. Why Julia? Julia is a new language designed from scratch for fast scientific calculations, and has a large and growing ecosystem, especially for astronomy, visualisation, and data science. If you are familiar with Python, using Julia will be similar and in AST4310 we will use simple interfaces to the language. Read our [guide to Julia when coming from Python](julia.md).

Besides the programming language used to run code, jupyter notebooks have a text component using the Markdown language. It will be used extensively in your project reports, and you are encouraged to get more familiar with it. The [Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) is a good reference. Familiarity with LaTeX is also important, but mostly to write equations, which can be included in Markdown using LaTeX syntax.

## JupyterHub

It may be possible to run Julia from [JupyterHub at UiO](https://jupyterhub.uio.no/). Check this page again later for updates.

## Installing Julia and necessary packages

Julia is open source and available for most operative systems. It is already installed in the linux machines at the Institute, and it is possible to connect to this system remotely with Jupyter in your browser (see below for details). We recommend you use Jupyter lab or VSCode/VSCodium for interacting with Julia notebooks.

Start by installing Julia itself:

=== "macOS or Linux laptop"
    Open a terminal and run:
    ```bash
    curl -fsSL https://install.julialang.org | sh
    ```

=== "Windows laptop"
    Install `juliaup` from the [Microsoft Store](https://www.microsoft.com/store/apps/9NJNWW8PVKMN). This can also be installed from the terminal:
    ```bash
    > winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
    ```

=== "Institute linux machines"
    Julia is already installed in the Institute's linux machines, but its module needs to be loaded from the terminal:
    ```bash
    $ module load julia
    ```

Then start `julia` in the terminal:

```bash
$ julia
```

This should start the Julia REPL, and you'll see a prompt with `julia> `. Press the `]` key to enter the Julia package manager, and the prompt should change to something like `(@v1.11) pkg> `. Here, you will install several basic packages with the command:

```
(@v1.11) pkg> add CairoMakie Colors HDF5 IJulia Interpolations NumericalIntegration PhysicalConstants  Unitful 
```

One of the packages above is `IJulia`, which will install Julia kernels for Jupyter. These should show when you restart Jupyter lab or Jupyter notebook. For VSCode or VSCodium, it is recommended you install the extensions `Julia` and `Jupyter`.


### Running Julia remotely from the Institute

If you are using Julia from the Institute's linux machines, it is possible to connect remotely via the browser in your laptop.

!!! warning 
    This procedure is technically more difficult and prone to failure, so it is recommended for advanced terminal users only.

Setting up Julia at the Institute is very similar to the setup of Python. Start by setting up your Python environment following the [documentation](https://www.mn.uio.no/astro/english/services/it/help/programming/using-python.html), and ensure that your `py313` environment is loaded. Follow [section 2](https://www.mn.uio.no/astro/english/services/it/help/programming/using-python.html#jupyter) to start jupyter lab remotely. When you start Jupyter, you should see a Julia kernel when you start a new notebook.

To ssh to your machine of choice, it is recommended to configure your ssh keys and ssh agent to avoid having to enter your password several times. Adding a proxy to your `~/.ssh/config` makes connecting even easier:

```bash
Host my_machine
    ForwardX11 yes
    ProxyJump username@login.astro.uio.no
```
