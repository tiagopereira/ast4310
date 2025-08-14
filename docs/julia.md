# Quick guide to Julia from Python

Julia is a programming language designed for speed and ease of use. It aims to be as fast as C and as easy to code as Python. In fact, writing Julia code can be remarkably similar to Python, but there are a few important differences to be aware. 

Here, we highlight how to perform some basic tasks with Julia, and how they compare with their Python equivalent. But first, a quick primer on what Julia is.

## Julia in nutshell

Julia is a dynamically-typed language that feels more like a scripting language (e.g. Python) than a compiled language (e.g. C, Fortran). However, code in Julia is indeed compiled, but this step is hidden from the user by using [just-in-time compilation](https://en.wikipedia.org/wiki/Just-in-time_compilation). It works similar to [Numba](https://numba.pydata.org/) in Python. A difference is that Numba can be used with a very limited set of libraries/functions, while Julia is much more general. 

Julia can be run interactively from a terminal or a Jupyter notebook. The Julia terminal programme is referred to as REPL. A julia programme can also be run as a script, or placed in a package. Julia has its own package manager that is very easy to use.

A key feature of Julia is called [multiple dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch). With multiple dispatch, Julia can split a function in several versions that are optimised for different inputs. For example, a multiplication function can have a version optimised for integers and another version optimised for floating point numbers. When a user calls that function, Julia automatically dispatches the call to the optimal function. In Julia, unlike Python, you can define the same function name multiple times and this will not cause problems, as long as each version has an input with a different data type.

## Julia for Python programmers

Read the manual on [important differences from Python](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python). Most importantly, Julia blocks are terminated by `end` keywords and indentation is not significant, and indexing starts at 1 instead of 0, and is column major (Fortran order) instead of row major (C/C++ order). 

All the array functionality and many mathematical functions are built-in into Julia, instead of being in separate packages like in Python (e.g. `math`, `numpy`, `scipy`).  

There are no classes in Julia. The closest thing are structures, which can be used to store data (similar to class attributes in Python), but have no methods. 

Here are a few translations of common operations:

Python | Julia
--------|-------
`numpy.shape()` | `size()`
`len()` | `length()` 
`type()` | `typeof()`
`range(1,5)` |   `1:4`
`numpy.linspace(0, 2, 10)` |   `range(0, 2, 10)`
`numpy.squeeze()` | `dropdims`
`import` | `using`
`my_array[-1]` | `my_array[end]`
`5 // 2` | `5 ÷ 2`
`2**3` | `2^3`
`1 + 5j` | `1 + 5im`

## Basic operations in Julia

### Functions


=== ":simple-julia: Julia"
    ```julia
    function my_sum(x, y)
        return x + y
    end
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    def my_sum(x, y):
        return x + y
    ```

or in one line:

=== ":simple-julia: Julia"
    ```julia
    my_sum(x, y) = x + y
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    my_sum = lambda x, y: x + y
    ```


### Control flow

Loops:

=== ":simple-julia: Julia"
    ```julia
    for i in 1:5
        println(i)
    end
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    for i in range(1,6):
        print(i)
    ```

`if`, `then`, `else`

=== ":simple-julia: Julia"
    ```julia
    if x < y
        println("x is less than y")
    elseif x > y
        println("x is greater than y")
    else
        println("x is equal to y")
    end
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    if x < y:
        print('x is less than y')
    elif x > y:
        print('x is greater than y')
    else:
        print('x is equal to y')
    ```

### Arrays

Creating arrays. For Python, assume we `import numpy as np`.

=== ":simple-julia: Julia"
    ```julia
    a = [1., 5, 6, 10]
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    a = np.array([1., 5, 6, 10])
    ```

Uninitiated arrays:

=== ":simple-julia: Julia"
    ```julia
    a = Array{Float32}(undef, 100, 100)
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    a = np.empty((100, 100), np.float32)
    ```

Slicing and indexing

=== ":simple-julia: Julia"
    ```julia
    a[1]
    a[end]
    a[end-1]
    # If a is a 2D array:
    a[5:9, 10:15]
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    a[0]
    a[-1]
    a[-2]
    # If a is a 2D array:
    a[4:9, 9:15]
    ```

Element-wise multiplication of arrays:


=== ":simple-julia: Julia"
    ```julia
    a .* b
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    a * b
    ```

Matrix multiplication:


=== ":simple-julia: Julia"
    ```julia
    a * b
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    a @ b
    ```

### Interpolation

Linear interpolation:

=== ":simple-julia: Julia"
    ```julia
    using Interpolations

    x = range(0, 4π, 50)
    y = cos.(x)
    interp_linear = linear_interpolation(x, y)

    new_x = [2, 4, 6, 8, 10, 12]
    new_y = interp_linear.(new_x)
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    import numpy as np

    x = np.linspace(0, 4*np.pi, 50)
    y = np.cos(x)

    new_x = [2, 4, 6, 8, 10, 12]
    new_y = np.interp(new_x, x, y)
    ```

### Plots

For Julia, using [Makie](https://makie.org/). For Python, using [Matplotlib](https://matplotlib.org/).


Simple plots:

=== ":simple-julia: Julia"
    In the REPL, use `GLMakie` for interactive plots. In the notebook, `CairoMakie` for static plots or `WGLMakie` for interactive plots (may not work in Jupyter).
    ```julia
    using CairoMakie

    x = range(0, 4π, 100)
    y = cos.(x)

    fig = Figure(size=(500,300))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Signal (mV)", title="Some data")
    lines!(ax, x, y, color=:red)
    scatter!(ax, x[1:10:end], y[1:10:end], color=:blue)
    fig
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    %matplotlib inline

    import matplotlib.pylab as plt
    import numpy as np

    x = np.linspace(0, 4*np.pi, 100)
    y = np.cos(x)

    fig, ax = plt.subplots(figsize=(5, 3), dpi=100)

    ax.plot(x, y, 'r-')
    ax.plot(x[::10], y[::10], 'bo')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Signal (mV)')
    ax.set_title('Some data')
    ```

Plotting images:


=== ":simple-julia: Julia"
    ```julia
    xs = range(0, 2π, length=100)
    zs = [sin(x * y) for x in xs, y in xs]

    fig, ax, hm = heatmap(xs, ys, zs; colormap=:inferno)
    Colorbar(fig[:, end+1], hm)
    fig
    ```

=== ":fontawesome-brands-python: Python"
    ```python
    tmp = np.linspace(0, 2*np.pi, 100)
    xs, ys = np.meshgrid(tmp, tmp)
    zs = np.sin(xs * ys)

    c = plt.imshow(zs, cmap='inferno', 
                   extent=(0, 2*np.pi, 0, 2*np.pi))
    plt.colorbar(c)
    ```



