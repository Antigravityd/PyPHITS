
# Table of Contents

1.  [Motivation](#org6bb0ea0)
2.  [Example](#orgdfe150a)
3.  [Result](#orgc1e05a4)
4.  [Code](#org93668d2)
5.  [More Code](#orga6c9b41)
6.  [TODOs](#org66607b0)
7.  [Assurances](#org67cdef9)
8.  [Installation](#org449a16f)
9.  [Docs-For-Docs](#orgfa81bae)
    1.  [In-REPL Documentation Features](#org8768f79)


<a id="org6bb0ea0"></a>

# Motivation

Many computational physics tools present Python interfaces these days to exploit the extensive data science tooling available in the language. Also, almost everyone knows it.

PHITS is still based on input cards, which make it both harder to use and less flexible: the card format emphasizes brevity over interpretability, and provides only a limited number of primitive options that require external tools to use them in moderately complex ways.

For example, if one wants to set the horizontal axis of a graph to anything other than linear or logarithmic, it's necessary to generate the points using an external tool and paste them into the input.

Specifying the input via an object system in a general-purpose programming language is both more readable, as objects can be named instead of numbered and "ownership" relationships can be explicit (e.g. a cell "has a" surface), and more expressive, e.g. using comprehension expressions to generate the points in the example below.


<a id="orgdfe150a"></a>

# Example

    
     # Construct objects via list comprehensions
    mats = [Material([("C", i/100), ("H", (2/3)*(100-i)/100), ("O", (1/3)*(100-i)/100)]) for i in range(55,100)]
    
    # Or loops
    cells = [Cell([Sphere(1)], mats[0], -1)]
    for i in range(1,11):
        new = Cell([Sphere(i+1), Sphere(i, inside=False)], mats[i], -1)
        cells.append(new)
    
    ocean = Void([Sphere(11, inside=False)])
    
    escape_e = DumpFluence(cells[-1], ocean, 4*pi/3, [1, 8], "flux")
    
    
    cells.append(ocean)
    
    source = Cylindrical("241Am", 2.2, fissile="neutrons", bounds=(-0.25,0.25), r_out=0.3)
    
    
    # Capture input for further analysis in Python as your choice of many common data formats,
    # or render and return the .eps
    # Pass raw input data in case of an error or unimplemented feature
    inp = make_input(cells, source, escape_e, raw="$ This could be e.g. an [Elastic Option] section text.\n",
    		 parameters={"negs": 1, "e-mode": 2})
    print(inp)


<a id="orgc1e05a4"></a>

# Result

The example produces the following result:

    [Title]
     2023-02-19 20:	09: 47.397746
    [Parameters]
       negs = 1
     e-mode = 2
    [Source]
     s-type = 1
       proj = 241Am
      ispfs = 2
         z0 = -0.25
         z1 = 0.25
         r0 = 0.3
         e0 = 2.2
    [Material]
    MAT[1]
     C 0.6 H 0.26666666666666666 O 0.13333333333333333
    
    MAT[2]
     C 0.63 H 0.24666666666666665 O 0.12333333333333332
    
    MAT[3]
     C 0.65 H 0.2333333333333333 O 0.11666666666666665
    
    MAT[4]
     C 0.55 H 0.3 O 0.15
    
    MAT[5]
     C 0.58 H 0.28 O 0.14
    
    MAT[6]
     C 0.61 H 0.26 O 0.13
    
    MAT[7]
     C 0.56 H 0.29333333333333333 O 0.14666666666666667
    
    MAT[8]
     C 0.57 H 0.2866666666666666 O 0.1433333333333333
    
    MAT[9]
     C 0.64 H 0.24 O 0.12
    
    MAT[10]
     C 0.59 H 0.2733333333333333 O 0.13666666666666666
    
    MAT[11]
     C 0.62 H 0.2533333333333333 O 0.12666666666666665
    
    [Surface]
     1  SO 7
     2  SO 10
     3  SO 8
     4  SO 4
     5  SO 1
     6  SO 6
     7  SO 2
     8  SO 11
     9  SO 3
     10 SO 9
     11 SO 5
    [Cell]
     1 7 -1 7 5
    
    
     2 3 -1 8 2
    
    
     3 1 -1 6 11
    
    
     4 8 -1 9 7
    
    
     5 2 -1 10 3
    
    
     6 -1 8
    
    
     7 9 -1 2 10
    
    
     8 11 -1 3 1
    
    
     9 10 -1 11 4
    
    
     10 6 -1 1 6
    
    
     11 4 -1 5
    
    
     12 5 -1 4 9
    
    
    [T-Cross]
     mesh = reg
     unit = 1
     axis = reg
     file = cross.dmp
    output_type = flux
      reg = 1
     mesh = reg
     unit = 1
     axis = reg
     file = cross.dmp
    output_type = flux
      reg = 1
     r-out r-in area
     2     6    4.1887902047863905
     dump = -2
     1 8
     $ This could be e.g. an [Elastic Option] section text.


<a id="org93668d2"></a>

# Code

The architecture of the code is relatively simple.

There's a factory base class, `PhitsObject`, that provides the common initialization, hashing, equality, and text definition methods for consequential objects.

This allows simple definitions of consequential objects, as follows:

    
    class Material(PhitsObject): # Composition is a list of pairs of (<element name string>, <ratio>) e.g. ("8Li", 0.5)
        name = "material"
        required = ["composition"]
        positional = ["composition"]
        optional = ["time_change", "data_max", "mat_name_color", "condensed", "conductive", "electron_step",
    		"neutron_lib", "proton_lib", "electron_lib", "photon_lib", "thermal_lib"]
        shape = (lambda self: f"MAT[{self.index}]",
    	     (lambda self: "".join(map(lambda tup: f"{tup[0]} {tup[1]} ", self.composition))), "condensed",
    	     "conductive", "electron_step", "neutron_lib", "proton_lib","electron_lib", "photon_lib",
    	     lambda self: f"MT{self.index} {self.thermal_lib}" if self.thermal_lib is not None else "")
        subobjects = ["time_change", "data_max", "mat_name_color"]
        ident_map = {"condensed": "GAS", "conductive": "COND", "electron_step": "ESTEP", "neutron_lib": "NLIB",
    		 "proton_lib": "HLIB", "electron_lib": "HLIB", "photon_lib": "PLIB"}

The `name` attribute indicates which section the object belongs to, `required` indicates the arguments to the initialization that are required, `positional` indicates those arguments which must be specified positionally, `optional` indicates those additional arguments which are optional and must be specified by keyword, `shape` indicates how the attributes are put into the input file, `subobjects` indicates those arguments which, if given, will be another `PhitsObject`, `ident_map` associates some of the object's attribute names with an alternative name so that `<alternative name> = <value>`  instead of `<attribute name> = <value>` is placed in the input, `value_map` does the same but for the `<value>` part of the equality, and `nones` maps attribute names to the value to place in the file when they're omitted instead of omitting them from the file altogether (e.g. the `non` skip operator).

The shape attribute is a tuple of strings, functions, and more tuples. If an element is a string, it's placed verbatim into the file if it's not an attribute name of the object (with a preceding backslash escaping this behavior, in the case of a collision), and otherwise `<attribute name> = <attribute value>\n` is placed in the file. If it's a function, that function is called with the object in question as a sole argument and the result of the function is placed into the file verbatim with a trailing newline. If it's a tuple, the same rules apply, but the `<attribute name> =` part is omitted and a trailing space instead of a trailing newline is added (for the function objects too). If any string ends with a backslash, the trailing newline is instead a space. Additionally, the two spaces on either side of the equals sign are removed, if an equals sign were otherwise to appear in the input file.


<a id="orga6c9b41"></a>

# More Code

The end user constructs cells, sources, and any non-cellwise tallies themselves, and passes them to `run_phits()`.  The work of constructing the input file is done in `make_input()`, which first puts all the passed objects in a Python `set` to ensure uniqueness (according to the definition of equality defined on `PhitsObject`), and then constructs the `type_divided` dictionary of lists that separates the objects by the section they are to end up in, also fixing the order in which they'll appear.

The objects are then given an index corresponding to their index in the list in `type_divided`, and replaces objects in the class hierarchy with the object in `type_divided` to which they are equal (since someone could initialize two `Surface` objects with identical arguments that are distinct in the eyes of Python).

It's then a simple matter of constructing the string representing the input file via each object's `definition()` method, which is done with the help of an `add_defs()` function to eliminate boilerplate.

`run_phits()` is quite simple: it creates a temporary folder, changes to it, creates the `.inp` file from the string returned by `make_input()`, and uses the `subprocess` library to run a shell command on that file to execute PHITS (by default, the `phits` shell script I wrote to make the calling syntax POSIX-compliant). A function (incomplete) that scrapes the result of the computation is then run, and that result is returned (ideally in the format of the user's choice: a `dict`, a Pandas dataframe, a Numpy array, or a `matplotlib` figure representing the same plot that PHITS would have constructed).


<a id="org66607b0"></a>

# TODOs

I originally intended to fully support all of the tallies available in PHITS. However, this would require developing a full-featured parser for the AnGeL files so-generated, as in order to support the Python integration, the contained data would need to be extracted. However, this proved to be rather involved; the decaying remains of this parser are still in the `angel` subdirectory, should it ever be useful in the future. Instead of doing this, I moved to allowing only the three tallies for which one can specify the `dump` option, as this produces a sensible, easily machine-readable file containing exactly the desired data. Many of the other tallies should be emulatable in post from these tallies alone.

The core object model currently doesn't support a few edge cases:

-   Non-tetrahedral `LAT` (just use a Python loop) and  `U` / `FILL` (stumped) options in `[Cell]`
-   Any tally besides `[T-Track]`, `[T-Product]`, or `[T-Time]` with the dump options (it'd require parsing that's too complicated; most of the functionality may be recovered *in post* anyway.)
-   Some of the grouped row-like data don't support optional arguments (e.g. `non` in the input) because PHITS requires that at least one of the rows in any given column be non-skipped, which is harder to implement checks of than it's worth.

PHITS's line length limit constrains machine-generated inputs' complexity. Writing a patch for PHITS's source that puts line-length in `param.inc` would greatly enhance the usability of this tool.

Getting semantic explanations from the PHITS documentation somehow would be nice.

The tests are *extremely* slow. I don't know if this is because I run them on a stone-age Thinkpad, if there are genuine optimizations to be made, or (most likely) some combination of those.

It would be nice if the Fortran-Python interop would make it possible to write the user-defined sources, cross-sections, functions, and tallies in Python. Similarly, defining the (electro)magnetic field files and other external files via Python would be feasible.


<a id="org67cdef9"></a>

# Assurances

The code is property-based tested via Hypothesis. Specifically, there is high confidence that any PHITS input accessible via these functions for which `icntl = 3` finds an error (more precisely: for which PHITS's output includes some case-altered version of "error") raises a more-descriptive Python `ValueError` before it is created.

Caveats:

-   The testing is random, so low- and zero-probability events may still occur (no assurance short of a formal proof mitigates this)
-   There is no guarantee that the program is free of semantics-altering errors—i.e. if a bug substitutes a correct, intended input file for a correct, but unintended input file.
-   There is no guarantee that the descriptive Python errors so-raised aren't too agressive, and prevent the user from accessing some valid PHITS input
-   Certain value-interdependence-based restrictions on `Parameters` objects aren't dealt with (and probably would catastrophically degrade test performance if they were).
-   Objects which require files to be present at runtime cannot be completely tested this way, as generating such files is more trouble than it's worth; there are checks which precede the check for the file's existence, and these can still be tested (one merely must verify that the PHITS error one eventually encounters is of a `FileNotFound` nature). The objects in question are `Tetrahedral`, `MappedMagneticField`, `MappedElectromagneticField`,
    `FragData`, and `TetrahedralSource`.
-   The test implementation is more unit than integration—the types are iterated through, generated, and minimal supporting infrastructure for each instance is used to produce an input that tests mostly the type itself. As a result, the implementation may not be robust to errors that arise due to combinations of sections. This is somewhat mitigated by the `Cell`-adjacent objects, whose subobjects are all (save the offenders above) generated and tested in this process.
-   There may be inputs which are over the PHITS line length limit but are otherwise in error. I expect these to be statistically independent, as I'm necessarily ignorant about the nature of "error," but they may well not be.


<a id="org449a16f"></a>

# Installation

TODO pending pip


<a id="orgfa81bae"></a>

# Docs-For-Docs

For documentation-by-example, see the `ex` folder in the source tree. A detailed API reference is available [here](https://antigravityd.github.io/PyPHITS). Note that this reference is in large part automatically generated, and so may have some quirky explanations at this stage. 


<a id="org8768f79"></a>

## In-REPL Documentation Features

Python's `help()` function will produce the top-level docstring for any object (not just this package's). Calling `syntax_desc()` on any subclass of `PhitsObject` or on `Parameters` (note: classmethod) will produce coarse-grained information about the arguments that function accepts. Look to the PHITS manual if this isn't sufficient. Similarly, `syntax_for` will produce the subset of this information corresponding to a particular (keyword) argument. If you know the PHITS identifier, but not the Python identifier, pass `phits=True`.

