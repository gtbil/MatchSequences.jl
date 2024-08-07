# MatchSequences.jl

Grant Billings 5/31/23

[![Build Status](https://github.com/gtbil/MatchSequences.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/gtbil/MatchSequences.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Create the package

```julia
t = Template(; user="gtbil", authors=["Grant Billings"], plugins=[License(name="MIT"), Git(), GitHubActions(),],)
t("MatchSequences")
```

This create a folder in `C:\Users\grant\.julia\dev\MatchSequences`. For getting git working, delete the default `.git` folder if it exists and copy in the one you get from the following (after authenticating):

```bash
gh repo clone gtbil/MatchSequences.jl
# then copy files into the same directory as the project above^
```

To load the environment in the Julia REPL in VSCode

```julia
import Pkg
Pkg.activate("MatchSequences")
# then to close
Pkg.activate()
```

To memory profile

```julia
using Profile
using PProf

Profile.Allocs.clear()

# Profile.Allocs.@profile sample_rate=0.01 main("../CokerTest/Coker312_A01.fasta")
Profile.Allocs.@profile sample_rate=1 MatchSequences.benchmark(t′, sa, l_vec, lms_vec, out_vec, heads, tails,
                               lms_subvec1, lms_subvec2, 
                               lms_boolvec1, lms_boolvec2, numlms; σ_in = col);

PProf.Allocs.pprof(from_c = false)
```

Enable the debugger:
```julia
ENV["JULIA_DEBUG"] = MatchSequences
```