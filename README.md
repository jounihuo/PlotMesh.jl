# PlotMesh.jl

[![Build Status](https://travis-ci.com/jounihuo/PlotMesh.jl.svg?branch=master)](https://travis-ci.com/jounihuo/PlotMesh.jl)
[![Coveralls](https://coveralls.io/repos/github/jounihuo/PlotMesh.jl/badge.svg?branch=master)](https://coveralls.io/github/jounihuo/PlotMesh.jl?branch=master)

Simple mesh plot for JuliaFEM h5-files.

```julia
julia> using PlotMesh

julia> plotmesh("test/test_geom.h5", outputfn="test_geom.png")
width:    400.0
height:   400.0
filename: test_geom.png
type:     png

```

![mesh image](https://raw.githubusercontent.com/jounihuo/PlotMesh.jl/master/test/test_geom.png)
