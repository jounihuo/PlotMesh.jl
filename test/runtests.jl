using PlotMesh
using Test
using VisualRegressionTests

istravis = "TRAVIS" ∈ keys(ENV)


function testfun(fname)
    plotmesh("test_geom.h5";outputfn=fname)
end

@testset "PlotMesh.jl" begin
     #@test test_images(VisualTest(testfun, "test_geom.png"), popup=!istravis) |> save_comparison |> success
     @visualtest testfun "test_geom.png" false #!istravis
end
