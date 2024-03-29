#r "nuget: MathNet.Numerics.FSharp, 5.0.0"
#r "nuget: Plotly.NET"
open MathNet.Numerics.Random
open MathNet.Numerics
open Plotly.NET

let rng = Random.shared
let buildLattice n = Array2D.init n n (fun _ _ ->  2 * rng.Next(2) - 1)

let magnetization lattice = 
    let length = 
        lattice |> Array2D.length1
    let r = 
        Array.init length (fun i -> lattice.[i,*] |> Array.sum)
    let mag:int = r |> Array.sum
    mag

let Energy lattice : int= 
    let mutable energy = 0
    let eVec vec = vec |> Array.pairwise|>Array.map (fun (a,b) -> a*b)|> Array.sum
    for n = 0 to (lattice |> Array2D.length1)-1 do
        energy <- energy + (eVec lattice.[n,*]) + (eVec lattice.[*,n])
    (-1) * energy

/// Slower with these functions
(*let neighbors lattice location = 
    let a,b = location
    let length = lattice |> Array2D.length1
    let edge = length - 1
    let horizontal =
        match b with
        | b when b=0 -> lattice.[a,b+1]
        | b when b=edge -> lattice.[a,b-1]
        | b when b<edge && b>0 -> (lattice.[a,b-1] + lattice.[a,b+1])
        |_ -> 0
    let vertical =
        match a with
        | a when a=0 -> lattice.[a+1,b]
        | a when a=edge -> lattice.[a-1,b]
        | a when a>0 && a<edge -> lattice.[a-1,b] + lattice.[a+1,b]
        |_ -> 0
    horizontal + vertical
let flipDE lattice location :int =
    let a,b = location
    let pairs = neighbors lattice location
    let initE = lattice.[a,b] * pairs
    let flipE = -1 * initE
    (-1) * (flipE - initE)*)
///
let pointEnergy lattice location :int =
    let (a,b) = location
    let length = (lattice|>Array2D.length1) - 1
    let pairs loc =
        match loc with
        | a,b when a=0 && b=0 -> lattice.[a,b] * (lattice.[1,0] + lattice.[0,1])
        | a,b when a=0 && b=length -> lattice.[a,b] * (lattice.[0,b-1] + lattice.[1,b])
        | a,b when a=length && b=0 -> lattice.[a,b] * (lattice.[a-1,0] + lattice.[a,1])
        | a,b when a=length && b=length -> lattice.[a,b] * (lattice.[a,b-1] + lattice.[a-1,b])
        | a,b when a = 0 && b>0 && b<length -> lattice.[a,b] * (lattice.[a,b+1] + lattice.[a,b-1] + lattice.[a+1,b])
        | a,b when a = length && b>0 && b<length -> lattice.[a,b] * (lattice.[a,b+1] + lattice.[a,b-1] + lattice.[a-1,b])
        | a,b when b = 0 && a>0 && a<length -> lattice.[a,b] * (lattice.[a,b+1] + lattice.[a-1,b] + lattice.[a+1,b])
        | a,b when b = length && a>0 && a<length -> lattice.[a,b] * (lattice.[a-1,b] + lattice.[a+1,b] + lattice.[a,b-1])
        | a,b when a>0 && a<length && b>0 && b<length ->
            lattice.[a,b] * (lattice.[a-1,b] + lattice.[a+1,b] + lattice.[a,b-1] + lattice.[a,b+1])
        | _,_ -> 0
    let initial = pairs location
    let flip = -1 * initial
    (-1) * (flip - initial)
    //2 * initial

let mcmc lat (beta:float) (nflips:int) (tsteps:int) =
    let lattice = lat |> Array2D.copy
    let length = lattice |> Array2D.length1
    let mutable energy:int = Energy lattice
    let av_energy = Array.zeroCreate tsteps
    let av_mag = Array.zeroCreate tsteps
    for step=0 to tsteps-1 do
        for flip=0 to nflips-1 do
            let flip_i,flip_j = (rng.Next(length), rng.Next(length))
            let flip_dE = pointEnergy lattice (flip_i,flip_j)
            if exp (-beta * (float flip_dE) ) > rng.NextDouble() then
                energy <- energy + flip_dE
                lattice.[flip_i,flip_j] <- -1*lattice.[flip_i,flip_j]
        
        av_energy[step] <- energy
        av_mag[step] <- lattice |> magnetization
    (av_energy,av_mag)

let lat20 = buildLattice 20

let timer = new System.Diagnostics.Stopwatch()
timer.Start()
let e1,m1 = mcmc lat20 0.9 100000 20
printfn "MCMC took %i milliseconds" timer.ElapsedMilliseconds

let mc_tlist = [1..1..20]

let ePlot = 
    Chart.Point(x=mc_tlist , y=e1)
    |> Chart.withYAxisStyle (TitleText = "Energy ", MinMax=(-820, 0))
    |> Chart.withXAxisStyle (TitleText = "Sample step")
    |> Chart.show
    

let mPlot = 
    Chart.Point(x=mc_tlist , y=m1)
    |> Chart.withYAxisStyle (TitleText = "Magnetization ")
    |> Chart.withXAxisStyle (TitleText = "Sample step")
    |> Chart.show