#r "nuget: MathNet.Numerics.FSharp, 5.0.0"
#r "nuget: Plotly.NET"
open MathNet.Numerics.Random
open MathNet.Numerics
open Plotly.NET

let rng = Random.shared
let buildLattice n = Array2D.init n n (fun _ _ ->  float (2 * rng.Next(2) - 1))
let coldStart n = Array2D.create n n 1.0

let magnetization (lattice: float array2d) :float = 
    let length = 
         (lattice |> Array2D.length1)
    let r = 
        Seq.init length (fun i -> lattice.[i,*] |> Array.sum)
    abs(float (r |> Seq.sum) / ( (float length)**2) )

let Energy lattice : float= 
    let mutable energy  = 0.0
    let length = lattice |> Array2D.length1
    let eVec vec = vec |> Seq.pairwise|>Seq.map (fun (a,b) -> a*b)|> Seq.sum
    for n = 0 to length-1 do
        energy <- energy + (eVec lattice.[n,*]) + (eVec lattice.[*,n])
    (-1.0) * energy/900.0

let pointEnergy lattice location :float =
    let (a,b) = location
    let length = (lattice|>Array2D.length1) - 1
    let pairs loc :float =
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
    let flip = -1.0 * initial
    (-1.0) * (flip - initial)
    //2 * initial

let mcmc (lat: float array2d) (beta:float) (nflips:int) (tsteps:int) : float array * float array =
    let lattice: float array2d = lat |> Array2D.copy
    let length = lattice |> Array2D.length1
    let mutable energy:float = Energy lattice
    let av_energy = Array.zeroCreate tsteps
    let av_mag:float array = Array.zeroCreate tsteps
    for step=0 to tsteps-1 do
        for flip=0 to nflips-1 do
            let flip_i,flip_j = (rng.Next(length), rng.Next(length))
            let flip_dE = float (pointEnergy lattice (flip_i,flip_j))
            if exp (-beta *  flip_dE ) > rng.NextDouble() then
                energy <- energy + flip_dE
                lattice.[flip_i,flip_j] <- -1.0*lattice.[flip_i,flip_j]
        
        av_energy[step] <- energy/(float length)**2
        av_mag[step] <- (lattice |> magnetization)
    (av_energy,av_mag)

let lat20 = buildLattice 20

let timer = new System.Diagnostics.Stopwatch()
timer.Start()
//let e1,m1 = mcmc lat20 0.6 100000 20
let e1,m1 = mcmc (30|>buildLattice) 0.7 (int (10.0**4)) 1000
//let e1,m1 = mcmc (30|>coldStart) 0.6 100000 200
printfn "MCMC took %i milliseconds" timer.ElapsedMilliseconds

let mc_tlist = [1..1..1000]

let ePlot = 
    Chart.Point(x=mc_tlist , y=e1)
    |> Chart.withYAxisStyle (TitleText = "Energy ", MinMax=(-2.2, 0))
    |> Chart.withXAxisStyle (TitleText = "Sample step")
    |> Chart.show

let mPlot = 
    Chart.Point(x=mc_tlist , y=m1)
    |> Chart.withYAxisStyle (TitleText = "Magnetization ", MinMax=(-1.1,1))
    |> Chart.withXAxisStyle (TitleText = "Sample step")
    |> Chart.show

let mHist =
    Chart.Histogram(Y=m1, NBinsY = 100)
    |> Chart.withYAxisStyle (TitleText = "Magnetization ", MinMax=(0,1))
    |> Chart.withXAxisStyle (TitleText = "Counts")
    |> Chart.show

let EvsT, MvsT, ChivsT = 
    //let lattice = 30|> buildLattice // hot start
    let lattice = 30 |> coldStart
    let bList = [0.2..0.01..0.6]
    let length = bList |> List.length
    let eList = Array.zeroCreate length
    let mList = Array.zeroCreate length
    let chiList = Array.zeroCreate length
    for i=0 to (length-1) do        
        let e1,m1 = (mcmc(lattice) bList.[i] 10000 1000)
        eList.[i] <- e1|>Array.average
        mList.[i] <- m1|>Array.average  
        chiList.[i] <- ((m1|>Array.map (fun mag -> mag**2)|> Array.average) - (mList.[i])**2)
                          
    eList,mList,chiList    

let etPlot = 
    Chart.Point(x=[0.2..0.01..0.6],y=EvsT)
    |> Chart.withYAxisStyle (TitleText = "Energy ", MinMax=(-2.2, 0))
    |> Chart.withXAxisStyle (TitleText = "Beta")
    |> Chart.show
let mtPlot = 
    Chart.Point(x=[0.2..0.01..0.6],y=MvsT)
    |> Chart.withYAxisStyle (TitleText = "Magnetization ", MinMax=(-0.1,1.05))
    |> Chart.withXAxisStyle (TitleText = "Beta")
    |> Chart.show
let chiTplot = 
    Chart.Point(x=[0.2..0.01..0.6], y=ChivsT)
    |> Chart.withYAxisStyle (TitleText = "Chi ")
    |> Chart.withXAxisStyle (TitleText = "Beta")
    |> Chart.show

let betaCrit = log(1.0 + sqrt(2.0))/2.0
