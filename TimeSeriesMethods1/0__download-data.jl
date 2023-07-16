using Distributed, ClusterManagers

if "SLURM_JOBID" ∈ keys(ENV)
    @info "Working on a slurm cluster"
    addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"])-1)
else
    @info "Working locally"
    addprocs(Threads.nthreads())
end


# using Dagger


# add relevant packages
@everywhere begin
    using Pkg
    println("Activating Environment")
    Pkg.activate(".")
    println("Instantiating")
    Pkg.instantiate()
    println("Precompiling")
    Pkg.precompile()

    println("Importing Packages")
    using Distributed
    using CSV
    using ProgressMeter

    # set up data directory
    outpath = "data/sharedair/raw"
end



if !ispath(outpath)
    mkpath(outpath)
end


@everywhere begin
    # load in functionality for AWSS3
    include("utils/osn_anonymous.jl")

    endpoint = "https://ncsa.osn.xsede.org"
    bucket = "s3://ees230012-bucket01"

    # configure AWS.jl to use our OSN configuration
    AWS.global_aws_config(AnonymousOSN(endpoint))

    nodes_to_use = ["Central_Hub_1", "Central_Hub_2", "Central_Hub_4"]
    years_to_use = ["2022", "2023"]
    months = [lpad(i, 2, "0") for i ∈ 1:12]
    days = [lpad(i, 2, "0") for i ∈ 1:31]
    sensors_to_use = ["IPS7100", "BME680"]
end


for node ∈ nodes_to_use
  for sensor ∈ sensors_to_use
      if !ispath(joinpath(outpath, node, sensor))
          mkpath(joinpath(outpath, node, sensor))
      end
  end
end



@everywhere begin
    function get_paths(p, sensors_to_use)
        df_paths = []

        if isdir(p)
            for f ∈ readdir(p)
                if any(occursin.(sensors_to_use, f))
                    idx = findfirst(occursin.(sensors_to_use, f))
                    push!(df_paths, (joinpath(p, f), sensors_to_use[idx]))
                end
            end
        end

        return df_paths
    end


    function process_paths(ptuples, outpath)
        for ptuple ∈ ptuples
            path, sensor = ptuple

            splitpath = split(convert(String, path), "/")
            node = splitpath[7]
            fname = splitpath[end]

            # download the file to our desired output directory
            # csv = Dagger.@spawn CSV.File(path)
            # Dagger.@spawn CSV.write(joinpath(outpath, node, sensor, fname), csv)
            csv = CSV.File(path)
            CSV.write(joinpath(outpath, node, sensor, fname), csv)
        end
    end
end


# @showprogress pmap(1:6) do x
#     sleep(5)
#     x^2
# end


# pmap([(i, j) for i ∈ 1:3 for j ∈ 1:5]) do (i,j)
#     println(i*j)
# end




# loop over node, year, month, day and
# 1. fetch the files matching sensor_list
# 2. download and place in data directory


flist = [(node, year, month, day) for node ∈ nodes_to_use for year ∈ years_to_use for month ∈ months for day ∈ days]

@showprogress pmap(flist) do (node, year, month, day)
    p = S3Path(joinpath(bucket, "AirQualityNetwork/data/raw", node, year, month, day*"/"))

    # ptuples = Dagger.@spawn get_paths(p, sensors_to_use)
    # Dagger.@spawn process_paths(ptuples, outpath)
    ptuples = get_paths(p, sensors_to_use)
    process_paths(ptuples, outpath)
end


# # task 1: for each month, walk directories and generate file lists. Afterwards, fetch and join results into single file list

# add1(value) = value + 1
# add2(value) = value + 2
# combine(a...) = sum(a)

# res1 = [Dagger.@spawn add1(i) for i ∈ 1:10]
# res = Dagger.@spawn combine(res1...)
# fetch(res)

# [i*j for i ∈ 1:3 for j∈1:5]

# [S3Path(joinpath(bucket, "AirQualityNetwork/data/raw", node, year, month, day*"/")) for node ∈]
