function gen_mzml(files::AbstractVector{<: AbstractString}; peak_picking_rule::Union{AbstractString, Nothing} = "1-", custom_args::AbstractVector{<: AbstractString} = String[], force::Bool = false)
    tool = "msconvert"
    args = ["--mzML", "--64", "--zlib", "--simAsSpectra"]

    if !isa(peak_picking_rule, Nothing)
        push!(args, "--filter \"peakPicking true $peak_picking_rule\"")
    end
    args = [args; custom_args]

    jobs = Channel(ctype = Tuple{Int, Cmd}, csize = length(files)) do channel
        for (i, file) in enumerate(files)
            file_split = splitext(file)
            if file_split[2] != ".raw"
                file = file_split[1] * file_split[2] * ".raw"
                file_mzml = file_split[1] * file_split[2] * ".mzML"
            else
                file_mzml = file_split[1] * ".mzML"
            end

            if !isfile(file)
                @warn "[$i] File $file is missing and will be skipped."
                continue
            end
            if isfile(file_mzml)
                if force
                    println("[$i] File $file_mzml is already present and will be overwritten.")
                else
                    println("[$i] File $file_mzml is already present and will be skipped.")
                    continue
                end
            end

            cmd = "$tool $file $(join(args, " "))"
            put!(channel, (i, `powershell $cmd`))
        end
    end

    function run_command()
        for (job_id, cmd) in jobs
            println("[$job_id] $cmd")
            try
                open(cmd) do process
                    for output in eachline(process)
                        println("[$job_id] $output")
                    end
                end
            catch error
                @warn "[$job_id] $(error.msg)"
            end
        end
    end

    n = Sys.CPU_THREADS
    task_list = map(_ -> schedule(Task(run_command)), 1:(n == 1 ? 1 : (n - 1)))

    while(!all(istaskdone.(task_list)))
        yield()
    end

    return
end

gen_mzml(file::AbstractString; peak_picking_rule::Union{AbstractString, Nothing} = "1-", custom_args::AbstractVector{<: AbstractString} = String[], force::Bool = false) = gen_mzml([file], peak_picking_rule = peak_picking_rule, custom_args = custom_args, force = force)
gen_mzml(files::AbstractVector{<: Pair{T, String}}; peak_picking_rule::Union{AbstractString, Nothing} = "1-", custom_args::AbstractVector{<: AbstractString} = String[], force::Bool = false) where {T} = gen_mzml(last.(files), peak_picking_rule = peak_picking_rule, custom_args = custom_args, force = force)

function msfiles(files::Vararg{AbstractString})
    [Symbol(file) => file for file in files]
end

function msfiles(; files...)
    collect(files)::Vector{Pair{Symbol, String}}
end
